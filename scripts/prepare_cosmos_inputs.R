library(readr)
library(decoupleR)

##Prepare metabolomic inputs
metabolomic_clean_vsn_common <- as.data.frame(read_csv("data/metabolomic/metabolomic_clean_vsn_common.csv"))
metabolite_matching <- as.data.frame(read_csv("support/metabolite_matching.csv"))
metabolomic_clean_vsn_common <- merge(metabolite_matching,metabolomic_clean_vsn_common, by = "metabolite")
metabolomic_clean_vsn_common <- metabolomic_clean_vsn_common[,-1]


SDs <- apply(metabolomic_clean_vsn_common[,-1],1,function(x){sd(x,na.rm = T)})
means <- rowMeans(metabolomic_clean_vsn_common[,-1], na.rm = T)

metabs_df_scaled <- metabolomic_clean_vsn_common
metabs_df_scaled[,-1] <- (metabs_df_scaled[,-1] - means) / SDs

row.names(metabs_df_scaled) <- metabs_df_scaled$hmdb
metabs_df_scaled<- metabs_df_scaled[,-1]

metabs_scaled <- apply(metabs_df_scaled, 2, function(x){
  x <- x[which(!is.na(x))]
  return(x)
},simplify = F)


## Prepare RNA and TF inputs
RNA_log2_FPKM_cleaned_common <- as.data.frame(read_csv("data/RNA/RNA_log2_FPKM_cleaned_common.csv"))
row.names(RNA_log2_FPKM_cleaned_common) <- RNA_log2_FPKM_cleaned_common$gene
RNA_log2_FPKM_cleaned_common <- RNA_log2_FPKM_cleaned_common[,-1]
RNA_log2_FPKM_cleaned_common <- RNA_log2_FPKM_cleaned_common[!(rowSums(is.na(RNA_log2_FPKM_cleaned_common)) > 20),]

SDs <- apply(RNA_log2_FPKM_cleaned_common,1,function(x){sd(x,na.rm = T)})
means <- rowMeans(RNA_log2_FPKM_cleaned_common, na.rm = T)

RNA_scaled <- (RNA_log2_FPKM_cleaned_common - means) / SDs



doro_regulons <- get_dorothea(levels = c("A","B"))

## This is were TF activities are estiamted from transcriptomic data, see https://github.com/saezlab/decoupler for more info
TF_activities <- apply(RNA_scaled,2,function(x){
                        x <- as.data.frame(x[which(!is.na(x))])
                        TFs <- run_wmean(as.matrix(x), network = doro_regulons, times = 1000, minsize = 20)
                        TFs <- as.data.frame(TFs)
                        TFs <- TFs[which(TFs$statistic == "norm_wmean"),c(2,4)]
                        as_input <- TFs[,2]
                        names(as_input) <- TFs[,1]
                        return(as_input)
                        })

RNA_scaled <- apply(RNA_scaled, 2, function(x){
  return(x[which(!is.na(x))])
}, simplify = F)

# ## This part is quick and dirty converting TF activities list into DF
# temp_1 <- data.frame(TF_activities[[1]])
# temp_1$TF <- row.names(temp_1)
# 
# temp_2 <- data.frame(TF_activities[[2]])
# temp_2$TF <- row.names(temp_2)
# 
# TF_activities_df <- merge(temp_1,temp_2, all = T)
# for(i in 3:length(TF_activities))
# {
#   temp <- data.frame(TF_activities[[i]])
#   temp$TF <- row.names(temp)
#   TF_activities_df <- merge(TF_activities_df,temp, all = T, by= "TF")
# }
# 
# rm(temp_1)
# rm(temp_2)
# rm(temp)
# 
# names(TF_activities_df)[-1] <- names(TF_activities)
# row.names(TF_activities_df) <- TF_activities_df$TF
# TF_activities_df <- TF_activities_df[,-1]
# 
# ## Now we can scale the TF activities too
# SDs <- apply(TF_activities_df,1,function(x){sd(x,na.rm = T)})
# means <- rowMeans(TF_activities_df, na.rm = T)
# 
# TF_activities_df_scaled <- (TF_activities_df - means) / SDs
# 
# TF_activities_scaled <- apply(TF_activities_df_scaled, 2, function(x){
#   return(x[which(!is.na(x))])
# }, simplify = F)
# 
# names(RNA_scaled)
# names(metabs_scaled)
# names(TF_activities_scaled)

#format the intputs for cosmos and save
cosmos_inputs <- sapply(names(TF_activities), function(cell_line, RNA_scaled, metabs_scaled, TF_activities_scaled){
  cosmos_input <- list()
  cosmos_input[["RNA"]] <- RNA_scaled[[cell_line]]
  cosmos_input[["metabolomic"]] <- metabs_scaled[[cell_line]]
  cosmos_input[["TF_scores"]] <- TF_activities[[cell_line]]
  return(cosmos_input)
}, RNA_scaled = RNA_scaled, metabs_scaled = metabs_scaled, TF_activities= TF_activities, USE.NAMES = T, simplify = F)

save(cosmos_inputs, file = "data/cosmos/cosmos_inputs.RData")
