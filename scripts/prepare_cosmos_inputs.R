library(readr)
library(decoupleR)

metabolomic_clean_vsn_common <- as.data.frame(read_csv("data/metabolomic/metabolomic_clean_vsn_common.csv"))
metabolite_matching <- as.data.frame(read_csv("support/metabolite_matching.csv"))
metabolomic_clean_vsn_common <- merge(metabolite_matching,metabolomic_clean_vsn_common, by = "metabolite")
metabolomic_clean_vsn_common <- metabolomic_clean_vsn_common[,-1]


SDs <- apply(metabolomic_clean_vsn_common[,-1],1,function(x){sd(x,na.rm = T)})
means <- rowMeans(metabolomic_clean_vsn_common[,-1], na.rm = T)

metabs_df_scaled <- metabolomic_clean_vsn_common
metabs_df_scaled[,-1] <- (metabs_df_scaled[,-1] - means) / SDs

RNA_log2_FPKM_cleaned_common <- as.data.frame(read_csv("data/RNA/RNA_log2_FPKM_cleaned_common.csv"))
row.names(RNA_log2_FPKM_cleaned_common) <- RNA_log2_FPKM_cleaned_common$gene
RNA_log2_FPKM_cleaned_common <- RNA_log2_FPKM_cleaned_common[,-1]

doro_regulons <- get_dorothea(levels = c("A","B"))

TF_activities <- apply(RNA_log2_FPKM_cleaned_common,2,function(x){
                        x <- as.data.frame(x[which(!is.na(x))])
                        TFs <- run_wmean(x, network = doro_regulons, times = 1000, minsize = 20)
                        TFs <- as.data.frame(TFs)
                        TFs <- TFs[which(TFs$statistic == "norm_wmean"),c(2,4)]
                        as_input <- TFs[,2]
                        names(as_input) <- TFs[,1]
                        return(as_input)
                        })

row.names(metabs_df_scaled) <- metabs_df_scaled$hmdb
metabs_df_scaled<- metabs_df_scaled[,-1]

metabs_scaled <- apply(metabs_df_scaled, 2, function(x){
  x <- x[which(!is.na(x))]
  return(x)
},simplify = F)

SDs <- apply(RNA_log2_FPKM_cleaned_common,1,function(x){sd(x,na.rm = T)})
means <- rowMeans(RNA_log2_FPKM_cleaned_common, na.rm = T)

RNA_scaled <- (RNA_log2_FPKM_cleaned_common - means) / SDs

RNA_scaled <- apply(RNA_scaled, 2, function(x){
  return(x[which(!is.na(x))])
}, simplify = F)

names(RNA_scaled)
names(metabs_scaled)
names(TF_activities)

cosmos_inputs <- sapply(names(TF_activities), function(cell_line, RNA_scaled, metabs_scaled, TF_activities){
  cosmos_input <- list()
  cosmos_input[["RNA"]] <- RNA_scaled[[cell_line]]
  cosmos_input[["metabolomic"]] <- metabs_scaled[[cell_line]]
  cosmos_input[["TF_scores"]] <- TF_activities[[cell_line]]
  return(cosmos_input)
}, RNA_scaled = RNA_scaled, metabs_scaled = metabs_scaled, TF_activities= TF_activities, USE.NAMES = T, simplify = F)

save(cosmos_inputs, file = "data/cosmos/cosmos_inputs.RData")
