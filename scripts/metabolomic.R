library(readr)
library(reshape2)
library(pheatmap)
library(vsn)

library(metaboliteIDmapping)

metabs <- as.data.frame(read_csv("data/metabolomic/WEB_DATA_METABOLON.TXT"))
metabs <- metabs[!grepl("Isobar",metabs$TITLE),]
metabs <- metabs[!grepl("^X[-]",metabs$TITLE),]
metabs <- metabs[!grepl("^Possible",metabs$TITLE),]

metabs_df <- dcast(metabs, formula = TITLE~cellname, value.var = "VALUE",fun.aggregate = mean)
row.names(metabs_df) <- metabs_df$TITLE
metabs_df <- metabs_df[,-1]

hist(as.numeric(unlist(log2(metabs_df))), breaks = 1000) #Things looks like they are getting funky around 5

pheatmap(log2(metabs_df), show_colnames = F, show_rownames = F) 

metabs_df[abs(log2(metabs_df)) > 5] <- NA 

fit <- vsnMatrix(as.matrix(metabs_df))
meanSdPlot(fit)
metabs_df <- as.data.frame(vsn::predict(fit,as.matrix(metabs_df)))

pheatmap(metabs_df[complete.cases(metabs_df),], show_colnames = T, show_rownames = F, cluster_rows = F) 

to_write <- metabs_df
to_write$metabolite <- row.names(to_write)
to_write <- to_write[,c(59,1:58)]
to_write$metabolite <- gsub(",","",to_write$metabolite)

write_csv(to_write, file = "data/metabolomic/metabolomic_clean_vsn.csv")

samples <- names(metabs_df)
metabs <- row.names(metabs_df)

SDs <- apply(metabs_df,1,function(x){sd(x,na.rm = T)})
means <- rowMeans(metabs_df, na.rm = T)

metabs_df_scaled <- (metabs_df - means) / SDs
pheatmap(metabs_df_scaled[complete.cases(metabs_df_scaled),], show_colnames = T, show_rownames = F, cluster_rows = F) 

hist(as.numeric(unlist(metabs_df)), breaks = 1000) 
