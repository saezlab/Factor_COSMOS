library(readr)

RNA_log2_FPKM <- as.data.frame(
  read_delim("data/RNA/RNA_log2_FPKM.csv", 
                            delim = ";", escape_double = FALSE, trim_ws = TRUE))

RNA_log2_FPKM <- RNA_log2_FPKM[,-c(2:6)]
names(RNA_log2_FPKM) <- gsub(".*[:]","",names(RNA_log2_FPKM))

RNA_log2_FPKM <- RNA_log2_FPKM[complete.cases(RNA_log2_FPKM),]

row.names(RNA_log2_FPKM) <- RNA_log2_FPKM$`Gene name d`
RNA_log2_FPKM <- RNA_log2_FPKM[,-1]
RNA_log2_FPKM <- apply(RNA_log2_FPKM,2,function(x){
  gsub(",",".",x)
}) 
RNA_log2_FPKM <- as.data.frame(RNA_log2_FPKM)

hist(as.numeric(unlist(RNA_log2_FPKM)), breaks = 1000) 

RNA_log2_FPKM[RNA_log2_FPKM < 1] <- NA
hist(as.numeric(unlist(RNA_log2_FPKM)), breaks = 1000) 
RNA_log2_FPKM <- RNA_log2_FPKM[rowSums(is.na(RNA_log2_FPKM)) < 45,]
# RNA_log2_FPKM <- RNA_log2_FPKM[rowSums(RNA_log2_FPKM == 0) < 45,]
hist(as.numeric(unlist(RNA_log2_FPKM)), breaks = 1000) 

metabolomic_clean_vsn <- as.data.frame(read_csv("data/metabolomic/metabolomic_clean_vsn.csv"))

cell_lines <- intersect(names(metabolomic_clean_vsn[,-1]),names(RNA_log2_FPKM))

metabolomic_clean_vsn <- metabolomic_clean_vsn[,c("metabolite",cell_lines)]

RNA_log2_FPKM <- RNA_log2_FPKM[,cell_lines]
RNA_log2_FPKM$gene <- row.names(RNA_log2_FPKM)

RNA_log2_FPKM <- RNA_log2_FPKM[,c(57,1:56)]

write_csv(RNA_log2_FPKM,file = "data/RNA/RNA_log2_FPKM_cleaned_common.csv")
write_csv(metabolomic_clean_vsn, file = "data/metabolomic/metabolomic_clean_vsn_common.csv")
