library(readxl)
library(readr)
library(pheatmap)

# Fragments per kilobase per million reads log2(FPKM + 1) --> RNAseq, composite expression from https://discover.nci.nih.gov/cellminer/datasets.do (https://discover.nci.nih.gov/cellminer/loadDownload.do)
RNA_log2_FPKM <- as.data.frame(read_excel("data/RNA/RNA__RNA_seq_composite_expression.xls", skip = 10, sheet = 1))

# Extract RNAseq measurements and align cell line names
RNA_log2_FPKM <- RNA_log2_FPKM[,-c(2:6)]
rownames(RNA_log2_FPKM) <- RNA_log2_FPKM[,1]
RNA_log2_FPKM <- RNA_log2_FPKM[,-1]
names(RNA_log2_FPKM) <- gsub(".*[:]","",names(RNA_log2_FPKM))
RNA_log2_FPKM <- RNA_log2_FPKM[,order(colnames(RNA_log2_FPKM))]
names(RNA_log2_FPKM) <- gsub("MDA-MB-231", "MDA-MB-231/ATCC", names(RNA_log2_FPKM))
names(RNA_log2_FPKM) <- gsub("RXF 393", "RXF-393", names(RNA_log2_FPKM))

# Filter dataset: Inspect distribution and remove outliers
hist(as.numeric(unlist(RNA_log2_FPKM)), breaks = 1000) 
hist(as.numeric(unlist(RNA_log2_FPKM[RNA_log2_FPKM>0])), breaks = 1000)
hist(as.numeric(unlist(RNA_log2_FPKM[RNA_log2_FPKM>=1])), breaks = 1000) 
RNA_log2_FPKM_clean <- RNA_log2_FPKM
RNA_log2_FPKM_clean[RNA_log2_FPKM_clean < 1] <- NA
hist(as.numeric(unlist(RNA_log2_FPKM_clean)), breaks = 1000, xlim = c(0,15)) 
RNA_log2_FPKM_clean <- RNA_log2_FPKM_clean[rowSums(is.na(RNA_log2_FPKM_clean))<45,]
hist(as.numeric(unlist(RNA_log2_FPKM_clean)), breaks = 1000, xlim = c(0,15)) 

# Visualize clean dataset
pheatmap(RNA_log2_FPKM_clean[complete.cases(RNA_log2_FPKM_clean),-1], show_colnames = T, show_rownames = F, cluster_rows = F) 

# Save overlapping clean or whole measurements
RNA_log2_FPKM <- cbind("Genes" = rownames(RNA_log2_FPKM), RNA_log2_FPKM)
RNA_log2_FPKM_clean <- cbind("Genes" = rownames(RNA_log2_FPKM_clean), RNA_log2_FPKM_clean)
write_csv(RNA_log2_FPKM,file = "data/RNA/RNA_log2_FPKM.csv")
write_csv(RNA_log2_FPKM_clean,file = "data/RNA/RNA_log2_FPKM_clean.csv")


