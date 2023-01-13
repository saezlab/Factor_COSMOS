library(readxl)
library(readr)

# Load metadata and reformat for consistency
RNA_metadata <- as.data.frame(read_xls("data/metadata/NCI60_CELL_LINE_METADATA.xls", skip = 7, n_max = 61))
RNA_metadata[,1] <- gsub(".*[:]","",RNA_metadata[,1])
RNA_metadata[,1] <- gsub("_","-",RNA_metadata[,1])
RNA_metadata[,1] <- gsub("MDA-MB-231", "MDA-MB-231/ATCC", RNA_metadata[,1])
RNA_metadata[,1] <- gsub("HS578T", "HS 578T", RNA_metadata[,1])
RNA_metadata[,1] <- gsub("T47D", "T-47D", RNA_metadata[,1])
RNA_metadata[,1] <- gsub("COLO205", "COLO 205", RNA_metadata[,1])
RNA_metadata[,1] <- gsub("HL-60", "HL-60(TB)", RNA_metadata[,1])
RNA_metadata[,1] <- gsub("LOXIMVI h", "LOX IMVI", RNA_metadata[,1])
RNA_metadata[,1] <- gsub("MDA-MB-435 i", "MDA-MB-435", RNA_metadata[,1])
RNA_metadata[,1] <- gsub("MDA-N j", "MDA-N", RNA_metadata[,1])
RNA_metadata[,1] <- gsub("A549", "A549/ATCC", RNA_metadata[,1])
RNA_metadata[,1] <- gsub("NCI-ADR-RES k", "NCI/ADR-RES", RNA_metadata[,1])

colnames(RNA_metadata)[1] <- "cell_line"

RNA_metadata <- RNA_metadata[order(RNA_metadata[,7]),]
RNA_metadata[c(1:19),7] <- "Adenocarcinoma"
RNA_metadata[c(20,21),7] <- "ALL"
RNA_metadata[c(23:26),7] <- "Carcinoma"
RNA_metadata[c(31,32),7] <- "Ductal carcinoma-mammary gland"
RNA_metadata[c(34:37),7] <- "Glioblastoma"
RNA_metadata[39,7] <- "Ductal carcinoma-mammary gland"
RNA_metadata[c(40,41),7] <- "Large Cell Carcinoma"
RNA_metadata[c(48:50),7] <- "Malignant melanotic melanoma"
RNA_metadata[52,7] <- "Ductal carcinoma-mammary gland"
RNA_metadata[53,7] <- "APL"
RNA_metadata[54,7] <- "Prostate carcinoma"
RNA_metadata[55:57,7] <- "Renal cell carcinoma"
RNA_metadata <- RNA_metadata[order(RNA_metadata[,2]),]

write_csv(RNA_metadata, file = "data/metadata/RNA_metadata.csv")

