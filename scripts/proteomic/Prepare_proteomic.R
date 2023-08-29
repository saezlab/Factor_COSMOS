library(readxl)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(vsn)

# SWATH values (Sequential Window Acquisition of all THeoretical Mass Spectra intensity values), log10 --> MassSpec data, protein (exclude peptides) from https://discover.nci.nih.gov/cellminer/datasets.do (https://discover.nci.nih.gov/cellminer/loadDownload.do)
Prot_log10_SWATH <- as.data.frame(read_excel("data/proteomic/Protein__SWATH_(Mass_spectrometry)_Protein.xls", skip = 10, sheet = 1))

# Extract Proteomics measurements and remove double entries by calculating the mean of both measurements, 
Prot_log10_SWATH <- Prot_log10_SWATH[,-c(1,3:9,43)]
Prot_log10_SWATH <- Prot_log10_SWATH[!Prot_log10_SWATH[,1]=="-",]
colnames(Prot_log10_SWATH)[1] <- "Gene" 

Prot_log10_SWATH <-  Prot_log10_SWATH %>%
  gather(., cellline, value, 2:60) %>%
  group_by(Gene, cellline) %>%
  summarize(mean=mean(as.numeric(value), na.rm = T)) %>%
  pivot_wider(id_cols = Gene ,names_from = cellline, values_from = mean) %>%
  as.data.frame()

# Align cell line names
rownames(Prot_log10_SWATH) <- Prot_log10_SWATH[,1]
Prot_log10_SWATH <- Prot_log10_SWATH[,-1]
names(Prot_log10_SWATH) <- gsub(".*[:]","",names(Prot_log10_SWATH))
Prot_log10_SWATH <- Prot_log10_SWATH[,order(colnames(Prot_log10_SWATH))]
names(Prot_log10_SWATH) <- gsub("MDA-MB-231", "MDA-MB-231/ATCC", names(Prot_log10_SWATH))
names(Prot_log10_SWATH) <- gsub("RXF 393", "RXF-393", names(Prot_log10_SWATH))

# Filter dataset: Inspect distribution and remove outliers
hist(as.numeric(unlist(Prot_log10_SWATH)), breaks = 1000) 
Prot_log10_SWATH_melted <- reshape2::melt(Prot_log10_SWATH)
ggplot(Prot_log10_SWATH_melted, aes(x = variable, y = value)) + geom_violin() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
# --> looks already pretty good; to transform back to original values: Prot_log10_SWATH <- (10^Prot_log10_SWATH) - 1 and then perform vsn

Prot_log10_SWATH_clean <- Prot_log10_SWATH
hist.data <- hist(as.numeric(unlist(Prot_log10_SWATH_clean)), breaks = 1000)
hist.data$counts = log10(hist.data$counts)
plot(hist.data, ylim = c(0,6), ylab='log10(Frequency)', xlab = "Prot log10(intensity)")


Prot_log10_SWATH_clean[Prot_log10_SWATH_clean==0] <- NA
dim(Prot_log10_SWATH_clean[rowSums(is.na(Prot_log10_SWATH_clean))<dim(Prot_log10_SWATH_clean)[2]*(1/3),])

Prot_log10_SWATH_clean_melted <- reshape2::melt(Prot_log10_SWATH_clean)
ggplot(Prot_log10_SWATH_clean_melted, aes(x = variable, y = value)) + geom_violin() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Visualize dataset
pheatmap(Prot_log10_SWATH_clean[complete.cases(Prot_log10_SWATH_clean),], show_colnames = T, show_rownames = F, cluster_rows = F) 

dim(Prot_log10_SWATH_clean)
# Save overlapping clean or whole measurements
Prot_log10_SWATH <- cbind("Proteins" = rownames(Prot_log10_SWATH), Prot_log10_SWATH)
Prot_log10_SWATH_clean <- cbind("Proteins" = rownames(Prot_log10_SWATH_clean), Prot_log10_SWATH_clean)
write_csv(Prot_log10_SWATH,file = "data/proteomic/Prot_log10_SWATH.csv")
write_csv(Prot_log10_SWATH_clean,file = "data/proteomic/Prot_log10_SWATH_clean.csv")



