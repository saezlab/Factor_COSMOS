library(readr)
library(tidyr)
library(reshape2)
library(pheatmap)
library(vsn)

# Metabolomic Data From Metabolon - data averaged from triplicate experiments (from https://wiki.nci.nih.gov/display/NCIDTPdata/Molecular+Target+Data)
#metabs <- as.data.frame(read_csv("data/metabolomic/WEB_DATA_METABOLON.TXT")) 

# Metabolomic Data From Metabolon - individual data from each of the triplicate experiments (from https://wiki.nci.nih.gov/display/NCIDTPdata/Molecular+Target+Data)
metabs <- as.data.frame(read_csv("data/metabolomic/WEB_DATA_METABOLON_ALL.TXT")) 

# Remove unclear metabolite measurements
metabs <- metabs[!grepl("Isobar",metabs$TITLE),]
metabs <- metabs[!grepl("^X[-]",metabs$TITLE),]
metabs <- metabs[!grepl("^Possible",metabs$TITLE),]

# Summarize triplicates by calculating mean of sample values 
metabs <- metabs %>%
  group_by(TITLE,cellname, panelnbr, cellnbr, pname) %>%
  summarize(mean=mean(Value), SD=sd(Value)) %>%
  as.data.frame()

# Inspect distribution of values
ggplot(metabs, aes(x = cellname, y = log2(mean))) + geom_violin() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Remove extreme values (outliers)
metabs[abs(log2(metabs$mean)) > 5,"mean"] <- NA  

ggplot(metabs, aes(x = cellname, y = log2(mean))) + geom_violin() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

plotNormalHistogram(log2(metabs$mean),breaks = 1000, xlim = c(-10,12))

# Create Sample x Metabolite dataframe
metabs_df <- reshape2::dcast(metabs, formula = TITLE~cellname, value.var = "mean", fill = NA)
rownames(metabs_df) <- metabs_df[,1]
metabs_df <- metabs_df[,-1]

# Normalize dataframe using vsn with glog transformation
fit <- vsnMatrix(as.matrix(metabs_df))
# Visually verify whether there is a dependence of the standard deviation (or variance) on the mean --> desired: no variance-mean dependence = horizontal line.
meanSdPlot(fit) 
# If decent, apply vsn transformation to data
metabs_df <- as.data.frame(vsn::predict(fit,as.matrix(metabs_df)))

# Check violins again
metabs_df_melted <- reshape2::melt(metabs_df)

ggplot(metabs_df_melted, aes(x = variable, y = value)) + geom_violin() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


# Metabolites mapping to have consistent naming (https://hmdb.ca HMDB IDs)
MetaboliteToHMDB <- data.frame(
  "common" = rownames(metabs_df),
  "HMDB" = c("HMDB0011747", "HMDB0000755", "HMDB0001294", "HMDB0001508", "HMDB0000012", "HMDB0001191", "HMDB0000355", "HMDB0000479", "HMDB0000807", "HMDB0000272", "HMDB0003464", "HMDB0000267", "HMDB0000076", "HMDB0003701", "HMDB0001173", "HMDB0001316", "HMDB0001494", "HMDB0000034", "HMDB0000050", "HMDB0000161", "HMDB0000462", "HMDB0001548", "HMDB0000208", "HMDB0001893", "HMDB0001123", "HMDB0000646", "HMDB0000517", "HMDB0000052", "HMDB0000168", "HMDB0000191", "HMDB0000056", "HMDB0001352", "HMDB0000902", "HMDB0000030", "HMDB0001847", "HMDB0000062", "HMDB0000067", "HMDB0000094", "HMDB0000562", "HMDB0000574", "HMDB0000089", "HMDB0000095", "HMDB0001151", "HMDB0000283", "HMDB0000099", "HMDB0000742", "HMDB0001321", "HMDB0000622", "HMDB0000134", "HMDB0001049", "HMDB0029147", "HMDB0011741", "HMDB0000625", "HMDB0001401", "HMDB0000148", "HMDB0000641", "HMDB0000125", "HMDB0001051", "HMDB0000139", "HMDB0000131", "HMDB0000123", "HMDB0000132", "HMDB0000133", "HMDB0001397", "HMDB0003466", "HMDB0000177", "HMDB0000130", "HMDB0000965", "HMDB0000157", "HMDB0000195", "HMDB0000211", "HMDB0000213", "HMDB0000172", "HMDB0000715", "HMDB0004041", "HMDB0000086", "HMDB0001851", "HMDB0002320", "HMDB0000190", "HMDB0000687", "HMDB0000156", "HMDB0000691", "HMDB0000169", "HMDB0001078", "HMDB0001389", "HMDB0000696", "HMDB0000212", "HMDB0000230", "HMDB0000638", "HMDB0001015", "HMDB0000220", "HMDB0001539", "HMDB0001325", "HMDB0000904", "HMDB0001406", "HMDB0004825", "HMDB0000207", "HMDB0000214", "HMDB0000218", "HMDB0000472", "HMDB0003229", "HMDB0000210", "HMDB0000159", "HMDB0001429", "HMDB0000263", "HMDB0002243", "HMDB0000245", "HMDB0000162", "HMDB0001414", "HMDB0001545", "HMDB0000250", "HMDB0000232", "HMDB0000618", "HMDB0001185", "HMDB0000939", "HMDB0000187", "HMDB0000126", "HMDB0000247", "HMDB0001257", "HMDB0000254", "HMDB0000956", "HMDB0000251", "HMDB0000806", "HMDB0000167", "HMDB0002231", "HMDB0000725", "HMDB0001124", "HMDB0000929", "HMDB0000306", "HMDB0000158", "HMDB0000300", "HMDB0000294", "HMDB0000289", "HMDB0000296", "HMDB0011720", "HMDB0000883", "HMDB0000292", "HMDB0000299", "HMDB0002917")
)
rownames(metabs_df) <- MetaboliteToHMDB$HMDB
write_csv(MetaboliteToHMDB, file = "data/metabolomic/MetaboliteToHMDB.csv")

# Visualizing clean dataset
pheatmap(metabs_df[complete.cases(metabs_df),], show_colnames = T, show_rownames = F, cluster_rows = F) 

# Save datasets
metabs_df <- cbind("metabolite" = rownames(metabs_df), metabs_df)
write_csv(metabs_df, file = "data/metabolomic/metabolomic_vsn.csv")

metabs_df_clean <- metabs_df[complete.cases(metabs_df),]
metabs_df_clean <- cbind("metabolite" = rownames(metabs_df_clean), metabs_df_clean)
write_csv(metabs_df_clean, file = "data/metabolomic/metabolomic_clean_vsn.csv")

