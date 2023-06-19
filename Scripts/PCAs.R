################################################ Create Line Plot of SST env variables & PCAs ##################################################

#PCAs with all loci, loci only in HWE, outlier loci excluded, and only outlier loci

################################################################################################################################################

######## Set-up ########

remove(list = ls())

#set working directory
setwd("C:/Users/Rene/Dropbox/Pinsky_Lab/Transcriptome_Proj/R_scripts/A_clarkii_transcriptomics/")
getwd()

#load libraries
library(tidyverse) #v.2.0.0

#read in data
data_PC_mac2 <- read_csv("Data/output.hicov2.snps.only.mac2_eigenvec.csv") #eigenvector for PCA w/mac >2 loci
data_PC_mac2_inhwe <- read_csv("Data/output.hicov2.snps.only.mac2.inhwe_eigenvec.csv") #eigenvector for PCA w/mac >2 and in HWE loci
data_PC_mac2_nooutliers <- read_csv("Data/output.hicov2.snps.only.mac2.nooutlierswRDA_eigenvec.csv") #eigenvector for PCA w/mac >2 w/no outlier loci
data_PC_mac2_outliersonly <- read_csv("Data/output.hicov2.snps.only.mac2.outliersonlywRDA_eigenvec.csv") #eigenvector for PCA w/only outlier loci

#################################################################################################################################################

######## mac2 PCAs ########

#PCA with all loci w/ mac >2
PCA_mac2 <- ggplot(data = data_PC_mac2, aes(x = PC1, y = PC2, color = Location)) + 
  geom_point(size = 8) + ggtitle("PCA with loci with mac > 2") + 
  scale_color_manual(values = c("dodgerblue4", "goldenrod1", "darkorange2"), labels = c("Japan", "Philippines", "Indonesia")) + 
  labs(x = "PC1 (explains 13.25% of total variance)", y = "PC2 (explains 9.8% of total variance)")
PCA_mac2_annotated <- PCA_mac2 + scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = c(1, 0), axis.line = element_line(size = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_text(size = 26), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26)) + guides(color = guide_legend(override.aes = list(size = 8)))
PCA_mac2_annotated

#PCA with loci w/ mac >2 and in HWE
PCA_mac2_inhwe <- ggplot(data = data_PC_mac2_inhwe, aes(x = PC1, y = PC2, color = Location)) + 
  geom_point(size = 8) + ggtitle("PCA with loci in HWE (mac > 2") + 
  scale_color_manual(values = c("dodgerblue4", "goldenrod1", "darkorange2"), labels = c("Japan", "Philippines", "Indonesia")) + 
  labs(x = "PC1 (explains 12.73% of total variance)", y = "PC2 (explains 9.17% of total variance)")
PCA_mac2_inhwe_annotated <- PCA_mac2_inhwe + scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = c(1, 0), axis.line = element_line(size = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_text(size = 26), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26)) + guides(color = guide_legend(override.aes = list(size = 8)))
PCA_mac2_inhwe_annotated

#PCA with loci w/ mac >2 and no outliers
PCA_mac2_nooutliers<- ggplot(data = data_PC_mac2_nooutliers, aes(x = PC1, y = PC2, color = Location)) + 
  geom_point(size = 8) + ggtitle("PCA without outlier loci (mac > 2)") + 
  scale_color_manual(values = c("dodgerblue4", "goldenrod1", "darkorange2"), labels = c("Japan", "Philippines", "Indonesia")) + 
  labs(x = "PC1 (explains 12.73% of total variance)", y = "PC2 (explains 9.52% of total variance)")
PCA_mac2_nooutliers_annotated <- PCA_mac2_nooutliers + scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = c(1, 0), axis.line = element_line(size = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_text(size = 26), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26)) + guides(color = guide_legend(override.aes = list(size = 8)))
PCA_mac2_nooutliers_annotated

#PCA with loci w/ mac >2 and only outliers
PCA_mac2_outliersonly <- ggplot(data = data_PC_mac2_outliersonly, aes(x = PC1, y = PC2, color = Location)) + 
  geom_point(size = 8) + ggtitle("PCA with only outlier loci (mac > 2)") + 
  scale_color_manual(values = c("dodgerblue4", "goldenrod1", "darkorange2"), labels = c("Japan", "Philippines", "Indonesia")) + 
  labs(x = "PC1 (explains 72.54% of total variance)", y = "PC2 (explains 7.41% of total variance)")
PCA_mac2_outliersonly_annotated <- PCA_mac2_outliersonly + scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = c(1, 0), axis.line = element_line(size = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_text(size = 26), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26)) + guides(color = guide_legend(override.aes = list(size = 8)))
PCA_mac2_outliersonly_annotated
