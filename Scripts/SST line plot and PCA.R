################################################ Create Line Plot of SST env variables & PCAs ##################################################

#ine plot with annual mean min and max SS temp in C
#PCAs with all loci, loci only in HWE, outlier loci excluded, and only outlier loci

################################################################################################################################################

######## Set-up ########

remove(list = ls())

#set working directory
setwd("C:/Users/Rene/Dropbox/Pinsky_Lab/Transcriptome_Proj/R_scripts/A_clarkii_transcriptomics/")
getwd()

#load libraries
library(tidyverse)

#read in data
temp_data <- read_csv("Data/SST_Lat_Data.csv", col_names = TRUE) #temp data for SST line plot
data_PC_mac2 <- read_csv("../../VCFs_and_PLINK/PCA_stuff/mac2/output.hicov2.snps.only.mac2_eigenvec.csv") #eigenvector for PCA w/mac >2 loci
data_PC_mac2_inhwe <- read_csv("../../VCFs_and_PLINK/PCA_stuff/mac2/output.hicov2.snps.only.mac2.inhwe_eigenvec.csv") #eigenvector for PCA w/mac >2 and in HWE loci
data_PC_mac2_nooutliers <- read_csv("../../VCFs_and_PLINK/PCA_stuff/mac2/output.hicov2.snps.only.mac2.nooutliers_eigenvec.csv") #eigenvector for PCA w/mac >2 w/no outlier loci
data_PC_mac2_outliersonly <- read_csv("../../VCFs_and_PLINK/PCA_stuff/mac2/output.hicov2.snps.only.mac2.outliersonly_eigenvec.csv") #eigenvector for PCA w/only outlier loci
data_PC_mac1 <- read_csv("../../VCFs_and_PLINK/PCA_stuff/mac1/output.hicov2.snps.only.mac1_eigenvec.csv")
data_PC_mac1_inhwe <- read_csv("../../VCFs_and_PLINK/PCA_stuff/mac1/output.hicov2.snps.only.mac1.inhwe_eigenvec.csv")
data_PC_mac1_nooutliers <- read_csv("../../VCFs_and_PLINK/PCA_stuff/mac1/output.hicov2.snps.only.mac1.nooutliers_eigenvec.csv")
data_PC_mac1_outliersonly <- read_csv("../../VCFs_and_PLINK/PCA_stuff/mac1/output.hicov2.snps.only.mac1.outliersonly_eigenvec.csv")

#################################################################################################################################################

######## Line plot of SST values ########

#temp_data$Population <- factor(temp_data$Population, levels = c("Indonesia", "Philippines", "Japan"))

#create line plot of SSTs (x-axis as population)
SST_line_plot <- ggplot(data = temp_data, aes(x = Population, y = Temperature, group = Temp)) + 
  geom_line(aes(color = Temp, linetype = Temp), size = 1.5) + 
  scale_color_manual(values = c("#990000", "#666666", "#0066CC")) + scale_linetype_manual(values = c("solid", "dashed", "dotted")) + 
  geom_point(aes(color = Temp, shape = Temp), size = 3.5) + scale_shape_manual(values = c(15, 16, 17))
SST_line_plot_annotated <- SST_line_plot + scale_x_discrete(labels = c("Japan", "Philippines", "Indonesia")) + 
  ggtitle("Sea Surface Temperature Across Sites") + ylab("Temperature (C)") + theme_minimal() + 
  theme(panel.border = element_blank(), panel.grid.major = element_line(size = 1.5), 
        legend.justification = c(1, 0), axis.line = element_blank(), plot.title = element_text(face = "bold"), 
        legend.position = c(0.83, 0.22), legend.text = element_text(size = 16), legend.title = element_blank(), 
        axis.ticks = element_line(size = 2), axis.ticks.length = unit(2, units = "mm"), 
        axis.text = element_text(size = 14, color = "black"), axis.title = element_text(size = 18, face = "bold")) 
SST_line_plot_annotated

#create line plot of SSTs (x-axis as SST)
SST_line_plot_2 <- ggplot(data = temp_data, aes(x = Temp, y = Temperature, group = Population)) + 
  geom_line(aes(color = Population, linetype = Population), size = 1.5) + 
  geom_point(aes(color = Population, shape = Population), size = 3.5) + 
  scale_color_manual(values = c("#0066CC", "#666666", "#990000"), labels = c("Japan", "Philippines", "Indonesia")) + 
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), labels = c("Japan", "Philippines", "Indonesia")) + 
  scale_shape_manual(values = c(15, 16, 17), labels = c("Japan", "Philippines", "Indonesia")) 
SST_line_plot_2_annotated <- SST_line_plot_2 + 
  ggtitle("Sea Surface Temperature Across Sites") + ylab("Temperature (C)") + theme_minimal() + 
  theme(panel.border = element_blank(), panel.grid.major = element_line(size = 1.5), 
        legend.justification = c(1, 0), axis.line = element_blank(), plot.title = element_text(face = "bold"),
        legend.position = c(0.85, 0.2), legend.text = element_text(size = 16), legend.title = element_blank(), 
        axis.ticks = element_line(size = 2), axis.ticks.length = unit(2, units = "mm"),
        axis.text = element_text(size = 14, color = "black"), axis.title.y = element_text(size = 18, face = "bold"), 
        axis.title.x = element_blank())

SST_bar_mean <- ggplot(data = temp_data, aes(x = Population, y = Temperature, group = Temp)) +
  geom_line()

#################################################################################################################################################

######## mac2 PCAs ########

#PCA with loci w/ mac >2
PCA_mac2 <- ggplot(data = data_PC_mac2, aes(x = PC1, y = PC2, color = Location)) + 
  geom_point(size = 8) + ggtitle("PCA with loci with mac > 2") + 
  scale_color_manual(values = c("dodgerblue4", "goldenrod1", "darkorange2"), labels = c("Japan", "Philippines", "Indonesia")) + 
  labs(x = "PC1 (explains 13.25% of total variance)", y = "PC2 (explains 9.8% of total variance)")
PCA_mac2_annotated <- PCA_mac2 + scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = c(1, 0), axis.line = element_line(size = 1), plot.title = element_text(size = 14, face = "bold"), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_text(size = 26), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26, face = "bold")) + guides(color = guide_legend(override.aes = list(size = 8)))
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
        axis.title = element_text(size = 26, face = "bold")) + guides(color = guide_legend(override.aes = list(size = 8)))
PCA_mac2_inhwe_annotated

#PCA with loci w/ mac >2 and no outliers
PCA_mac2_nooutliers<- ggplot(data = data_PC_mac2_nooutliers, aes(x = PC1, y = PC2, color = Location)) + 
  geom_point(size = 8) + ggtitle("PCA without outlier loci (mac > 2)") + 
  scale_color_manual(values = c("dodgerblue4", "goldenrod1", "darkorange2"), labels = c("Japan", "Philippines", "Indonesia")) + 
  labs(x = "PC1 (explains 12.75% of total variance)", y = "PC2 (explains 9.54% of total variance)")
PCA_mac2_nooutliers_annotated <- PCA_mac2_nooutliers + scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = c(1, 0), axis.line = element_line(size = 1), plot.title = element_text(size = 14, face = "bold"), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_text(size = 26), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26, face = "bold")) + guides(color = guide_legend(override.aes = list(size = 8)))
PCA_mac2_nooutliers_annotated

#PCA with loci w/ mac >2 and only outliers
PCA_mac2_outliersonly <- ggplot(data = data_PC_mac2_outliersonly, aes(x = PC1, y = PC2, color = Location)) + 
  geom_point(size = 8) + ggtitle("PCA with only outlier loci (mac > 2)") + 
  scale_color_manual(values = c("dodgerblue4", "goldenrod1", "darkorange2"), labels = c("Japan", "Philippines", "Indonesia")) + 
  labs(x = "PC1 (explains 65.66% of total variance)", y = "PC2 (explains 8.14% of total variance)")
PCA_mac2_outliersonly_annotated <- PCA_mac2_outliersonly + scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = c(1, 0), axis.line = element_line(size = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_text(size = 26), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26, face = "bold")) + guides(color = guide_legend(override.aes = list(size = 8)))
PCA_mac2_outliersonly_annotated

#################################################################################################################################################

######## mac1 PCAs ########

#PCA with loci w/ mac >1
PCA_mac1 <- ggplot(data = data_PC_mac1, aes(x = PC1, y = PC2, color = Location)) + 
  geom_point(size = 8) + ggtitle("PCA with loci with mac > 1") + 
  scale_color_manual(values = c("dodgerblue4", "goldenrod1", "darkorange2"), labels = c("Japan", "Philippines", "Indonesia")) + 
  labs(x = "PC1 (explains 11.53% of total variance)", y = "PC2 (explains 7.66% of total variance)")
PCA_mac1_annotated <- PCA_mac1 + scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = c(1, 0), axis.line = element_line(size = 1), plot.title = element_text(size = 14, face = "bold"), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_text(size = 26), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26, face = "bold")) + guides(color = guide_legend(override.aes = list(size = 8)))
PCA_mac1_annotated

#PCA with loci w/ mac >1 and in HWE
PCA_mac1_inhwe <- ggplot(data = data_PC_mac1_inhwe, aes(x = PC1, y = PC2, color = Location)) + 
  geom_point(size = 8) + ggtitle("PCA with loci in HWE (mac > 1") + 
  scale_color_manual(values = c("dodgerblue4", "goldenrod1", "darkorange2"), labels = c("Japan", "Philippines", "Indonesia")) + 
  labs(x = "PC1 (explains 11.28% of total variance)", y = "PC2 (explains 7.21% of total variance)")
PCA_mac1_inhwe_annotated <- PCA_mac1_inhwe + scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = c(1, 0), axis.line = element_line(size = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_text(size = 26), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26, face = "bold")) + guides(color = guide_legend(override.aes = list(size = 8)))
PCA_mac1_inhwe_annotated

#PCA with loci w/ mac >1 and no outliers
PCA_mac1_nooutliers<- ggplot(data = data_PC_mac1_nooutliers, aes(x = PC1, y = PC2, color = Location)) + 
  geom_point(size = 8) + ggtitle("PCA without outlier loci (mac > 1)") + 
  scale_color_manual(values = c("dodgerblue4", "goldenrod1", "darkorange2"), labels = c("Japan", "Philippines", "Indonesia")) + 
  labs(x = "PC1 (explains 11.48% of total variance)", y = "PC2 (explains 7.52% of total variance)")
PCA_mac1_nooutliers_annotated <- PCA_mac1_nooutliers + scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = c(1, 0), axis.line = element_line(size = 1), plot.title = element_text(size = 14, face = "bold"), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_text(size = 26), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26, face = "bold")) + guides(color = guide_legend(override.aes = list(size = 8)))
PCA_mac1_nooutliers_annotated

#PCA with loci w/mac >1 and only outliers
PCA_mac1_outliersonly <- ggplot(data = data_PC_mac1_outliersonly, aes(x = PC1, y = PC2, color = Location)) + 
  geom_point(size = 8) + ggtitle("PCA with only outlier loci (mac > 1)") + 
  scale_color_manual(values = c("dodgerblue4", "goldenrod1", "darkorange2"), labels = c("Japan", "Philippines", "Indonesia")) + 
  labs(x = "PC1 (explains 67.10% of total variance)", y = "PC2 (explains 8.27% of total variance)")
PCA_mac1_outliersonly_annotated <- PCA_mac1_outliersonly + scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = c(1, 0), axis.line = element_line(size = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_text(size = 26), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26, face = "bold")) + guides(color = guide_legend(override.aes = list(size = 8)))
PCA_mac1_outliersonly_annotated
