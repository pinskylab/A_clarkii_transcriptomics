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
temp_data <- read_csv("SST_Lat_Data.csv", col_names = TRUE) #temp data for SST line plot
data_PC_all <- read_csv("../../VCFs_and_PLINK/PCA_stuff/output.hicov2.snps.only eigenvectors.csv") #eigenvector for PCA w/all loci 
data_PC_inhwe <- read_csv("../../VCFs_and_PLINK/PCA_stuff/output.hicov2.snps.only.inhwe_eigenvec.csv") #eigenvector for PCA w/loci only in HWE
data_PC_outliersonly <- read_csv("../../VCFs_and_PLINK/PCA_stuff/output.hicov2.snps.only.outliersonly.strict_eigenvec.csv") #eigenvector for PCA w/only outlier loci
data_PC_mac2 <- read_csv("../../VCFs_and_PLINK/PCA_stuff/output.hicov2.snps.only.mac2_eigenvec.csv") #eigenvector for PCA w/mac >2 loci
data_PC_mac2_inhwe <- read_csv("../../VCFs_and_PLINK/PCA_stuff/output.hicov2.snps.only.mac2.inhwe_eigenvec.csv") #eigenvector for PCA w/mac >2 and in HWE loci
data_PC_mac2_nooutliers <- read_csv("../../VCFs_and_PLINK/PCA_stuff/output.hicov2.snps.only.mac2.nooutliers.strict_eigenvec.csv") #eigenvector for PCA w/mac >2 w/no outlier loci

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

######## PCAs ########

#PCA with all loci
PCA_all <- ggplot(data = data_PC_all, aes(x = PC1, y = PC2, color = Location)) + 
  geom_point(aes(size = 1)) + ggtitle("PCA with all loci") + 
  scale_color_manual(values = c("#999999", "#E69F00", "#0072B2")) + 
  labs(x = "PC1 (explains 11.53% of total variance)", y = "PC2 (explains 7.66% of total variance)")
PCA_all_annotated <- PCA_all + scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = c(1, 0), axis.line = element_line(size = 1), plot.title = element_text(size = 14, face = "bold"), 
        legend.position = c(1, 0.75), legend.text = element_text(size = 20), legend.title = element_text(size = 20), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 14, color = "black"), 
        axis.title = element_text(size = 18, face = "bold")) + guides(color = guide_legend(override.aes = list(size = 4)))
PCA_all_annotated

#PCA with loci in HWE
PCA_inhwe <- ggplot(data = data_PC_inhwe, aes(x = PC1, y = PC2, color = Location)) + 
  geom_point(aes(size = 1)) + ggtitle("PCA with loci in HWE") + 
  scale_color_manual(values = c("#999999", "#E69F00", "#0072B2")) + 
  labs(x = "PC1 (explains 11.29% of total variance)", y = "PC2 (explains 7.39% of total variance)")
PCA_inhwe_annotated <- PCA_inhwe + scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = c(1, 0), axis.line = element_line(size = 1), plot.title = element_text(size = 14, face = "bold"), 
        legend.position = c(1, 0.75), legend.text = element_text(size = 20), legend.title = element_text(size = 20), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 14, color = "black"), 
        axis.title = element_text(size = 18, face = "bold")) + guides(color = guide_legend(override.aes = list(size = 4)))
PCA_inhwe_annotated

#PCA with only outlier loci
PCA_outliersonly <- ggplot(data = data_PC_outliersonly, aes(x = PC1, y = PC2, color = Location)) + 
  geom_point(size = 8) + ggtitle("PCA with only outlier loci") + 
  scale_color_manual(values = c("darkorange2", "dodgerblue4", "goldenrod1")) + 
  labs(x = "PC1 (explains 64.62% of total variance)", y = "PC2 (explains 7.72% of total variance)")
PCA_outliersonly_annotated <- PCA_outliersonly + scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = c(1, 0), axis.line = element_line(size = 1), plot.title = element_blank(), 
        legend.position = c(1, 0.75), legend.text = element_text(size = 26), legend.title = element_text(size = 26), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26, face = "bold")) + guides(color = guide_legend(override.aes = list(size = 8)))
PCA_outliersonly_annotated

#PCA with loci w/ mac >2
PCA_mac2 <- ggplot(data = data_PC_mac2, aes(x = PC1, y = PC2, color = Location)) + 
  geom_point(size = 8) + ggtitle("PCA with loci with mac > 2") + 
  scale_color_manual(values = c("darkorange2", "dodgerblue4", "goldenrod1")) + 
  labs(x = "PC1 (explains 13.25% of total variance)", y = "PC2 (explains 9.8% of total variance)")
PCA_mac2_annotated <- PCA_mac2 + scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = c(1, 0), axis.line = element_line(size = 1), plot.title = element_text(size = 14, face = "bold"), 
        legend.position = c(1, 0.75), legend.text = element_text(size = 26), legend.title = element_text(size = 26), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26, face = "bold")) + guides(color = guide_legend(override.aes = list(size = 8)))
PCA_mac2_annotated

#PCA with loci w/ mac >2 and in HWE
PCA_mac2_inhwe <- ggplot(data = data_PC_mac2_inhwe, aes(x = PC1, y = PC2, color = Location)) + 
  geom_point(size = 8) + 
  scale_color_manual(values = c("dodgerblue4", "goldenrod1", "darkorange2"), labels = c("Japan", "Philippines", "Indonesia")) + 
  labs(x = "PC1 (explains 12.73% of total variance)", y = "PC2 (explains 9.17% of total variance)")
PCA_mac2_inhwe_annotated <- PCA_mac2_inhwe + scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = c(1, 0), axis.line = element_line(size = 1), plot.title = element_blank(), 
        legend.position = c(1, 0.7), legend.text = element_text(size = 26), legend.title = element_text(size = 26), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26, face = "bold")) + guides(color = guide_legend(override.aes = list(size = 8)))
PCA_mac2_inhwe_annotated

#PCA with loci w/ mac >2 and no outliers
PCA_mac2_nooutliers<- ggplot(data = data_PC_mac2_nooutliers, aes(x = PC1, y = PC2, color = Location)) + 
  geom_point(size = 8) + ggtitle("PCA without outlier loci (mac>2)") + 
  scale_color_manual(values = c("dodgerblue4", "goldenrod1", "darkorange2")) + 
  labs(x = "PC1 (explains 12.76% of total variance)", y = "PC2 (explains 9.53% of total variance)")
PCA_mac2_nooutliers_annotated <- PCA_mac2_nooutliers + scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = c(1, 0), axis.line = element_line(size = 1), plot.title = element_text(size = 14, face = "bold"), 
        legend.position = c(1, 0.75), legend.text = element_text(size = 26), legend.title = element_text(size = 26), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26, face = "bold")) + guides(color = guide_legend(override.aes = list(size = 8)))
PCA_mac2_nooutliers_annotated
