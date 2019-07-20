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
data_PC_all <- read_csv("output.hicov2.snps.only eigenvectors.csv") #eigenvector for PCA w/all loci 
data_PC_inhwe <- read_csv("output.hicov2.snps.only.inhwe_eigenvec.csv") #eigenvector for PCA w/loci only in HWE
data_PC_nooutliers <- read_csv("output.hicov2.snps.only.nooutliers.strict_eigenvec.csv") #eigenvector for PCA w/no outlier loci
data_PC_outliersonly <- read_csv("output.hicov2.snps.only.outliersonly.strict_eigenvec.csv") #eigenvector for PCA w/only outlier loci
data_PC_mac2 <- read_csv("output.hicov2.snps.only.mac2_eigenvec.csv") #eigenvector for PCA w/mac >2 loci
data_PC_mac2_inhwe <- read_csv("output.hicov2.snps.only.mac2.inhwe_eigenvec.csv") #eigenvector for PCA w/mac >2 and in HWE loci

#################################################################################################################################################

######## Line plot of SST values ########

#temp_data$Population <- factor(temp_data$Population, levels = c("Indonesia", "Philippines", "Japan"))

#create line plot of SSTs
SST_line_plot <- ggplot(data = temp_data, aes(x = Population, y = Temperature, group = Temp)) + 
  geom_line(aes(color = Temp, size = Temp)) + 
  geom_point(aes(color = Temp, size = Temp)) + scale_size_manual(values = c(1.5, 1.5, 1.5))
SST_line_plot_annotated <- SST_line_plot + scale_x_discrete(labels = c("Indonesia", "Philippines", "Japan")) + 
  ggtitle("Sea Surface Temperature Across Sites") + ylab("Temperature (C)") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = c(1, 0), axis.line = element_line(size = 1), plot.title = element_text(face = "bold"), 
        legend.position = c(0.33, 0.25), legend.text = element_text(size = 16), legend.title = element_text(size = 16), 
        axis.text = element_text(size = 14, color = "black"), axis.title = element_text(size = 18, face = "bold")) 
SST_line_plot_annotated

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

#PCA without outlier loci
PCA_nooutliers <- ggplot(data = data_PC_nooutliers, aes(x = PC1, y = PC2, color = Location)) + 
  geom_point(aes(size = 1)) + ggtitle("PCA without outlier loci") + 
  scale_color_manual(values = c("#999999", "#E69F00", "#0072B2")) + 
  labs(x = "PC1 (explains 11.5% of total variance)", y = "PC2 (explains 7.51% of total variance)")
PCA_nooutliers_annotated <- PCA_nooutliers + scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = c(1, 0), axis.line = element_line(size = 1), plot.title = element_text(size = 14, face = "bold"), 
        legend.position = c(1, 0.75), legend.text = element_text(size = 20), legend.title = element_text(size = 20), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 14, color = "black"), 
        axis.title = element_text(size = 18, face = "bold")) + guides(color = guide_legend(override.aes = list(size = 4)))
PCA_nooutliers_annotated

#PCA with only outlier loci
PCA_outliersonly <- ggplot(data = data_PC_outliersonly, aes(x = PC1, y = PC2, color = Location)) + 
  geom_point(aes(size = 1)) + ggtitle("PCA with only outlier loci") + 
  scale_color_manual(values = c("#999999", "#E69F00", "#0072B2")) + 
  labs(x = "PC1 (explains 64.62% of total variance)", y = "PC2 (explains 7.72% of total variance)")
PCA_outliersonly_annotated <- PCA_outliersonly + scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = c(1, 0), axis.line = element_line(size = 1), plot.title = element_text(size = 14, face = "bold"), 
        legend.position = c(1, 0.75), legend.text = element_text(size = 20), legend.title = element_text(size = 20), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 14, color = "black"), 
        axis.title = element_text(size = 18, face = "bold")) + guides(color = guide_legend(override.aes = list(size = 4)))
PCA_outliersonly_annotated

#PCA with loci w/ mac >2
PCA_mac2 <- ggplot(data = data_PC_mac2, aes(x = PC1, y = PC2, color = Location)) + 
  geom_point(aes(size = 1)) + ggtitle("PCA with loci with mac > 2") + 
  scale_color_manual(values = c("#999999", "#E69F00", "#0072B2")) + 
  labs(x = "PC1 (explains 13.25% of total variance)", y = "PC2 (explains 9.8% of total variance)")
PCA_mac2_annotated <- PCA_mac2 + scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = c(1, 0), axis.line = element_line(size = 1), plot.title = element_text(size = 14, face = "bold"), 
        legend.position = c(1, 0.75), legend.text = element_text(size = 20), legend.title = element_text(size = 20), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 14, color = "black"), 
        axis.title = element_text(size = 18, face = "bold")) + guides(color = guide_legend(override.aes = list(size = 4)))
PCA_mac2_annotated

#PCA with loci w/ mac >2 and in HWE
PCA_mac2_inhwe <- ggplot(data = data_PC_mac2_inhwe, aes(x = PC1, y = PC2, color = Location)) + 
  geom_point(aes(size = 1)) + ggtitle("PCA with loci with mac > 2 and in HWE") + 
  scale_color_manual(values = c("#999999", "#E69F00", "#0072B2")) + 
  labs(x = "PC1 (explains 12.73% of total variance)", y = "PC2 (explains 9.17% of total variance)")
PCA_mac2_inhwe_annotated <- PCA_mac2_inhwe + scale_size(guide = "none") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = c(1, 0), axis.line = element_line(size = 1), plot.title = element_text(size = 14, face = "bold"), 
        legend.position = c(1, 0.75), legend.text = element_text(size = 20), legend.title = element_text(size = 20), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 14, color = "black"), 
        axis.title = element_text(size = 18, face = "bold")) + guides(color = guide_legend(override.aes = list(size = 4)))
PCA_mac2_inhwe_annotated
