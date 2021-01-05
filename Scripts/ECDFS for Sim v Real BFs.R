############################################## Create ECDFs for Sim v Real BayPass Output #################################################

#ECDFs for simulation vs real BayPass output by environmental variable
#includes Komolgorov-Smirnov test

################################################################################################################################################

######## Set-up ########

remove(list = ls())

#Set working directory
setwd("C:/Users/rclar/Dropbox/Pinsky_Lab/Transcriptome_Proj/R_scripts/A_clarkii_transcriptomics/")
getwd()

#load libraries
library(tidyverse)
library(ggpubr)

#read in data
BFs_KS <- read_csv("Data/K-S_test_BF_data.csv", col_names = TRUE)

################################################################################################################################################

######## ECDFs & K-S Tests ########
#using average BFs from all 10 runs (not median) per locus

#SSS_mean analysis
#subset data for plot
SSS_mean_df <- BFs_KS[which(BFs_KS$cov == "SSS_mean"), ]

#subset data for K-S test
SSS_mean_permuted <- SSS_mean_df[which(SSS_mean_df$dataset == "permuted"), ] 
SSS_mean_permuted <- SSS_mean_permuted$BF

SSS_mean_real <- SSS_mean_df[which(SSS_mean_df$dataset == "real"), ] 
SSS_mean_real <- SSS_mean_real$BF

#run K-S test for SSS mean dataset
SSS_mean_kstest <- ks.test(SSS_mean_permuted, SSS_mean_real)

#plot distributions
SSS_ecdf_plot <- ggplot(data = SSS_mean_df, aes(x = BF, color = dataset)) + stat_ecdf(geom = "step", size = 1.5) + 
  geom_hline(yintercept = 1, size = 1, linetype = "dashed", color = "#666666") + 
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", color = "#666666") + 
  annotate("text", x = 40, y = 0.65, label = "D = 0.622", size = 12) + 
  annotate("text", x = 40, y = 0.57, label = "p < 0.001", size = 12) + 
  annotate("text", x = 0.0, y = 0.95, label = "A", size = 18)
SSS_ecdf_plot_annotated <- SSS_ecdf_plot + theme_bw() + labs(x = "SSS Mean BF", y = "Proportion of BFs") + 
  scale_x_continuous(limits = c(0, 55)) + scale_color_manual(values = c("#999999", "#E69F00"), labels = c("permuted", "observed")) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(size = 1), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        axis.title = element_text(size = 24, face = "bold"), legend.position = "top", 
        legend.text = element_text(size = 20), legend.title = element_blank())
SSS_ecdf_plot_annotated

#SST_mean analysis
#subset data for plot
SST_mean_df <- BFs_KS[which(BFs_KS$cov == "SST_mean"), ]

#subset data for K-S test
SST_mean_permuted <- SST_mean_df[which(SST_mean_df$dataset == "permuted"), ] 
SST_mean_permuted <- SST_mean_permuted$BF

SST_mean_real <- SST_mean_df[which(SST_mean_df$dataset == "real"), ] 
SST_mean_real <- SST_mean_real$BF

#run K-S test for SST mean dataset
SST_mean_kstest <- ks.test(SST_mean_permuted, SST_mean_real)

#plot distributions
SST_mean_ecdf_plot <- ggplot(data = SST_mean_df, aes(x = BF, color = dataset)) + stat_ecdf(geom = "step", size = 1.5) + 
  geom_hline(yintercept = 1, size = 1, linetype = "dashed", color = "#666666") + 
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", color = "#666666") + 
  annotate("text", x = 40, y = 0.65, label = "D = 0.868", size = 12) + 
  annotate("text", x = 40, y = 0.57, label = "p < 0.001", size = 12) + 
  annotate("text", x = 0.0, y = 0.95, label = "B", size = 18)
SST_mean_ecdf_plot_annotated <- SST_mean_ecdf_plot + theme_bw() + labs(x = "SST Mean BF", y = "Proportion of BFs") + 
  scale_x_continuous(limits = c(0, 55)) + scale_color_manual(values = c("#999999", "#E69F00"), labels = c("permuted", "observed")) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(size = 1), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        axis.title = element_text(size = 24, face = "bold"), legend.position = "top", 
        legend.text = element_text(size = 20), legend.title = element_blank())
SST_mean_ecdf_plot_annotated

#SST_min analysis
#subset data for plot
SST_min_df <- BFs_KS[which(BFs_KS$cov == "SST_min"), ]

#subset data for K-S test
SST_min_permuted <- SST_min_df[which(SST_min_df$dataset == "permuted"), ] 
SST_min_permuted <- SST_min_permuted$BF

SST_min_real <- SST_min_df[which(SST_min_df$dataset == "real"), ] 
SST_min_real <- SST_min_real$BF

#run K-S test for SST min dataset
SST_min_kstest <- ks.test(SST_min_permuted, SST_min_real)

#plot distributions
SST_min_ecdf_plot <- ggplot(data = SST_min_df, aes(x = BF, color = dataset)) + stat_ecdf(geom = "step", size = 1.5) + 
  geom_hline(yintercept = 1, size = 1, linetype = "dashed", color = "#666666") + 
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", color = "#666666") + 
  annotate("text", x = 40, y = 0.65, label = "D = 0.994", size = 12) + 
  annotate("text", x = 40, y = 0.57, label = "p < 0.001", size = 12) + 
  annotate("text", x = 0.0, y = 0.95, label = "C", size = 18)
SST_min_ecdf_plot_annotated <- SST_min_ecdf_plot + theme_bw() + labs(x = "SST Min BF", y = "Proportion of BFs") + 
  scale_x_continuous(limits = c(0, 55)) + scale_color_manual(values = c("#999999", "#E69F00"), labels = c("permuted", "observed")) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(size = 1), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        axis.title = element_text(size = 24, face = "bold"), legend.position = "top", 
        legend.text = element_text(size = 20), legend.title = element_blank())
SST_min_ecdf_plot_annotated

#SST_max analysis
#subset data for plot
SST_max_df <- BFs_KS[which(BFs_KS$cov == "SST_max"), ]

#subset data for K-S test
SST_max_permuted <- SST_max_df[which(SST_max_df$dataset == "permuted"), ] 
SST_max_permuted <- SST_max_permuted$BF

SST_max_real <- SST_max_df[which(SST_max_df$dataset == "real"), ] 
SST_max_real <- SST_max_real$BF

#run K-S test for SST max dataset
SST_max_kstest <- ks.test(SST_max_permuted, SST_max_real)

#plot distributions
SST_max_ecdf_plot <- ggplot(data = SST_max_df, aes(x = BF, color = dataset)) + stat_ecdf(geom = "step", size = 1.5) + 
  geom_hline(yintercept = 1, size = 1, linetype = "dashed", color = "#666666") + 
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", color = "#666666") + 
  annotate("text", x = 40, y = 0.65, label = "D = 0.225", size = 12) + 
  annotate("text", x = 40, y = 0.57, label = "p < 0.001", size = 12) + 
  annotate("text", x = 0.0, y = 0.95, label = "D", size = 18)
SST_max_ecdf_plot_annotated <- SST_max_ecdf_plot + theme_bw() + labs(x = "SST Max BF", y = "Proportion of BFs") + 
  scale_x_continuous(limits = c(0, 55)) + scale_color_manual(values = c("#999999", "#E69F00"), labels = c("permuted", "observed")) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(size = 1), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        axis.title = element_text(size = 24, face = "bold"), legend.position = "top", 
        legend.text = element_text(size = 20), legend.title = element_blank())
SST_max_ecdf_plot_annotated

#Lat analysis
#subset data for plot
Lat_df <- BFs_KS[which(BFs_KS$cov == "Lat"), ]

#subset data for K-S test
Lat_permuted <- Lat_df[which(Lat_df$dataset == "permuted"), ] 
Lat_permuted <- Lat_permuted$BF

Lat_real <- Lat_df[which(Lat_df$dataset == "real"), ] 
Lat_real <- Lat_real$BF

#run K-S test for Lat dataset
Lat_kstest <- ks.test(Lat_permuted, Lat_real)

#plot distributions
Lat_ecdf_plot <- ggplot(data = Lat_df, aes(x = BF, color = dataset)) + stat_ecdf(geom = "step", size = 1.5) + 
  geom_hline(yintercept = 1, size = 1, linetype = "dashed", color = "#666666") + 
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", color = "#666666") + 
  annotate("text", x = 40, y = 0.65, label = "D = 0.794", size = 12) + 
  annotate("text", x = 40, y = 0.57, label = "p < 0.001", size = 12) + 
  annotate("text", x = 0.0, y = 0.95, label = "E", size = 18)
Lat_ecdf_plot_annotated <- Lat_ecdf_plot + theme_bw() + labs(x = "Lat BF", y = "Proportion of BFs") + 
  scale_x_continuous(limits = c(0, 55)) + scale_color_manual(values = c("#999999", "#E69F00"), labels = c("permuted", "observed")) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(size = 1), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 20, color = "black"), 
        axis.title = element_text(size = 24, face = "bold"), legend.position = "top", 
        legend.text = element_text(size = 20), legend.title = element_blank())
Lat_ecdf_plot_annotated

#arrange all plots on same page
ggarrange(SSS_ecdf_plot_annotated, SST_mean_ecdf_plot_annotated, 
          SST_min_ecdf_plot_annotated, SST_max_ecdf_plot_annotated, 
          Lat_ecdf_plot_annotated, ncol = 2, nrow =3)
