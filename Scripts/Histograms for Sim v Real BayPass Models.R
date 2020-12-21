############################################## Create Histograms for Sim v Real BayPass Output #################################################

#histograms for simulation vs real BayPass output by environmental variable
#includes Fisher's exact test for mean and median BF for each loci (averaging across 10 simulation v 10 real runs)

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
bfs_mean_all <- read_csv("Data/mac1_sim_v_real_mean.csv", col_names = TRUE) #mean BF across 10 runs (sim and real)
bfs_median_all <- read_csv("Data/mac1_sim_v_real_median.csv", col_names = TRUE) #median BF across 10 runs (sim and real)
BFs <- read_csv("Data/BF_total_mac1.csv", col_names = TRUE)
BFs_KS <- read_csv("Data/K-S_test_BF_data_mac1.csv", col_names = TRUE)

################################################################################################################################################

######## Fisher's exact test with mean data ########

#create 2x2 contingency tables for mean tests
SSS_mean_mean_table <- matrix(c(5030.9, 687.1, 4340, 1378), nrow = 2, ncol = 2) #create matrix with: col 1 = simulated data, col 2 = real data, row 1 = # loci w/BF<10, row 2 = # loci w/BF>10
SST_mean_mean_table <- matrix(c(5168, 550, 4252, 1466), nrow = 2, ncol = 2)
SST_min_mean_table <- matrix(c(5160.8, 557.2, 3950, 1768), nrow = 2, ncol = 2)
SST_max_mean_table <- matrix(c(5541.7, 176.3, 5682, 36), nrow = 2, ncol = 2)
Lat_mean_table <- matrix(c(5075.6, 642.4, 4103, 1615), nrow = 2, ncol = 2)

#Fisher's exact test by env. covariable
#conduct Fisher's Exact Test rounding to nearest integer
FT_SSS_mean_mean <- fisher.test(SSS_mean_mean_table) #p < 0.001 Odds ratio = 2.324 95%CI (2.101, 2.574)
FT_SST_mean_mean <- fisher.test(SST_mean_mean_table)
FT_SST_min_mean <- fisher.test(SST_min_mean_table)
FT_SST_max_mean <- fisher.test(SST_max_mean_table)
FT_lat_mean <- fisher.test(Lat_mean_table)

######## Fisher's exact test with median data ########

#create 2x2 contingency tables for median tests
SSS_mean_median_table <- matrix(c(5095.5, 693.01, 4340, 1378), nrow = 2, ncol = 2)
SST_mean_median_table <- matrix(c(5197, 521, 4252, 1466), nrow = 2, ncol = 2)
SST_min_median_table <- matrix(c(5078, 640, 3950, 1768), nrow = 2, ncol = 2)
SST_max_median_table <- matrix(c(5702.5, 15.5, 5682, 36), nrow = 2, ncol = 2)
Lat_median_table <- matrix(c(5146.5, 571.5, 4103, 1615), nrow = 2, ncol = 2)

#Fisher's exact test by env. covariable
FT_SSS_mean_median <- fisher.test(SSS_mean_median_table) #p < 0.001 Odds ration = 2.335 95%CI (2.111, 2.584)
FT_SST_mean_median <- fisher.test(SST_mean_median_table)
FT_SST_min_median <- fisher.test(SST_min_median_table)
FT_SST_max_median <- fisher.test(SST_max_median_table)
FT_lat_median <- fisher.test(Lat_median_table)

################################################################################################################################################

######## Histograms with mean data ########

#create histograms (1 histogram = 1 environmental covariate with mean for sim & real models) --> 5 histograms total

#Sort dataframe by environmental covariate
SSS_mean_bfs_mean <- bfs_mean_all[1:12, -2]
SST_mean_bfs_mean <- bfs_mean_all[13:24, -2]
SST_min_bfs_mean <- bfs_mean_all[25:36, -2]
SST_max_bfs_mean <- bfs_mean_all[37:48, -2]
lat_bfs_mean <- bfs_mean_all[49:60, -2]

#SSS mean histogram of means
sss_mean_bfs_mean_plot <- ggplot(data = SSS_mean_bfs_mean) + 
  geom_bar(mapping = aes(x = BF_Bin, y = Count, fill = Class), stat = "identity", position = "dodge") #specify that data are already counted (stat="identity") and want two bars for each value (position = "dodge") separated by Class (fill = Class)
sss_mean_mean_plot_annotated <- sss_mean_bfs_mean_plot + theme(legend.position = "top") + 
  annotate("text", x = 5, y = c(4000, 3650, 3300), label = c("p < 0.001", "Odds Ratio = 2.324", "95% CI = (2.101, 2.574)")) + 
  geom_vline(xintercept = 4.5 ) + ggtitle("SSS mean") + 
  scale_x_discrete(breaks = c("bf1", "bf2", "bf3", "bf4", "bf5", "bf6"), 
                   labels = c("BF<0", "0<BF<5", "5<BF<10", "10<BF<15", "15<BF<20", "BF>20")) #relabel X-axis, move legend to top of chart, add p-value and odds ratio, add vertical line dividing two categories and add title

#SST mean histogram of means
sst_mean_bfs_mean_plot <- ggplot(data = SST_mean_bfs_mean) + 
  geom_bar(mapping = aes(x= BF_Bin, y = Count, fill = Class), stat = "identity", position = "dodge")
sst_mean_mean_plot_annotated <- sst_mean_bfs_mean_plot + theme(legend.position = "top") + 
  annotate("text", x = 5, y = c(4000, 3650, 3300), label = c("p < 0.001", "Odds Ratio = 3.239", "95% CI = 2.911, 3.609")) + 
  geom_vline(xintercept = 4.5) + ggtitle("SST mean") + 
  scale_x_discrete(breaks = c("bf1", "bf2", "bf3", "bf4", "bf5", "bf6"), 
                   labels = c("BF<0", "0<BF<5", "5<BF<10", "10<BF<15", "15<BF<20", "BF>20"))

#SST min histogram of means
sst_min_bfs_mean_plot <- ggplot(data = SST_min_bfs_mean) + 
  geom_bar(mapping = aes(x = BF_Bin, y = Count, fill = Class), stat = "identity", position = "dodge")
sst_min_mean_plot_annotated <- sst_min_bfs_mean_plot + theme(legend.position = "top") + 
  annotate("text", x = 5, y = c(4000, 3650, 3300), label = c("p < 0.001", "Odds Ratio = 4.147", "95% CI = (3.734, 4.610)")) + 
  geom_vline(xintercept = 4.5) + ggtitle("SST min") + 
  scale_x_discrete(breaks = c("bf1", "bf2", "bf3", "bf4", "bf5", "bf6"), 
                   labels = c("BF<0", "0<BF<5", "5<BF<10", "10<BF<15", "15<BF<20", "BF>20"))

#SST max histogram of means
sst_max_bfs_mean_plot <- ggplot(data = SST_max_bfs_mean) + 
  geom_bar(mapping = aes(x= BF_Bin, y = Count, fill = Class), stat = "identity", position = "dodge")
sst_max_mean_plot_annotated <- sst_max_bfs_mean_plot + theme(legend.position = "top") + 
  annotate("text", x = 5, y = c(4000, 3650, 3300), label = c("p < 0.001", "Odds Ratio = 0.120", "95% CI = (0.135, 0.288)")) + 
  geom_vline(xintercept = 4.5) + ggtitle("SST max") + 
  scale_x_discrete(breaks = c("bf1", "bf2", "bf3", "bf4", "bf5", "bf6"), 
                   labels = c("BF<0", "0<BF<5", "5<BF<10", "10<BF<15", "15<BF<20", "BF>20"))

#Lat histogram of means
lat_bfs_mean_plot <- ggplot(data = lat_bfs_mean) + 
  geom_bar(mapping = aes(x= BF_Bin, y = Count, fill = Class), stat = "identity", position = "dodge")
lat_mean_plot_annotated <- lat_bfs_mean_plot + theme(legend.position = "top") + 
  annotate("text", x = 5, y = c(4000, 3650, 3300), label = c("p < 0.01", "Odds Ratio = 3.112", "95% CI = (2.812, 3.446)")) + 
  geom_vline(xintercept = 3.5) + ggtitle("Latitude") + 
  scale_x_discrete(breaks = c("bf1", "bf2", "bf3", "bf4", "bf5", "bf6"), 
                   labels = c("BF<0", "0<BF<5", "5<BF<10", "10<BF<15", "15<BF<20", "BF>20"))

#arrange all histograms on same page
ggarrange(sss_mean_mean_plot_annotated, sst_mean_mean_plot_annotated, sst_min_mean_plot_annotated, sst_max_mean_plot_annotated, lat_mean_plot_annotated, ncol = 2, nrow =3)
ggarrange(sst_mean_mean_plot_annotated, sst_min_mean_plot_annotated, sst_max_mean_plot_annotated, ncol = 3, nrow = 1)

######## Histograms with median data ########

#create histograms (1 histogram = 1 environmental covariate with median for sim & real models) --> 5 histograms total

#sort dataframe by environmental covariate
SSS_mean_bfs_median <- bfs_median_all[1:12, -2]
SST_mean_bfs_median <- bfs_median_all[13:24, -2]
SST_min_bfs_median <- bfs_median_all[25:36, -2]
SST_max_bfs_median <- bfs_median_all[37:48, -2]
lat_bfs_median <- bfs_median_all[49:60, -2]

#SSS mean histogram of medians
sss_mean_bfs_median_plot <- ggplot(data = SSS_mean_bfs_median) + 
  geom_bar(mapping = aes(x= BF_Bin, y = Count, fill = Class), stat = "identity", position = "dodge")
sss_mean_median_plot_annotated <- sss_mean_bfs_median_plot + theme(legend.position = "top") + 
  annotate("text", x = 5, y = c(4000, 3650, 3300), label = c("p < 0.001", "Odds Ratio = 2.335", "95% CI = (2.111, 2.584)")) + 
  geom_vline(xintercept = 3.5) + ggtitle("SSS mean") + 
  scale_x_discrete(breaks = c("bf1", "bf2", "bf3", "bf4", "bf5", "bf6"), 
                   labels = c("BF<0", "0<BF<5", "5<BF<10", "10<BF<15", "15<BF<20", "BF>20"))

#SST mean histogram of medians
sst_mean_bfs_median_plot <- ggplot(data = SST_mean_bfs_median) + 
  geom_bar(mapping = aes(x= BF_Bin, y = Count, fill = Class), stat = "identity", position = "dodge")
sst_mean_median_plot_annotated <- sst_mean_bfs_median_plot + theme(legend.position = "top") + 
  annotate("text", x = 5, y = c(4000, 3750, 3500), label = c("p < 0.001", "Odds Ratio = 3.439", "95% CI = (3.084, 3.838)")) + 
  geom_vline(xintercept = 3.5) + ggtitle("SST mean") + 
  scale_x_discrete(breaks = c("bf1", "bf2", "bf3", "bf4", "bf5", "bf6"), 
                   labels = c("BF<0", "0<BF<5", "5<BF<10", "10<BF<15", "15<BF<20", "BF>20"))

#SST min histogram of medians
sst_min_bfs_median_plot <- ggplot(data = SST_min_bfs_median) + 
  geom_bar(mapping = aes(x= BF_Bin, y = Count, fill = Class), stat = "identity", position = "dodge")
sst_min_median_plot_annotated <- sst_min_bfs_median_plot + theme(legend.position = "top") + 
  annotate("text", x = 5, y = c(4000, 3750, 3500), label = c("p < 0.001", "Odds Ratio = 3.551", "95% CI = (3.212, 3.929)")) + 
  geom_vline(xintercept = 3.5) + ggtitle("SST min") + 
  scale_x_discrete(breaks = c("bf1", "bf2", "bf3", "bf4", "bf5", "bf6"), 
                   labels = c("BF<0", "0<BF<5", "5<BF<10", "10<BF<15", "15<BF<20", "BF>20"))

#SST max histogram of medians
sst_max_bfs_median_plot <- ggplot(data = SST_max_bfs_median) + 
  geom_bar(mapping = aes(x= BF_Bin, y = Count, fill = Class), stat = "identity", position = "dodge")
sst_max_median_plot_annotated <- sst_max_bfs_median_plot + theme(legend.position = "top") + 
  annotate("text", x = 5, y = c(4000, 3750, 3500), label = c("p < 0.001", "Odds Ratio = 2.258", "95% CI = (1.220, 4.363)")) + 
  geom_vline(xintercept = 3.5) + ggtitle("SST max") + 
  scale_x_discrete(breaks = c("bf1", "bf2", "bf3", "bf4", "bf5", "bf6"), 
                   labels = c("BF<0", "0<BF<5", "5<BF<10", "10<BF<15", "15<BF<20", "BF>20"))

#Lat histogram of medians
lat_bfs_median_plot <- ggplot(data = lat_bfs_median) + 
  geom_bar(mapping = aes(x= BF_Bin, y = Count, fill = Class), stat = "identity", position = "dodge")
lat_median_plot_annotated <- lat_bfs_median_plot + theme(legend.position = "top") + 
  annotate("text", x = 5, y = c(4000, 3650, 3300), label = c("p < 0.001", "Odds Ratio = 3.541", "95% CI = (3.189, 3.935)")) + 
  geom_vline(xintercept = 3.5) + ggtitle("Latitude") + 
  scale_x_discrete(breaks = c("bf1", "bf2", "bf3", "bf4", "bf5", "bf6"), 
                   labels = c("BF<0", "0<BF<5", "5<BF<10", "10<BF<15", "15<BF<20", "BF>20"))

#arrange all histograms on same page
ggarrange(sss_mean_median_plot_annotated, sst_mean_median_plot_annotated, sst_min_median_plot_annotated, sst_max_median_plot_annotated, lat_median_plot_annotated, ncol = 2, nrow =3)

###########################################################################################################################################

######## ECDFs & K-S Tests ########

#using average BFs from all 10 runs (not median) per locus

#SSS_mean analysis
#subset data for plot
SSS_mean_df <- BFs_KS[which(BFs_KS$cov == "SSS_mean"), ]

#subset data for K-S test
SSS_mean_permuted <- SSS_mean_df[which(SSS_mean_df$dataset == "permuted"), ] 
SSS_mean_permuted <- SSS_mean_permuted$BFs

SSS_mean_real <- SSS_mean_df[which(SSS_mean_df$dataset == "real"), ] 
SSS_mean_real <- SSS_mean_real$BFs

#run K-S test for SSS mean dataset
SSS_mean_kstest <- ks.test(SSS_mean_permuted, SSS_mean_real)

#plot distributions
SSS_ecdf_plot <- ggplot(data = SSS_mean_df, aes(x = BFs, color = dataset)) + stat_ecdf(geom = "step", size = 1.5) + 
  geom_hline(yintercept = 1, size = 1, linetype = "dashed", color = "#666666") + 
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", color = "#666666") + 
  annotate("text", x = 40, y = 0.65, label = "D = 0.539", size = 12) + 
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
SST_mean_permuted <- SST_mean_permuted$BFs

SST_mean_real <- SST_mean_df[which(SST_mean_df$dataset == "real"), ] 
SST_mean_real <- SST_mean_real$BFs

#run K-S test for SST mean dataset
SST_mean_kstest <- ks.test(SST_mean_permuted, SST_mean_real)

#plot distributions
SST_mean_ecdf_plot <- ggplot(data = SST_mean_df, aes(x = BFs, color = dataset)) + stat_ecdf(geom = "step", size = 1.5) + 
  geom_hline(yintercept = 1, size = 1, linetype = "dashed", color = "#666666") + 
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", color = "#666666") + 
  annotate("text", x = 40, y = 0.65, label = "D = 0.812", size = 12) + 
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
SST_min_permuted <- SST_min_permuted$BFs

SST_min_real <- SST_min_df[which(SST_min_df$dataset == "real"), ] 
SST_min_real <- SST_min_real$BFs

#run K-S test for SST min dataset
SST_min_kstest <- ks.test(SST_min_permuted, SST_min_real)

#plot distributions
SST_min_ecdf_plot <- ggplot(data = SST_min_df, aes(x = BFs, color = dataset)) + stat_ecdf(geom = "step", size = 1.5) + 
  geom_hline(yintercept = 1, size = 1, linetype = "dashed", color = "#666666") + 
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", color = "#666666") + 
  annotate("text", x = 40, y = 0.65, label = "D = 0.871", size = 12) + 
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
SST_max_permuted <- SST_max_permuted$BFs

SST_max_real <- SST_max_df[which(SST_max_df$dataset == "real"), ] 
SST_max_real <- SST_max_real$BFs

#run K-S test for SST max dataset
SST_max_kstest <- ks.test(SST_max_permuted, SST_max_real)

#plot distributions
SST_max_ecdf_plot <- ggplot(data = SST_max_df, aes(x = BFs, color = dataset)) + stat_ecdf(geom = "step", size = 1.5) + 
  geom_hline(yintercept = 1, size = 1, linetype = "dashed", color = "#666666") + 
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", color = "#666666") + 
  annotate("text", x = 40, y = 0.65, label = "D = 0.657", size = 12) + 
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
Lat_permuted <- Lat_permuted$BFs

Lat_real <- Lat_df[which(Lat_df$dataset == "real"), ] 
Lat_real <- Lat_real$BFs

#run K-S test for Lat dataset
Lat_kstest <- ks.test(Lat_permuted, Lat_real)

#plot distributions
Lat_ecdf_plot <- ggplot(data = Lat_df, aes(x = BFs, color = dataset)) + stat_ecdf(geom = "step", size = 1.5) + 
  geom_hline(yintercept = 1, size = 1, linetype = "dashed", color = "#666666") + 
  geom_hline(yintercept = 0, size = 1, linetype = "dashed", color = "#666666") + 
  annotate("text", x = 40, y = 0.65, label = "D = 0.761", size = 12) + 
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
