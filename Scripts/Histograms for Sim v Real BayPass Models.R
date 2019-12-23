############################################## Create Histograms for Sim v Real BayPass Output #################################################

#histograms for simulation vs real BayPass output by environmental variable
#includes Fisher's exact test for mean and median BF for each loci (averaging across 10 simulation v 10 real runs)

################################################################################################################################################

######## Set-up ########

remove(list = ls())

#Set working directory
setwd("C:/Users/Rene/Dropbox/Pinsky_Lab/Transcriptome_Proj/R_scripts/A_clarkii_transcriptomics/")
getwd()

#load libraries
library(tidyverse)
library(ggpubr)

#read in data
bfs_mean_all <- read_csv("Data/sim_v_real_mean.csv", col_names = TRUE) #mean BF across 10 runs (sim and real)
bfs_median_all <- read_csv("Data/sim_v_real_median.csv", col_names = TRUE) #median BF across 10 runs (sim and real)

################################################################################################################################################

######## Fisher's exact test with mean data ########

#create 2x2 contingency tables for mean tests
SSS_mean_mean_table <- matrix(c(5636.8, 92.2, 5729, 0), nrow = 2, ncol = 2) #create matrix with: col 1 = simulated data, col 2 = real data, row 1 = # loci w/BF<10, row 2 = # loci w/BF>10
SST_mean_mean_table <- matrix(c(5492.6, 236.4, 4740, 989), nrow = 2, ncol = 2)
SST_min_mean_table <- matrix(c(5545.4, 183.6, 5054, 675), nrow = 2, ncol = 2)
SST_max_mean_table <- matrix(c(3593.7, 2135.3, 2973, 2756), nrow = 2, ncol = 2)
Lat_mean_table <- matrix(c(5697.4, 31.6, 5729, 0), nrow = 2, ncol = 2)

#Fisher's exact test by env. covariable
FT_SSS_mean_mean <- fisher.test(SSS_mean_mean_table) #conduct Fisher's Exact Test rounding to nearest integer and print to consule
FT_SST_mean_mean <- fisher.test(SST_mean_mean_table)
FT_SST_min_mean <- fisher.test(SST_min_mean_table)
FT_SST_max_mean <- fisher.test(SST_max_mean_table)
FT_lat_mean <- fisher.test(Lat_mean_table)

######## Fisher's exact test with median data ########

#create 2x2 contingency tables for median tests
SSS_mean_median_table <- matrix(c(5682, 47, 5729, 0), nrow = 2, ncol = 2)
SST_mean_median_table <- matrix(c(5433.5, 295.5, 4740, 989), nrow = 2, ncol = 2)
SST_min_median_table <- matrix(c(5560, 169, 5054, 675), nrow = 2, ncol = 2)
SST_max_median_table <- matrix(c(3476, 2253, 2973, 2756), nrow = 2, ncol = 2)
Lat_median_table <- matrix(c(5729, 0, 5729, 0), nrow = 2, ncol = 2)

#Fisher's exact test by env. covariable
FT_SSS_mean_median <- fisher.test(SSS_mean_median_table)
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
  annotate("text", x = 5, y = c(4000, 3650, 3300), label = c("p < 0.01", "Odds Ratio = 0.0", "95% CI = (0.0, 0.0403)")) + 
  geom_vline(xintercept = 3.5) + ggtitle("SSS mean") + 
  scale_x_discrete(breaks = c("bf1", "bf2", "bf3", "bf4", "bf5", "bf6"), 
                   labels = c("BF<0", "0<BF<5", "5<BF<10", "10<BF<15", "15<BF<20", "BF>20")) #relabel X-axis, move legend to top of chart, add p-value and odds ratio, add vertical line dividing two categories and add title

#SST mean histogram of means
sst_mean_bfs_mean_plot <- ggplot(data = SST_mean_bfs_mean) + 
  geom_bar(mapping = aes(x= BF_Bin, y = Count, fill = Class), stat = "identity", position = "dodge")
sst_mean_mean_plot_annotated <- sst_mean_bfs_mean_plot + theme(legend.position = "top") + 
  annotate("text", x = 5, y = c(4000, 3750, 3500), label = c("p < 0.01", "Odds Ratio = 4.8558", "95% CI = (4.185, 5.651)")) + 
  geom_vline(xintercept = 3.5) + ggtitle("SST mean") + 
  scale_x_discrete(breaks = c("bf1", "bf2", "bf3", "bf4", "bf5", "bf6"), 
                   labels = c("BF<0", "0<BF<5", "5<BF<10", "10<BF<15", "15<BF<20", "BF>20"))

#SST min histogram of means
sst_min_bfs_mean_plot <- ggplot(data = SST_min_bfs_mean) + 
  geom_bar(mapping = aes(x = BF_Bin, y = Count, fill = Class), stat = "identity", position = "dodge")
sst_min_mean_plot_annotated <- sst_min_bfs_mean_plot + theme(legend.position = "top") + 
  annotate("text", x = 5, y = c(4000, 3750, 3500), label = c("p < 0.01", "Odds Ratio = 4.0225", "95% CI = (3.3982, 4.7846)")) + 
  geom_vline(xintercept = 3.5) + ggtitle("SST min") + 
  scale_x_discrete(breaks = c("bf1", "bf2", "bf3", "bf4", "bf5", "bf6"), 
                   labels = c("BF<0", "0<BF<5", "5<BF<10", "10<BF<15", "15<BF<20", "BF>20"))

#SST max histogram of means
sst_max_bfs_mean_plot <- ggplot(data = SST_max_bfs_mean) + 
  geom_bar(mapping = aes(x= BF_Bin, y = Count, fill = Class), stat = "identity", position = "dodge")
sst_max_mean_plot_annotated <- sst_max_bfs_mean_plot + theme(legend.position = "top") + 
  annotate("text", x = 5, y = c(4000, 3750, 3500), label = c("p < 0.01", "Odds Ratio = 1.5604", "95% CI = (1.4474, 1.6824)")) + 
  geom_vline(xintercept = 3.5) + ggtitle("SST max") + 
  scale_x_discrete(breaks = c("bf1", "bf2", "bf3", "bf4", "bf5", "bf6"), 
                   labels = c("BF<0", "0<BF<5", "5<BF<10", "10<BF<15", "15<BF<20", "BF>20"))

#Lat histogram of means
lat_bfs_mean_plot <- ggplot(data = lat_bfs_mean) + 
  geom_bar(mapping = aes(x= BF_Bin, y = Count, fill = Class), stat = "identity", position = "dodge")
lat_mean_plot_annotated <- lat_bfs_mean_plot + theme(legend.position = "top") + 
  annotate("text", x = 5, y = c(4000, 3650, 3300), label = c("p < 0.01", "Odds Ratio = 0", "95% CI = (0.0, 0.1216)")) + 
  geom_vline(xintercept = 3.5) + ggtitle("Latitude") + 
  scale_x_discrete(breaks = c("bf1", "bf2", "bf3", "bf4", "bf5", "bf6"), 
                   labels = c("BF<0", "0<BF<5", "5<BF<10", "10<BF<15", "15<BF<20", "BF>20"))

#arrange all histograms on same page
ggarrange(sss_mean_mean_plot_annotated, sst_mean_mean_plot_annotated, sst_min_mean_plot_annotated, sst_max_mean_plot_annotated, lat_mean_plot_annotated, ncol = 2, nrow =3)

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
  annotate("text", x = 5, y = c(4000, 3650, 3300), label = c("p < 0.01", "Odds Ratio = 0.0", "95% CI = (0.0, 0.0810)")) + 
  geom_vline(xintercept = 3.5) + ggtitle("SSS mean") + 
  scale_x_discrete(breaks = c("bf1", "bf2", "bf3", "bf4", "bf5", "bf6"), 
                   labels = c("BF<0", "0<BF<5", "5<BF<10", "10<BF<15", "15<BF<20", "BF>20"))

#SST mean histogram of medians
sst_mean_bfs_median_plot <- ggplot(data = SST_mean_bfs_median) + 
  geom_bar(mapping = aes(x= BF_Bin, y = Count, fill = Class), stat = "identity", position = "dodge")
sst_mean_median_plot_annotated <- sst_mean_bfs_median_plot + theme(legend.position = "top") + 
  annotate("text", x = 5, y = c(4000, 3750, 3500), label = c("p < 0.01", "Odds Ratio = 3.83", "95% CI = (3.3402, 4.4019)")) + 
  geom_vline(xintercept = 3.5) + ggtitle("SST mean") + 
  scale_x_discrete(breaks = c("bf1", "bf2", "bf3", "bf4", "bf5", "bf6"), 
                   labels = c("BF<0", "0<BF<5", "5<BF<10", "10<BF<15", "15<BF<20", "BF>20"))

#SST min histogram of medians
sst_min_bfs_median_plot <- ggplot(data = SST_min_bfs_median) + 
  geom_bar(mapping = aes(x= BF_Bin, y = Count, fill = Class), stat = "identity", position = "dodge")
sst_min_median_plot_annotated <- sst_min_bfs_median_plot + theme(legend.position = "top") + 
  annotate("text", x = 5, y = c(4000, 3750, 3500), label = c("p < 0.01", "Odds Ratio = 4.3936", "95% CI = (3.6898, 5.2538)")) + 
  geom_vline(xintercept = 3.5) + ggtitle("SST min") + 
  scale_x_discrete(breaks = c("bf1", "bf2", "bf3", "bf4", "bf5", "bf6"), 
                   labels = c("BF<0", "0<BF<5", "5<BF<10", "10<BF<15", "15<BF<20", "BF>20"))

#SST max histogram of medians
sst_max_bfs_median_plot <- ggplot(data = SST_max_bfs_median) + 
  geom_bar(mapping = aes(x= BF_Bin, y = Count, fill = Class), stat = "identity", position = "dodge")
sst_max_median_plot_annotated <- sst_max_bfs_median_plot + theme(legend.position = "top") + 
  annotate("text", x = 5, y = c(4000, 3750, 3500), label = c("p < 0.01", "Odds Ratio = 1.4301", "95% CI = (1.3271, 1.5414)")) + 
  geom_vline(xintercept = 3.5) + ggtitle("SST max") + 
  scale_x_discrete(breaks = c("bf1", "bf2", "bf3", "bf4", "bf5", "bf6"), 
                   labels = c("BF<0", "0<BF<5", "5<BF<10", "10<BF<15", "15<BF<20", "BF>20"))

#Lat histogram of medians
lat_bfs_median_plot <- ggplot(data = lat_bfs_median) + 
  geom_bar(mapping = aes(x= BF_Bin, y = Count, fill = Class), stat = "identity", position = "dodge")
lat_median_plot_annotated <- lat_bfs_median_plot + theme(legend.position = "top") + 
  annotate("text", x = 5, y = c(4000, 3650, 3300), label = c("p = 1", "Odds Ratio = 0", "95% CI = (0, Infinity)")) + 
  geom_vline(xintercept = 3.5) + ggtitle("Latitude") + 
  scale_x_discrete(breaks = c("bf1", "bf2", "bf3", "bf4", "bf5", "bf6"), 
                   labels = c("BF<0", "0<BF<5", "5<BF<10", "10<BF<15", "15<BF<20", "BF>20"))

#arrange all histograms on same page
ggarrange(sss_mean_median_plot_annotated, sst_mean_median_plot_annotated, sst_min_median_plot_annotated, sst_max_median_plot_annotated, lat_median_plot_annotated, ncol = 2, nrow =3)
