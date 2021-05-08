################################ Script for Creating Tajima's D v. Pi Plots  ########################################################

#Uses Tajima's D estimates from VCFtools --> only synonymous sites incorporated into calculations (ID'd from mapping to A. frenatus)
#Uses pi estimates (across contig) from VCFtools --> only synonymous sites as well

#################################################################################################################################################

######## Set-up ########

#set working directory
setwd("C:/Users/Rene/Dropbox/Pinsky_Lab/Transcriptome_Proj/R_scripts/A_clarkii_transcriptomics/")
getwd()

remove(list = ls())

#load libraries
library(tidyverse)

#read in data
TD_only_all <- read.csv("Data/SYN_mac1_TajimasD_combined_full.csv", header = TRUE, row.names = 1)
TD_outlier_only_all <- read.csv("Data/SYN_mac1_TajimasD_outlier_combined_full.csv", header = TRUE, row.names = 1)
outlierseq <- read.csv("Data/outlier_sequences_TD_mac1_SYN.csv", header = TRUE) #only if running separately to visualize data
pi_only_all <- read.csv("Data/pi_combined_site_full_mac1_SYN_aggregate.csv", header = TRUE, row.names = 1)

data_all <- cbind(TD_only_all, pi_only_all$pi_avg)
  colnames(data_all) <- c("CHROM", "BIN_START", "N_SNPS", "TajimaD", "uniqseq", "Pop", "NUM", "pi")

#################################################################################################################################################
  
######## Format data ########

#add outlier status
outliers <- data_all[data_all$uniqseq %in% outlierseq$CONTIG, ] #df with only outliers
nonoutliers <- data_all[!(data_all$uniqseq %in% outlierseq$CONTIG), ] #df with all other loci

#set seq status
outliers$Status <- c(rep("Outlier", times = 96))
nonoutliers$Status <- c(rep("Not_Outlier", times = 1896))
data_all$Status <- c(rep("All", times = 1992))

#merge together
data_all <- rbind(outliers, nonoutliers, data_all)
data_all$Status <- factor(data_all$Status, levels = c("All", "Not_Outlier", "Outlier"))

#subset to sampling sites
all_only <- subset(data_all, Pop == "All" & Status != "All")
J_only <- subset(data_all, Pop == "Japan" & Status != "All")
N_only <- subset(data_all, Pop == "Indonesia" & Status != "All")
P_only <- subset(data_all, Pop == "Philippines" & Status != "All")

#subset to outliers only
outliers_only <- subset(data_all, Status == "Outlier")

#################################################################################################################################################

######## Visualize data ########
#Scatterplots of Tajima's D v. pi

#TD v pi for pooled individuals
all_plot <- ggplot(data = all_only %>% arrange(Status), aes(x = pi, y = TajimaD, color = Status)) + 
  geom_point(size = 8) + geom_hline(yintercept = 0, color = "black", size = 3)
all_plot_annotated <- all_plot + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(size = 2), 
        axis.ticks = element_line(color = "black", size = 2), 
        axis.text = element_text(size = 34, color = "black"), 
        axis.title = element_text(size = 34, face = "bold"), legend.position = "top", 
        legend.text = element_text(size = 32), legend.title = element_text(size = 32))
all_plot_annotated

#TD v pi for Japan
J_plot <- ggplot(data = J_only %>% arrange(Status), aes(x = pi, y = TajimaD, color = Status)) + 
  geom_point(size = 8) + geom_hline(yintercept = 0, color = "black", size = 3)
J_plot_annotated <- J_plot + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(size = 2), 
        axis.ticks = element_line(color = "black", size = 2), 
        axis.text = element_text(size = 34, color = "black"), 
        axis.title = element_text(size = 34, face = "bold"), legend.position = "top", 
        legend.text = element_text(size = 32), legend.title = element_text(size = 32))
J_plot_annotated

#TD v pi for Philippines
P_plot <- ggplot(data = P_only %>% arrange(Status), aes(x = pi, y = TajimaD, color = Status)) + 
  geom_point(size = 8) + geom_hline(yintercept = 0, color = "black", size = 3)
P_plot_annotated <- P_plot + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(size = 2), 
        axis.ticks = element_line(color = "black", size = 2), 
        axis.text = element_text(size = 34, color = "black"), 
        axis.title = element_text(size = 34, face = "bold"), legend.position = "top", 
        legend.text = element_text(size = 32), legend.title = element_text(size = 32))
P_plot_annotated

#TD v pi for Indonesia
N_plot <- ggplot(data = N_only %>% arrange(Status), aes(x = pi, y = TajimaD, color = Status)) + 
  geom_point(size = 8) + geom_hline(yintercept = 0, color = "black", size = 3)
N_plot_annotated <- N_plot + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(size = 2), 
        axis.ticks = element_line(color = "black", size = 2), 
        axis.text = element_text(size = 34, color = "black"), 
        axis.title = element_text(size = 34, face = "bold"), legend.position = "top", 
        legend.text = element_text(size = 32), legend.title = element_text(size = 32))
N_plot_annotated

#TD v pi with only outliers
outliers_plot <- ggplot(data = outliers_only, aes(x = pi, y = TajimaD, color = Pop)) + 
  geom_point(size = 8) + geom_hline(yintercept = 0, color = "black", size = 3)
outliers_plot_annotated <- outliers_plot + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(size = 2), 
        axis.ticks = element_line(color = "black", size = 2), 
        axis.text = element_text(size = 34, color = "black"), 
        axis.title = element_text(size = 34, face = "bold"), legend.position = "top", 
        legend.text = element_text(size = 32), legend.title = element_text(size = 32))
outliers_plot_annotated