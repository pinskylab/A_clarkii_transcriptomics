############################################# Script for creating allele frequency line plots ###################################################

#allele frequencies at each location for different sets of loci (loci associated with diff SST env variables)
#includes plot of mean allele frequencies for all 4212 loci
#allele frequencies are polarized so J is highest

#################################################################################################################################################

######## Set-up ########

#set working directory
setwd("C:/Users/rclar/Dropbox/Pinsky_Lab/Transcriptome_Proj/R_scripts/A_clarkii_transcriptomics/")
getwd()

remove(list = ls())

#loading libraries
library(readr) #v.2.1.4
library(tidyverse) #v.2.0.0
library(gridExtra) #v.2.3

#read in data
allele_freqs <- read_csv("Data/polarized_allele_freqs.csv", col_names = TRUE)
  colnames(allele_freqs) <- c("Number", "Contig_bp", "p", "Status", "Pop")

#################################################################################################################################################

######## Allele frequency plots ########

allele_freqs$Pop2 <- factor(allele_freqs$Pop, levels = c("Japan", "Philippines", "Indonesia")) #ordering X-axis
allele_freqs$Status <- factor(allele_freqs$Status, c("Outliers-Mean", "NonOutliers-Mean", "Outlier", "Non-Outlier"))

#no outlier loci AF plot
allele_freqs_NO_only <- subset(allele_freqs, Status == "Non-Outlier" | Status == "NonOutliers-Mean")
allele_freqs_NO_plot <- ggplot(data = allele_freqs_NO_only, aes(x = Pop2, y = p, group = Contig_bp, color = Status, size = Status)) + 
  geom_line(data = subset(allele_freqs_NO_only, Status == "Non-Outlier"), alpha = 0.2) + geom_point(data = subset(allele_freqs_NO_only, Status == "Non-Outlier"), alpha = 0.2) + 
  geom_line(data = subset(allele_freqs_NO_only, Status == "NonOutliers-Mean")) + geom_point(data = subset(allele_freqs_NO_only, Status == "NonOutliers-Mean"))
allele_freqs_NO_plot_annotated <- allele_freqs_NO_plot + theme_bw() + 
  labs(title = "Allele Frequencies", y = "p", x = "Sampling Location") + 
  scale_size_manual(values = c(1.5, 4)) + scale_color_manual(values = c("#CCCCCC", "#333333")) + 
  annotate("text", x = 0.5, y = 1, label = "B", size = 16) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(size = 2), 
        axis.ticks = element_line(color = "black", size = 2), 
        axis.text = element_text(size = 28, color = "black"), 
        axis.title = element_text(size = 30), legend.position = "top",
        plot.title = element_blank(), plot.margin = unit(c(.5,.5,.5,.5), "cm"),
        legend.text = element_text(size = 30), legend.title = element_text(size = 30))
allele_freqs_NO_plot_annotated

#outlier only AF plot
allele_freqs_outlier_only <- subset(allele_freqs, Status == "Outlier" | Status == "Outliers-Mean")
allele_freqs_outlier_plot <- ggplot(data = allele_freqs_outlier_only, aes(x = Pop2, y = p, group = Contig_bp, color = Status, size = Status)) + 
  geom_line(data = subset(allele_freqs_outlier_only, Status == "Outlier"), alpha = 0.2) + 
  geom_line(data = subset(allele_freqs_outlier_only, Status == "Outliers-Mean"))
allele_freqs_outlier_annotated <- allele_freqs_outlier_plot + theme_bw() + 
  labs(title = "Allele Frequencies", y = "p", x = "Sampling Location") + 
  scale_size_manual(values = c(1.5, 4)) + scale_color_manual(values = c("#3333CC", "#000033")) + 
  annotate("text", x = 0.5, y = 1, label = "A", size = 16) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(size = 2), 
        axis.ticks = element_line(color = "black", size = 2), 
        axis.text = element_text(size = 28, color = "black"), 
        axis.title = element_text(size = 30), legend.position = "top",
        plot.title = element_blank(), plot.margin = unit(c(.5,.5,.5,.5), "cm"),
        legend.text = element_text(size = 30), legend.title = element_text(size = 30))
allele_freqs_outlier_annotated

#all loci AF plot
allele_freqs_plot <- ggplot(data = allele_freqs, aes(x = Pop2, y = p, group = Contig_bp, color = Status, size = Status)) + 
  geom_line(data = subset(allele_freqs, Status == "Non-Outlier")) + geom_point(data = subset(allele_freqs, Status == "Non-Outlier")) + 
  geom_line(data = subset(allele_freqs, Status == "Outlier")) + geom_point(data = subset(allele_freqs, Status == "Outlier")) + 
  geom_line(data = subset(allele_freqs, Status == "NonOutliers-Mean")) + geom_point(data = subset(allele_freqs, Status == "NonOutliers-Mean")) + 
  geom_line(data = subset(allele_freqs, Status == "Outliers-Mean")) + geom_point(data = subset(allele_freqs, Status == "Outliers-Mean")) 
allele_freqs_plot_annotated <- allele_freqs_plot + theme_bw() + 
  labs(title = "Allele Frequencies", y = "p", x = "Sampling Location") + 
  scale_size_manual(values = c(4, 0.5, 0.5, 4)) + scale_color_manual(values = c("#333333", "#CCCCCC", "#3333FF", "#0000CC")) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(size = 2), 
        axis.ticks = element_line(color = "black", size = 2), 
        axis.text = element_text(size = 28, color = "black"), 
        axis.title = element_text(size = 30), legend.position = "top", 
        legend.text = element_text(size = 30), legend.title = element_text(size = 30))
allele_freqs_plot_annotated

#arrange plots together
AF_all_plot <- grid.arrange(allele_freqs_outlier_annotated, allele_freqs_NO_plot_annotated, ncol = 1)
AF_all_plot
