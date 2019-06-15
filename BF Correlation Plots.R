############################################# Script for creating Bayes Factors Correlation Plots #################################################

#includes BFs v. BFs plots for different environmental factors (SST variables)
#includes BFs. v. XtX plots for SST variables

#################################################################################################################################################

######## Set-up ########

#set working directory
setwd("C:/Users/Rene/Dropbox/Pinsky_Lab/Transcriptome_Proj/R_scripts/A_clarkii_transcriptomics/")
getwd()

remove(list = ls())

#load libraries
library(readr)
library(tidyverse)

#read in data
BFs <- read_csv("BF_total.csv", col_names = TRUE)

################################################################################################################################################

######## BFs scatter plots ########

#correlation statistics
SSTMean_SSTMin_r <- cor(BFs$SST_Mean, BFs$SST_Min)
SSTMean_SSTMax_r <- cor(BFs$SST_Max, BFs$SST_Mean)
SSTMax_SSTMin_r <- cor(BFs$SST_Max, BFs$SST_Min)

#SSTmean v. SSTmin plot
SSTmean_SSTmin_plot <- ggplot(data = BFs, aes(x = SST_Mean, y = SST_Min)) + geom_point(aes(size = 0.5)) + 
  geom_hline(yintercept = 15, size = 2, linetype = "dashed", color = "#666666") + 
  geom_vline(xintercept = 15, size = 2, linetype = "dashed", color = "#666666") + 
  geom_abline(intercept = 0, slope = 1, size = 2, color = "#999999") + 
  annotate("text", x = 3, y = 20, label = "r = 0.873", size = 10)
SSTmean_SSTmin_annotated_plot <- SSTmean_SSTmin_plot + ggtitle("SST Mean BF v. SST Min BF") + labs(x = "SST Mean BF", y = "SST Min BF") + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.line = element_line(size = 1), plot.title = element_text(size = 14, face = "bold"), 
                     axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 30, color = "black"), 
                     axis.title = element_text(size = 18, face = "bold")) + scale_y_continuous(limits = c(0, 25)) + scale_x_continuous(limits = c(0, 25))
SSTmean_SSTmin_annotated_plot

#SSTmean v. SSTmax plot
SSTmean_SSTmax_plot <- ggplot(data = BFs, aes(x = SST_Mean, y = SST_Max)) + geom_point(aes(size = 0.5)) + 
  geom_hline(yintercept = 15, size = 2, linetype = "dashed", color = "#666666") + 
  geom_vline(xintercept = 15, size = 2, linetype = "dashed", color = "#666666") + 
  geom_abline(intercept = 0, slope = 1, size = 2, color = "#999999") + 
  annotate("text", x = 20, y = 7, label = "r = 0.215", size = 10)
SSTmean_SSTmax_annotated_plot <- SSTmean_SSTmax_plot + ggtitle("SST Mean BF v. SST Max BF") + labs(x = "SST Mean BF", y = "SST Max BF") + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.line = element_line(size = 1), plot.title = element_text(size = 14, face = "bold"), 
                     axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 30, color = "black"), 
                     axis.title = element_text(size = 18, face = "bold")) + scale_y_continuous(limits = c(0, 25)) + scale_x_continuous(limits = c(0, 25))
SSTmean_SSTmax_annotated_plot

#SSTmin v. SSTmax plot
SSTmin_SSTmax_plot <- ggplot(data = BFs, aes(x = SST_Min, y = SST_Max)) + geom_point(aes(size = 0.5)) + 
  geom_hline(yintercept = 15, size = 1.5, linetype = "dashed", color = "#666666") + 
  geom_vline(xintercept = 15, size = 1.5, linetype = "dashed", color = "#666666") + 
  geom_abline(intercept = 0, slope = 1, size = 1.5, color = "#999999") + 
  annotate("text", x = 20, y = 7, label = "r = 0.293", size = 4)
SSTmin_SSTmax_annotated_plot <- SSTmin_SSTmax_plot + ggtitle("SST Min BF v. SST Max BF") + labs(x = "SST Min BF", y = "SST Max BF") + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.line = element_line(size = 1), plot.title = element_text(size = 14, face = "bold"), 
                     axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 14, color = "black"), 
                     axis.title = element_text(size = 18, face = "bold")) + scale_y_continuous(limits = c(0, 25)) + scale_x_continuous(limits = c(0, 25))
SSTmin_SSTmax_annotated_plot

###############################################################################################################################################

######## BF v. XtX scatter plots ########

#SSTmean v. XtX plot
SSTmean_XtX_plot <- ggplot(data = BFs, aes(x = SST_Mean, y = M_XtX)) + geom_point(aes(size = 0.5)) + 
  geom_hline(yintercept = 5.78, size = 2, color = "#666666", linetype = "dashed") + 
  geom_vline(xintercept = 15, size = 2, color = "#666666", linetype = "dashed")
SSTmean_XtX_annotated_plot <- SSTmean_XtX_plot + ggtitle("SST Mean BF v. XtX") + labs(x = "SST Mean BF", y = "XtX") + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.line = element_line(size = 1), plot.title = element_text(size = 14, face = "bold"), 
                     axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 30, color = "black"), 
                     axis.title = element_text(size = 18, face = "bold"), legend.justification = c(1, 0), legend.position = c(0.85, 0.1), 
                     legend.text = element_text(size = 20), legend.title = element_text(size = 20)) + scale_y_continuous(limits = c(0, 10)) + scale_x_continuous(limits = c(0, 30)) + scale_color_manual(values = c("#CCCCCC", "#000000"), name = "Outlier", breaks = c("No", "Yes"), labels = c("No", "Yes")) + scale_size(guide = "none") + guides(color = guide_legend(override.aes = list(size = 4)))
SSTmean_XtX_annotated_plot

#SSTmin v. XtX plot
SSTmin_XtX_plot <- ggplot(data = BFs, aes(x = SST_Min, y = M_XtX)) + geom_point(aes(size = 0.5)) + 
  geom_hline(yintercept = 5.78, size = 2, color = "#666666", linetype = "dashed") + 
  geom_vline(xintercept = 15, size = 2, color = "#666666", linetype = "dashed")
SSTmin_XtX_annotated_plot <- SSTmin_XtX_plot + ggtitle("SST Min BF v. XtX") + labs(x = "SST Min BF", y = "XtX") + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.line = element_line(size = 1), plot.title = element_text(size = 14, face = "bold"), 
                     axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 30, color = "black"), 
                     axis.title = element_text(size = 18, face = "bold"), legend.justification = c(1, 0), legend.position = c(0.85, 0.1), 
                     legend.text = element_text(size = 20), legend.title = element_text(size = 20)) + scale_y_continuous(limits = c(0, 10)) + scale_x_continuous(limits = c(0, 30)) + scale_color_manual(values = c("#CCCCCC", "#000000"), name  = "Outlier", breaks = c("No", "Yes"), labels = c("No", "Yes")) + scale_size(guide = "none") + guides(color = guide_legend(override.aes = list(size = 4)))
SSTmin_XtX_annotated_plot

#SSTmax v. XtX plot
SSTmax_XtX_plot <- ggplot(data = BFs, aes(x = SST_Max, y = M_XtX)) + geom_point(aes(size = 0.5)) + 
  geom_hline(yintercept = 5.78, size = 2, color = "#666666", linetype = "dashed") + 
  geom_vline(xintercept = 15, size = 2, color = "#666666", linetype = "dashed")
SSTmax_XtX_annotated_plot <- SSTmax_XtX_plot + ggtitle("SST Max BF v. XtX") + labs(x = "SST Max BF", y = "XtX") + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.line = element_line(size = 1), plot.title = element_text(size = 14, face = "bold"), 
                     axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 30, color = "black"), 
                     axis.title = element_text(size = 18, face = "bold"), legend.justification = c(1, 0), legend.position = c(0.85,0.1), 
                     legend.text = element_text(size = 20), legend.title = element_text(size = 20)) + scale_y_continuous(limits = c(0, 10)) + scale_x_continuous(limits = c(0, 30)) + scale_color_manual(values = c("#CCCCCC", "#000000"), name = "Outlier", breaks = c("No", "Yes"), labels = c("No", "Yes")) + scale_size(guide = "none") + guides(color = guide_legend(override.aes = list(size = 4)))
SSTmax_XtX_annotated_plot

#SSTmean v. fst plot
#SSTmean_fst_plot <- ggplot(data=BFs, aes(x=SST_Mean, y=fst)) + geom_point() + geom_hline(yintercept=0.15) + geom_vline(xintercept=20) + geom_vline(xintercept=15, linetype="dashed") + theme(text = element_text(size=20))
#SSTmean_fst_plot

#SSTmin v. fst plot
#SSTmin_fst_plot <- ggplot(data=BFs, aes(x=SST_Min, y=fst)) + geom_point() + geom_hline(yintercept=0.15) + geom_vline(xintercept=20) + geom_vline(xintercept=15, linetype="dashed") + theme(text = element_text(size=20))
#SSTmin_fst_plot

#SSTmax v. fst plot
#SSTmax_fst_plot <- ggplot(data=BFs, aes(x=SST_Max, y=fst)) + geom_point() + geom_hline(yintercept=0.15) + geom_vline(xintercept=20) + geom_vline(xintercept=15, linetype="dashed") + theme(text = element_text(size=20))
#SSTmax_fst_plot
