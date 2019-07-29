##################################################### Script for Plotting Sampling Sites  #########################################################

#Created for transcriptome project with sites in Japan, Philippines, & Indonesia

#################################################################################################################################################

######## Set-up ########

#set working directory
setwd("C:/Users/Rene/Dropbox/Pinsky_Lab/Transcriptome_Proj/R_scripts/A_clarkii_transcriptomics/")
getwd()

remove(list = ls())

#load libraries
library(maps)
library(mapdata)
library(tidyverse)
#library(ggsn)

#read in data
JNP <- c("Japan", "Philippines", "Indonesia")
JNP_map <- map_data('world', region = JNP)

################################################################################################################################################

######## Create map ########

#create data frame for sampling site lat long values
country <- c("Japan", "Philippines", "Indonesia")
lat <- c(33.00223, 10.87304, 0.652217)
long <- c(132.5026, 124.7122, 119.739)
sample_size <- c(8, 10, 7)
location <- c("edge", "core", "core")
ssites_df <- data.frame(country, lat, long, sample_size, location)
ssites_df$sample_size <- factor(ssites_df$sample_size) #make sample size a factor incase want to scale point size by sample size
str(ssites_df) #check df structure

#create annotated map of SE Pacific
JNP_plot <- ggplot(JNP_map, aes(x = long, y = lat, group = group)) + 
  geom_polygon(fill = "gray", color = "white") +  
  coord_cartesian(xlim = c(115, 140), ylim = c(-5, 40)) + 
  geom_point(data = ssites_df, aes(x = long, y = lat, color = location), size = 8, inherit.aes = FALSE) + 
  scale_color_manual(values = c("darkorange2", "dodgerblue4"))
JNP_plot
JNP_plot_annotated <- JNP_plot + xlab("Longitude (°)") + ylab("Latitude (°)") + 
  annotate("text", x = 134, y = 31, label = "Japan", size = 8) + 
  annotate("text", x = 129, y = 14, label = "Philippines", size = 8) + 
  annotate("text", x = 122, y = 3, label = "Indonesia", size = 8) + theme_bw() + 
  annotate("text", x = 134, y = 29, label = "n = 8", size = 8) + 
  annotate("text", x = 130, y = 12, label = "n = 10", size = 8) + 
  annotate("text", x = 121, y = 4.75, label = "n = 7", size = 8) + 
  guides(color = guide_legend(override.aes = list(size = 8))) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        axis.line = element_line(size = 1), axis.ticks = element_line(color = "black", size = 1),
        axis.text = element_text(size = 20, color = "black"), axis.title = element_text(size = 24, face = "bold"), 
        legend.position = c(0.15, 0.95), legend.text = element_text(size = 20), legend.title = element_blank())
JNP_plot_annotated
