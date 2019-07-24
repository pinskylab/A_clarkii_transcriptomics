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
ssites_df <- data.frame(country, lat, long, sample_size)
ssites_df$sample_size <- factor(ssites_df$sample_size) #make sample size a factor incase want to scale point size by sample size
str(ssites_df) #check df structure

#create annotated map of SE Pacific
JNP_plot <- ggplot(JNP_map, aes(x = long, y = lat, group = group)) + 
  geom_polygon(fill = "lightgray", color = "black") +  
  coord_cartesian(xlim = c(115, 140), ylim = c(-5, 40)) + 
  geom_point(data = ssites_df, aes(x = long, y = lat, fill = "black", size = 4), inherit.aes = FALSE)
JNP_plot
JNP_plot_annotated <- JNP_plot + xlab("Longitude (?)") + ylab("Latitude (?)") + 
  annotate("text", x = 134, y = 31, label = "Japan", size = 6) + 
  annotate("text", x = 129, y = 14, label = "Philippines", size = 6) + 
  annotate("text", x = 122, y = 3, label = "Indonesia", size = 6) + theme_bw() + 
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(size = 1), axis.ticks = element_line(color = "black", size = 1),
        axis.text = element_text(size = 16, color = "black"), axis.title = element_text(size = 18, face = "bold"))
JNP_plot_annotated