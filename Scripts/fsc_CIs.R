#################################### Script to calculate CIs for fsc runs ##########################################

#Script to calculate CIs from fastsimcoal2 runs (boostrapped)
#Modifying Jennifer Hoey's script (https://github.com/pinskylab/NePADE/blob/master/best_lhoods.R)

#################################################################################################################################################

######## Set-up ########

#set working directory
setwd("C:/Users/rclar/Dropbox/Pinsky_Lab/Transcriptome_Proj/R_scripts/A_clarkii_transcriptomics/")
getwd()

remove(list = ls())

#load libraries
library(tidyverse)

#read in data
fsc_bootstraps <- read.csv("Data/fsc_maxLhood_CI_summary.csv", header = TRUE) #read in data from fsc runs

#################################################################################################################################################

#calculate 2.5 & 97.5 quantiles
quantiles <- data.frame(t(apply(fsc_bootstraps, MARGIN = 2, FUN = quantile, c(0.025, 0.50, 0.975)))) #calculates quantiles for each column in diversity_boot and transform df so in tidy format
colnames(quantiles) <- c("2.5_per", "median", "97.5_per")
