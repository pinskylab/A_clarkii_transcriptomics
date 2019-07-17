################################################### Script for Genetic Diversity Estimates  #######################################################

#Created for transcriptome project
#Uses adegenet package and dependencies

#################################################################################################################################################

######## Set-up ########

#set working directory
setwd("C:/Users/Rene/Dropbox/Pinsky_Lab/Transcriptome_Proj/R_scripts/A_clarkii_transcriptomics/")
getwd()

remove(list = ls())

#load libraries
library(adegenet)
library(pegas)
library(hierfstat)

#read in data
SNPs <- read.genepop("output.hicov2.snps.only_genepop.gen", ncode = 2) #converts genepop file to genind object for use in R
SNPs #get information on genind object to make sure is correct

################################################################################################################################################

######## Compute diversity statistics ########

#set population names
popNames(SNPs) <- c("Japan", "Philippines", "Indonesia")
popNames(SNPs) #check to make sure popnames are right

#subset dataframe to each population
J_inds <- SNPs[1:8, ]
P_inds <- SNPs[9:18, ]
I_inds <- SNPs[19:25, ]

#summary information for each population
sum_stats <- basic.stats(SNPs) #hierfstat function --> gives Ho, He, Fis on locus basis/pop as well as some pop-wide summary stats

#subset data
J_info <- data.frame(sum_stats$Ho[,1], sum_stats$Hs[,1], sum_stats$Fis[,1]) #pull out within-pop stats for Japanese pop
P_info <- data.frame(sum_stats$Ho[,2], sum_stats$Hs[,2], sum_stats$Fis[,2])
I_info <- data.frame(sum_stats$Ho[,3], sum_stats$Hs[,3], sum_stats$Fis[,3])
header <- c("Ho", "Hs", "Fis")
colnames(J_info) <- header #rename columns
colnames(P_info) <- header
colnames(I_info) <- header

#write out
write.table(J_info, file = "Japan_diversity.csv", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(P_info, file = "Philippines_diversity.csv", sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(I_info, file = "Indonesia_diversity.csv", sep = "\t", row.names = TRUE, col.names = TRUE)

#inbreeding measure
J_inbred <- inbreeding(J_inds) #creates likelihood distribution for Fis of each locus
J_inbred_mean <- sapply(J_inbred, mean) #mean likelihood
hist(J_inbred_mean) #plot

P_inbred <- inbreeding(P_inds)
P_inbred_mean <- sapply(P_inbred, mean)
hist(P_inbred_mean)

I_inbred <- inbreeding(I_inds)
I_inbred_mean <- sapply(I_inbred, mean)
hist(I_inbred_mean)