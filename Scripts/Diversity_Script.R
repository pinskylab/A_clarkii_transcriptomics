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
library(diveRsity)
library(tidyverse)

#read in data
SNPs_mac2 <- read.genetix("../../VCFs_and_PLINK/output.hicov2.snps.only.mac2_genetix.gtx") #convert genetix file to genind object for use in R (filtered SNPs for those with mac >/= 2) (maf >/= 0.05))
allele_counts <- read_table2("Data/clownfish_mac2_geno.txt", col_names = TRUE) #read in allele count input data for BayPass
locnames <- read_table2("Data/Loc_Names_mac2.txt", col_names = TRUE) #read in contig & bp for loci passing mac2 filter
dim(locnames) #800 x 2 (Name & Position)

################################################################################################################################################

######## Compute diversity statistics w/only SNPs w/mac >/= 2 ########

#format data
snps <- c(1:4212) #create column for indexing
locnames$snp_id <- snps
allele_counts$snp_id <- snps

allele_counts <- merge(allele_counts, locnames, by = "snp_id") #merge allele_counts and locnames so know which loci corresponds with what allele frequencies
allele_counts <- allele_counts[,2:8] #remove snp_id index

#calculate Ho per locus in each pop
sum_stats_mac2 <- basic.stats(SNPs_mac2) #hierfstat function --> gives Ho, Hs, Fis on locus basis/pop as well as some pop-wide summary stats (Hs from Nei & Chesser 1983)
J_Ho <- data.frame(sum_stats_mac2$Ho[,1]) #pull out Ho for J pop
I_Ho <- data.frame(sum_stats_mac2$Ho[,2])
P_Ho <- data.frame(sum_stats_mac2$Ho[,3])

#calculate Hs (~He) per locus in each pop
J_Hs <- apply(allele_counts[,1:2], 1, heterozygosity) #pegas heterozygosity function to calculate Nei's average gene diversity w/correction for small sample sizes
P_Hs <- apply(allele_counts[,3:4], 1, heterozygosity)
I_Hs <- apply(allele_counts[,5:6], 1, heterozygosity)

#calculate Fis per locus in each pop
J_header <- c("J_Ho", "J_Hs")
J_diversity <- data.frame(J_Ho, J_Hs) #combine het metrics into same data frame
colnames(J_diversity) <- J_header #set header
J_diversity$J_Fis <- 1 - (J_diversity$J_Ho/J_diversity$J_Hs) #calculate Fis per locus

P_header <- c("P_Ho", "P_Hs")
P_diversity <- data.frame(P_Ho, P_Hs)
colnames(P_diversity) <- P_header
P_diversity$P_Fis <- 1 - (P_diversity$P_Ho/P_diversity$P_Hs)

I_header <- c("I_Ho", "I_Hs")
I_diversity <- data.frame(I_Ho, I_Hs)
colnames(I_diversity) <- I_header
I_diversity$I_Fis <- 1 - (I_diversity$I_Ho/I_diversity$I_Hs)

#combine to one data frame and write out
diversity_tot <- data.frame(J_diversity, P_diversity, I_diversity)
diversity_tot$snp_id <- snps
diversity_tot <- merge(diversity_tot, locnames, by="snp_id")
diversity_tot <- diversity_tot[2:11]
write.csv(diversity_tot, file = "Data/mac2_diversity.csv")

###############################################################################################################################################

######## (Optional) inbreeding measure (Fis) with adegenet ########

#set population names
popNames(SNPs_mac2) <- c("Japan", "Indonesia", "Philippines")
popNames(SNPs_mac2) #check to make sure popnames are right

#subset dataframe to each population
J_inds <- SNPs_mac2[1:8, ]
J_sum <- summary(J_inds) #get summary to make sure subsetting correctly
I_inds <- SNPs_mac2[9:15, ]
I_inds <- SNPs_mac2[16:25, ]

#inbreeding measure
J_inbred <- inbreeding(J_inds) #creates likelihood distribution for Fis of each locus
J_inbred_mean_mac2 <- sapply(J_inbred_mac2, mean) #mean likelihood
hist(J_inbred_mean_mac2) #plot

P_inbred <- inbreeding(P_inds)
P_inbred_mean <- sapply(P_inbred, mean)
hist(P_inbred_mean)

I_inbred <- inbreeding(I_inds)
I_inbred_mean <- sapply(I_inbred, mean)
hist(I_inbred_mean)

Genepop_mac2 <- readGenepop(infile = "Data/output.hicov2.snps.only.mac2_genpop.gen", gp = 3, bootstrap = FALSE)
basicstats <- divBasic(infile = "Data/output.hicov2.snps.only.mac2_genpop.gen", outfile = NULL, gp = 3, bootstraps = 0)
allele <- Genepop_mac2$all_alleles

