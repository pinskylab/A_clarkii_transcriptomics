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
library(tidyverse)
library(boot)

#read in data
SNPs_mac2 <- read.genetix("../../VCFs_and_PLINK/output.hicov2.snps.only.mac2_genetix.gtx") #convert genetix file to genind object for use in R (filtered SNPs for those with mac >/= 2) (maf >/= 0.05))
allele_counts <- read_table2("Data/clownfish_mac2_geno.txt", col_names = TRUE) #read in allele count input data for BayPass
locnames <- read_table2("Data/Loc_Names_mac2.txt", col_names = TRUE) #read in contig & bp for loci passing mac2 filter
  dim(locnames) #4212 x 1

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
  diversity_tot <- diversity_tot[2:11] #remove snp_id indexing
write.csv(diversity_tot, file = "Data/mac2_diversity.csv")

###############################################################################################################################################

######## Bootstrapping for Fis ########

#write function to calculate mean
samp_mean <- function(x, i) {
  mean(x[i])
}

#bootstrap for J pop Fis
J_Fis <- data.frame(diversity_tot$J_Fis)
  colnames(J_Fis) <- c("Fis")
J_Fis <- subset(J_Fis, !is.na(J_Fis$Fis)) #remove NAs
J_Fismean <- mean(J_Fis$Fis) #calculate mean (-0.06498  )

boot_J_Fis <- boot(data = J_Fis$Fis, statistic = samp_mean, R = 1000) #1000 permutations of pi
J_95ci_Fis <- boot.ci(boot_J_Fis, conf = 0.95, type = "norm") #get 95% CI for pi
J_95ci_Fis_normal <- J_95ci_Fis$normal #pull out normal distribution  2.5 & 97.5 percentiles for pi

#bootstrap for P pop Fis
P_Fis <- data.frame(diversity_tot$P_Fis)
  colnames(P_Fis) <- c("Fis")
P_Fis <- subset(P_Fis, !is.na(P_Fis$Fis)) #remove NAs
P_Fismean <- mean(P_Fis$Fis) #calculate mean (0.013314)

boot_P_Fis <- boot(data = P_Fis$Fis, statistic = samp_mean, R = 1000) #1000 permutations of pi
P_95ci_Fis <- boot.ci(boot_P_Fis, conf = 0.95, type = "norm") #get 95% CI for pi
P_95ci_Fis_normal <- P_95ci_Fis$normal #pull out normal distribution  2.5 & 97.5 percentiles for pi

#bootstrap for N pop Fis
N_Fis <- data.frame(diversity_tot$I_Fis)
  colnames(N_Fis) <- c("Fis")
N_Fis <- subset(N_Fis, !is.na(N_Fis$Fis)) #remove NAs
N_Fismean <- mean(N_Fis$Fis) #calculate mean (-0.002233)

boot_N_Fis <- boot(data = N_Fis$Fis, statistic = samp_mean, R = 1000) #1000 permutations of pi
N_95ci_Fis <- boot.ci(boot_N_Fis, conf = 0.95, type = "norm") #get 95% CI for pi
N_95ci_Fis_normal <- N_95ci_Fis$normal #pull out normal distribution  2.5 & 97.5 percentiles for pi

#create summary tables
Fis_mean <- as.data.frame(c(J_Fismean, P_Fismean, N_Fismean))
rownames(Fis_mean) <- c("Japan", "Philippines", "Indonesia")

Fis_norm_ci <- rbind(J_95ci_Fis_normal, P_95ci_Fis_normal, N_95ci_Fis_normal) #combine df w/ci info for each pop into one dataframe
Fis_sum <- cbind(Fis_mean, Fis_norm_ci)
  colnames(Fis_sum) <- c("mean", "ci", "2.5_per", "97.5_per")
Fis_sum$Pop <- c("Japan", "Philippines", "Indonesia")

Fis_sum$diff_lower <- Fis_sum$mean - Fis_sum$`2.5_per` #calculate diff btwn sample mean and 2.5 percentile for CI visualization
Fis_sum$diff_upper <- Fis_sum$`97.5_per` - Fis_sum$mean # calculate diff btwn sample mean and 97.5 percentile for CI visualization

#write out data
write.csv(Fis_sum, "Data/Fis_cis_mac2.csv")