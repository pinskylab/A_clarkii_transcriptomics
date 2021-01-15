############################################# Script for pulling candidate loci from BayPass runs #################################################

#assumes candidate loci have a BF > 15 for at least one env variable tested

#################################################################################################################################################

######## Set-up ############
remove(list = ls())

#set working directory
setwd("C:/Users/rclar/Dropbox/Pinsky_Lab/Transcriptome_Proj/R_scripts/A_clarkii_transcriptomics/")
getwd()

#load libraries
library(readr)

#read in BFs
BFs_ssmean <- read_table2("Data/mac2aux1_ssmean_summary_betai.txt", col_names = TRUE) #use read_table2 because uneven amount of white spaces between columns
BFs_sstmean <- read_table2("Data/mac2aux1_sstmean_summary_betai.txt", col_names = TRUE) #use read_table2 because uneven amount of white spaces between columns
BFs_sstmin <- read_table2("Data/mac2aux1_sstmin_summary_betai.txt", col_names = TRUE) #use read_table2 because uneven amount of white spaces between columns
BFs_sstmax <- read_table2("Data/mac2aux1_sstmax_summary_betai.txt", col_names = TRUE) #use read_table2 because uneven amount of white spaces between columns
BFs_lat <- read_table2("Data/mac2aux1_lat_summary_betai.txt", col_names = TRUE) #use read_table2 because uneven amount of white spaces between columns

#read in loc names
locnames <- read_table2("Data/Loc_Names_mac2.txt", col_names = FALSE)
  dim(locnames) #4213

################################################################################################################################################

######## Pull info out from BF dataframe ########

#pull BFs for each covariable
BFs_SSSmean <- c(BFs_ssmean[1:5718, 'BF(dB)'])
BFs_SSTmean <- c(BFs_sstmean[1:5718, 'BF(dB)'])
BFs_SSTmin <- c(BFs_sstmin[1:5718, 'BF(dB)'])
BFs_SSTmax <- c(BFs_sstmax[1:5718, 'BF(dB)'])
BFs_lat <- c(BFs_lat[1:5718, 'BF(dB)'])

#pull marker numbers
#markers <- c(BFs[1:5729, 'MRK'])
markers <- c(BFs_ssmean[1:5718, 'MRK'])

#Create data frame with markers & BFs for each covariable as rows
BFs_ALL <- data.frame(markers, BFs_SSSmean, BFs_SSTmean, BFs_SSTmin, BFs_SSTmax, BFs_lat)
  names(BFs_ALL) <- c("snp_id", "SSS_Mean", "SST_Mean", "SST_Min", "SST_Max", "Lat")

######## Pull out putative candidate loci by env variable ########

#pull loci # with BF > 20 for each environmental covariate
#pull for col 2: SSS Mean
BFs_greater_than20_SSSmean <- c(which(BFs_ALL[, 2] > 20))

#Pull for col 3: SST Mean
BFs_greater_than20_SSTmean <- c(which(BFs_ALL[, 3] > 20))

#Pull for col 4: SST Min
BFs_greater_than20_SSTmin <- c(which(BFs_ALL[, 4] > 20))

#Pull for col 5: SST Max
BFs_greater_than20_SSTmax <- c(which(BFs_ALL[, 5] > 20))

#Pull for col 6: Latitude
BFs_greater_than20_Latitude <- c(which(BFs_ALL[, 6] > 20))

######## Pull candidate loci names ########

#create vector of potential candidate loci
potential_loci_list <- BFs_greater_than20_SSSmean

#add candidate loci flagged from other covariates
potential_loci_list <- append(potential_loci_list, BFs_greater_than20_SSTmean)
potential_loci_list <- append(potential_loci_list, BFs_greater_than20_SSTmin)
potential_loci_list <- append(potential_loci_list, BFs_greater_than20_SSTmax)
potential_loci_list <- append(potential_loci_list, BFs_greater_than20_Latitude)

######## Grab ID names of the contigs containing candidate SNPs ########

#snps <- c(1:4213) #create column to match "marker" column pulled earlier from aux1_summary_betai.txt
snps <- c(1:4213)
locnames$snp_id <- snps

#filter potential_loci_list so only one instance of each potential can loci
can_loci <- unique(sort(potential_loci_list)) #keep only unique loci (ignore repeats)
can_names <- locnames[c(can_loci),] #create dataframe of can loci that have been flagged with BF > 20
can_BFs <- BFs_ALL[c(can_loci),] #create dataframe of BFs for can loci for each env variable

#merge dataframes
can_all <- merge(can_names, can_BFs, by = "snp_id")
all <- merge(locnames, BFs_ALL, by = "snp_id")

#write out
write.csv(can_all, file = "Data/mac2aux1_candidate_loci.csv")
write.csv(all, file = "Data/mac2aux1_all_loci.csv")