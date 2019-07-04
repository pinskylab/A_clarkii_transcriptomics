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
library(readr)
BFs <- read_table2("aux1_all_summary_betai.txt", col_names = TRUE) #use read_table2 because uneven amount of white spaces between columns

################################################################################################################################################

######## Pull info out from BF dataframe ########

#pull BFs for each covariable
BFs_SSSmean <- c(BFs[1:5729, 'BF(dB)'])
BFs_SSTmean <- c(BFs[5730:11458, 'BF(dB)'])
BFs_SSTmin <- c(BFs[11459:17187, 'BF(dB)'])
BFs_SSTmax <- c(BFs[17188:22916, 'BF(dB)'])
BFs_lat <- c(BFs[22917:28645, 'BF(dB)'])

#pull marker numbers
markers <- c(BFs[1:5729, 'MRK'])

#Create data frame with markers & BFs for each covariable as rows
BFs_ALL <- data.frame(markers, BFs_SSSmean, BFs_SSTmean, BFs_SSTmin, BFs_SSTmax, BFs_lat)
names(BFs_ALL) <- c("snp_id", "SSS_Mean", "SST_Mean", "SST_Min", "SST_Max", "Lat")

######## Pull out putative candidate loci by env variable ########

#pull loci # with BF > 15 for each environmental covariate
#pull for col 2: SSS Mean
BFs_greater_than_15_SSSmean <- c(which(BFs_ALL[, 2] > 15))

#Pull for col 3: SST Mean
BFs_greater_than15_SSTmean <- c(which(BFs_ALL[, 3] > 15))

#Pull for col 4: SST Min
BFs_greater_than15_SSTmin <- c(which(BFs_ALL[, 4] > 15))

#Pull for col 5: SST Max
BFs_greater_than15_SSTmax <- c(which(BFs_ALL[, 5] > 15))

#Pull for col 6: Latitude
BFs_greater_than15_Latitude <- c(which(BFs_ALL[, 6] > 15))

######## Pull candidate loci names ########

#create vector of potential candidate loci
potential_loci_list <- BFs_greater_than_15_SSSmean

#add candidate loci flagged from other covariates
potential_loci_list <- append(potential_loci_list, BFs_greater_than15_SSTmean)
potential_loci_list <- append(potential_loci_list, BFs_greater_than15_SSTmin)
potential_loci_list <- append(potential_loci_list, BFs_greater_than15_SSTmax)
potential_loci_list <- append(potential_loci_list, BFs_greater_than15_Latitude)

######## Grab ID names of the contigs containing candidate SNPs ########

#read in loc names
locnames <- read_table2("Loc_Names_All.txt", col_names = TRUE)
dim(locnames) #800 x 2 (Name & Position)

snps <- c(1:5729) #create column to match "marker" column pulled earlier from aux1_summary_betai.txt
locnames$snp_id <- snps

#filter potential_loci_list so only one instance of each potential can loci
can_loci <- unique(sort(potential_loci_list)) #keep only unique loci (ignore repeats)
can_names <- locnames[c(can_loci),] #create dataframe of can loci that have been flagged with BF > 15
can_BFs <- BFs_ALL[c(can_loci),] #create dataframe of BFs for can loci for each env variable

#merge dataframes
can_all <- merge(can_names, can_BFs, by = "snp_id")
all <- merge(locnames, BFs_ALL, by = "snp_id")

#write out
write.table(can_loci, file = "aux1_all_candidate_loci.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(all, file = "aux1_all_BFs.txt", sep = "\t", row.names = FALSE, col.names = TRUE)


