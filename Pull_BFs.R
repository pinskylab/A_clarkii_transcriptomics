#### Pulling loci out that have BF > 10 OR BF > 15####
remove(list = ls())

#Read in BF split by environmental covariate
library(readr)
BFs_sim <- read_table2("aux1_summary_betai.txt", col_names = TRUE) #use read_table2 because uneven amount of white spaces between columns

# Pull BFs for each covariable
BFs_SSSmean <- c(BFs_sim[1:5729, 'BF(dB)'])
BFs_SSTmean <- c(BFs_sim[5730:11458, 'BF(dB)'])
BFs_Lat <- c(BFs_sim[11459:17187, 'BF(dB)'])
#BFs_SSTmax <- c(BFs_sim[17188:22916, 'BF(dB)'])
#BFs_lat <- c(BFs_sim[22917:28645, 'BF(dB)'])

#Pull marker numbers
markers <- c(BFs_sim[1:5729, 'MRK'])

#Create data frame with markers & BFs for each covariable as rows
BFs_ALL <- data.frame(markers, BFs_SSSmean, BFs_SSTmean, BFs_Lat)
names(BFs_ALL) <- c("snp_id", "SSS Mean", "SST Mean", "Lat")

#Write out
#write.table(BFs_ALL, file = "aux1_all_BFs.txt", sep="\t", row.names = FALSE, col.names = TRUE)

#Pull loci # with BF > 10 OR BF > 15 for each environmental covariate
#Pull for col 2: SSS Mean
BFs_greater_than_10_SSSmean <- c(which(BFs_ALL[,2] > 10))
BFs_greater_than_15_SSSmean <- c(which(BFs_ALL[,2] > 15))

#Pull for col 3: SST Mean
BFs_greater_than10_SSTmean <- c(which(BFs_ALL[,3] > 10))
BFs_greater_than15_SSTmean <- c(which(BFs_ALL[,3] > 15))

#Pull for col 4: SST Min
BFs_greater_than10_SSTmin <- c(which(BFs_ALL[,4] > 10))
BFs_greater_than15_SSTmin <- c(which(BFs_ALL[,4] > 15))

#Pull for col 5: SST Max
#BFs_greater_than10_SSTmax <- c(which(BFs_ALL[,5] > 10))
#BFs_greater_than15_SSTmax <- c(which(BFs_ALL[,5] > 15))

#Pull for col 6: Latitude
#BFs_greater_than10_Latitude <- c(which(BFs_ALL[,6] > 10))
#BFs_greater_than15_Latitude <- c(which(BFs_ALL[,6] > 15))

#### Pulling Candidate Loci Names (loci w/BF > 10 for any covariate) ####
#Create vector of potential candidate loci
potential_loci_list <- BFs_greater_than_15_SSSmean

#Add candidate loci flagged from other covariates
potential_loci_list <- append(potential_loci_list, BFs_greater_than15_SSTmean)
potential_loci_list <- append(potential_loci_list, BFs_greater_than15_SSTmin)
#potential_loci_list <- append(potential_loci_list, BFs_greater_than15_SSTmax)
#potential_loci_list <- append(potential_loci_list, BFs_greater_than15_Latitude)

#### Grabbing the ID names of the contigs containing candidate SNPs ####
#Read in loc names
locnames <- read_table2("Loc_Names_All.txt", col_names = TRUE)
dim(locnames) #800 x 2 (Name & Position)
#rownames(locnames)
#locnames["snp"] <- rownames(locnames)
#dim(locnames) # 800 x 3 (Name, Position, snp)
snps <- c(1:5729)
locnames$snp_id <- snps

can_loci <- unique(sort(potential_loci_list)) #keep only unique loci (ignore repeats)
can_names <- locnames[c(can_loci),] #create matrix of contigs that have been flagged with BF > 10
can_BFs <- BFs_ALL[c(can_loci),]

can_all <- merge(can_names, can_BFs, by = "snp_id")

#Write out
write.table(can_loci, file = "aux1_all_candidate_loci.txt", sep="\t", row.names = FALSE, col.names = TRUE)

#merge BFs_ALL & locnames data frames
all <- merge(locnames, BFs_ALL, by = "snp_id")

#write out
write.table(all, file = "aux1_all_BFs.txt", sep="\t", row.names = FALSE, col.names = TRUE)
