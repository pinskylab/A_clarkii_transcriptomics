################################################### Script for Tajima D  #######################################################

#Created for transcriptome project
#Uses Tajima's D estimates from VCFtools

#################################################################################################################################################

######## Set-up ########

#set working directory
setwd("C:/Users/rclar/Dropbox/Pinsky_Lab/Transcriptome_Proj/R_scripts/A_clarkii_transcriptomics/")
getwd()

remove(list = ls())

#load libraries
library(tidyverse)

#read in data
J_TD <- read.csv("Data/JTajimasD.csv", header = TRUE) #read in data from vcftools
P_TD <- read.csv("Data/PTajimasD.csv", header =TRUE)
N_TD <- read.csv("Data/NTajimasD.csv", header = TRUE)
All_TD <- read.csv("Data/allTajimasD.csv", header = TRUE)
outlierseq <- read.csv("Data/outlier_sequences.csv", header = TRUE)

#clean data
JTD_data_clean <- J_TD[which(J_TD$TajimaD!='NA'),]
PTD_data_clean <- P_TD[which(P_TD$TajimaD!='NA'),]
NTD_data_clean <- N_TD[which(N_TD$TajimaD!='NA'),]
AllTD_data_clean <- All_TD[which(All_TD$TajimaD!='NA'),] #bc some windows with no SNPs at all

#coercing Tajima D to numeric
JTD_data_clean$TajimaD <- as.numeric(JTD_data_clean$TajimaD)
PTD_data_clean$TajimaD <- as.numeric(PTD_data_clean$TajimaD)
NTD_data_clean$TajimaD <- as.numeric(NTD_data_clean$TajimaD)
AllTD_data_clean$TajimaD <- as.numeric(AllTD_data_clean$TajimaD)

#################################################################################################################################################

######## Adding monomorphic sequences back in ########

#make sure outlierseqs are correct
AllTD_data_clean$uniqseq <- paste(AllTD_data_clean$CHROM, "-", AllTD_data_clean$BIN_START) #create uniq column to compare dfs by
AllTD_data_clean$Pop <- c(rep("All", times = 2255))
AllTD_data_clean$NUM <- c(1:2255)
setdiff(outlierseq$CONTIG, AllTD_data_clean$uniqseq) #should = 0, as all should be present in All_TD

######## J pop calculations ########

#number of monomorphic sequences not included
nrow(AllTD_data_clean) - nrow(JTD_data_clean) #531 --> need to add this many instances of 0 into dataframe

#create uniq column to compare dfs by
JTD_data_clean$uniqseq <- paste(JTD_data_clean$CHROM, "-", JTD_data_clean$BIN_START)

#create dataframe of sequences to add back
contigs_needed_J <- as.data.frame(setdiff(AllTD_data_clean$uniqseq, JTD_data_clean$uniqseq))
colnames(contigs_needed_J) <- c("uniqseq")

contigs_needed_J <- separate(contigs_needed_J, uniqseq, c("CHROM", "BIN_START"), sep = "-", remove = FALSE) #separate uniqseq out to original columns
contigs_needed_J$N_SNPS <- c(rep(0, times = 531)) #add N_SNPS column
contigs_needed_J$TajimaD <- c(rep("NA", times = 531)) #add TajimaD column

#create full Tajima D df
JTD_nonum <- JTD_data_clean[, 2:6]
full_JTD <- rbind(JTD_nonum, contigs_needed_J)
full_JTD <- full_JTD[order(full_JTD$uniqseq), ] #reorder by contig
nrow(full_JTD) #check equals 2255 (# total sequences)
full_JTD$Pop <- c(rep("Japan", times = 2255))
full_JTD$NUM <- c(1:2255)

#pull outlier seq & stats for J pop
J_outlierseq <- full_JTD[full_JTD$uniqseq %in% outlierseq$CONTIG, ]
J_outlierseq <- J_outlierseq[, 1:6]
J_outlierseq$NUM <- c(1:73) #re-numbering for plotting
J_outlierseq$TajimaD <- as.numeric(J_outlierseq$TajimaD)

######## P pop calculations ########

#number of monomorphic sequences not included
nrow(AllTD_data_clean) - nrow(PTD_data_clean) #262 --> need to add this many instances of 0 into dataframe

#create uniq column to compare dfs by
PTD_data_clean$uniqseq <- paste(PTD_data_clean$CHROM, "-", PTD_data_clean$BIN_START)

#create dataframe of sequences to add back
contigs_needed_P <- as.data.frame(setdiff(AllTD_data_clean$uniqseq, PTD_data_clean$uniqseq))
colnames(contigs_needed_P) <- c("uniqseq")

contigs_needed_P <- separate(contigs_needed_P, uniqseq, c("CHROM", "BIN_START"), sep = "-", remove = FALSE)
contigs_needed_P$N_SNPS <- c(rep(0, times = 262))
contigs_needed_P$TajimaD <- c(rep("NA", times = 262))

#create full Tajima D df
PTD_nonum <- PTD_data_clean[, 2:6]
full_PTD <- rbind(PTD_nonum, contigs_needed_P)
full_PTD <- full_PTD[order(full_PTD$uniqseq), ]
nrow(full_PTD) #check equals 2255 (# total sequences)
full_PTD$Pop <- c(rep("Philippines", times = 2255))
full_PTD$NUM <- c(1:2255)

#pull outlier seq & stats for P pop
P_outlierseq <- full_PTD[full_PTD$uniqseq %in% outlierseq$CONTIG, ]
P_outlierseq <- P_outlierseq[, 1:6]
P_outlierseq$NUM <- c(1:73)
P_outlierseq$TajimaD <- as.numeric(P_outlierseq$TajimaD)

######## N pop calculations ########

#number of monomorphic sequences not included
nrow(AllTD_data_clean) - nrow(NTD_data_clean) #472 --> need to add this many instances of 0 into dataframe

#create uniq column to compare dfs by
NTD_data_clean$uniqseq <- paste(NTD_data_clean$CHROM, "-", NTD_data_clean$BIN_START)

#create dataframe of sequences to add back
contigs_needed_N <- as.data.frame(setdiff(AllTD_data_clean$uniqseq, NTD_data_clean$uniqseq))
colnames(contigs_needed_N) <- c("uniqseq")

contigs_needed_N <- separate(contigs_needed_N, uniqseq, c("CHROM", "BIN_START"), sep = "-", remove = FALSE)
contigs_needed_N$N_SNPS <- c(rep(0, times = 472))
contigs_needed_N$TajimaD <- c(rep("NA", times = 472))

#create full Tajima D df
NTD_nonum <- NTD_data_clean[, 2:6]
full_NTD <- rbind(NTD_nonum, contigs_needed_N)
full_NTD <- full_NTD[order(full_NTD$uniqseq), ]
nrow(full_NTD) #check equals 2255 (# total sequences)
full_NTD$Pop <- c(rep("Indonesia", times = 2255))
full_NTD$NUM <- c(1:2255)

#pull outlier seq & stats for N pop
N_outlierseq <- full_NTD[full_NTD$uniqseq %in% outlierseq$CONTIG, ]
N_outlierseq <- N_outlierseq[, 1:6]
N_outlierseq$NUM <- c(1:73)
N_outlierseq$TajimaD <- as.numeric(N_outlierseq$TajimaD)

######## Combine pop data ########

#pull outlier seq & stats for pops combined
all_outlierseq <- AllTD_data_clean[AllTD_data_clean$uniqseq %in% outlierseq$CONTIG, ]
all_outlierseq <- all_outlierseq[order(all_outlierseq$uniqseq), ] #reorder by contig
all_outlierseq <- all_outlierseq[, 2:7]
all_outlierseq$NUM <- c(1:73)
all_outlierseq$TajimaD <- as.numeric(all_outlierseq$TajimaD)

#combine TD for all seqs
AllTD_data_clean <- AllTD_data_clean[order(AllTD_data_clean$uniqseq), ] #reorder by contig
TD_only_all <- rbind(full_JTD, full_PTD, full_NTD, AllTD_data_clean)
lapply(TD_only_all, class) #check character class for columns
TD_only_all$TajimaD <- as.numeric(TD_only_all$TajimaD)

#combine TD for outlier seqs only
TD_outlier_only_all <- rbind(J_outlierseq, P_outlierseq, N_outlierseq, all_outlierseq)
lapply(TD_outlier_only_all, class) #check character class for columns

#write out
write.csv(TD_only_all, "Data/TajimasD_combined_full.csv")
write.csv(TD_outlier_only_all, "Data/TajimasD_outlier_combined_full.csv")

#################################################################################################################################################

######## Calculate means by pop ########

######## J pop calculations ########

#calculate mean Tajimas D for all seq
J_mean_TD_all <- mean(JTD_data_clean$TajimaD) #-0.0167

#calculate mean Tajimas D for outlier seq
J_outlierseq_clean <- J_outlierseq[which(J_outlierseq$TajimaD!='NA'),]
J_mean_TD_outlierseq <- mean(J_outlierseq_clean$TajimaD) #0.4231

######## P pop calculations ########

#calculate mean Tajimas D for all seq
P_mean_TD_all <- mean(PTD_data_clean$TajimaD) #-0.0003

#calculate mean Tajimas D for outlier seq
P_outlierseq_clean <- P_outlierseq[which(P_outlierseq$TajimaD!='NA'),]
P_mean_TD_outlierseq <- mean(P_outlierseq_clean$TajimaD) #0.7385

######## N pop calculations ########

#calculate mean Tajimas D for all seq
N_mean_TD_all <- mean(NTD_data_clean$TajimaD) #0.0559

#calculate mean Tajimas D for outlier seq
N_outlierseq_clean <- N_outlierseq[which(N_outlierseq$TajimaD!='NA'),]
N_mean_TD_outlierseq <- mean(N_outlierseq_clean$TajimaD) #0.6321

######## All pop calculations ########

#calculate mean Tajimas D for all seq
Tot_mean_TD_all <- mean(AllTD_data_clean$TajimaD) #0.1109

#calculate mean Tajimas D for outlier seq
all_outlierseq_clean <- all_outlierseq[which(all_outlierseq$TajimaD!='NA'),]
Tot_mean_TD_outlierseq <- mean(all_outlierseq_clean$TajimaD) #1.102

#################################################################################################################################################

######## Visualize data ########

#all seqs by pop
all_plot <- ggplot(data = TD_only_all, aes(x = NUM, y = TajimaD, color = Pop)) + 
  geom_point() + geom_hline(yintercept = 0, color = "black", size = 1)
all_plot_annotated <- all_plot + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(size = 1), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 14, color = "black"), 
        axis.title = element_text(size = 14, face = "bold"), legend.position = "top", 
        legend.text = element_text(size = 12), legend.title = element_text(size = 12))
all_plot_annotated

#outliers only by pop
outliers_plot <- ggplot(data = TD_outlier_only_all, aes(x = NUM, y = TajimaD, color = Pop)) + 
  geom_point() + geom_hline(yintercept = 0, color = "black", size = 1)
outliers_plot_annotated <- outliers_plot + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(size = 1), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 14, color = "black"), 
        axis.title = element_text(size = 14, face = "bold"), legend.position = "top", 
        legend.text = element_text(size = 12), legend.title = element_text(size = 12))
outliers_plot_annotated

