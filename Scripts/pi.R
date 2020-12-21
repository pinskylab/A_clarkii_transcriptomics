################################################### Script for Pi  #######################################################

#Created for transcriptome project
#Uses window pi estimates from VCFtools and bootstraps to get 95% CIs 
#Section at bottom uses site pi estimates from VCFtools

#################################################################################################################################################

######## Set-up ########

#set working directory
setwd("C:/Users/rclar/Dropbox/Pinsky_Lab/Transcriptome_Proj/R_scripts/A_clarkii_transcriptomics/")
getwd()

remove(list = ls())

#load libraries
library(tidyverse)
library(boot)

#read in data
#window pi all calculated with mac2
all_pi <- read.csv("Data/allpi.csv", header = TRUE) #read in data from vcftools
J_pi <- read.csv("Data/Jwindowpi.csv", header = TRUE)
P_pi <- read.csv("Data/Pwindowpi.csv", header = TRUE)
N_pi <- read.csv("Data/Nwindowpi.csv", header = TRUE)

#################################################################################################################################################

######## Adding monomorphic sequences back in ########

######## J pop calculations ########

#number of monomorphic sequences not included
nrow(all_pi) - nrow(J_pi) #526 --> need to add this many instances of 0 into dataframe

#create uniq column to compare dfs by
J_pi$uniqseq <- paste(J_pi$CHROM, "-", J_pi$BIN_START, "-", J_pi$BIN_END)
all_pi$uniqseq <- paste(all_pi$CHROM, "-", all_pi$BIN_START, "-", all_pi$BIN_END)

#create dataframe of sequences to add back
contigs_needed_J <- as.data.frame(setdiff(all_pi$uniqseq, J_pi$uniqseq))
colnames(contigs_needed_J) <- c("uniqseq")

contigs_needed_J <- separate(contigs_needed_J, uniqseq, c("CHROM", "BIN_START", "BIN_END"), sep = "-", remove = FALSE) #separate uniqseq out to original columns
contigs_needed_J$N_VARIANTS <- c(rep(0, times = 526)) #add N_VARIANTS column
contigs_needed_J$PI <- c(rep(0, times = 526)) #add PI column

#create full pi df
full_Jpi <- rbind(J_pi, contigs_needed_J)
full_Jpi <- full_Jpi[order(full_Jpi$uniqseq), ] #reorder by contig
nrow(full_Jpi) #check equals 2253 (# total sequences)
full_Jpi$Pop <- c(rep("Japan", times = 2253))
full_Jpi$NUM <- c(1:2253)

######## P pop calculations ########

#number of monomorphic sequences not included
nrow(all_pi) - nrow(P_pi) #259 --> need to add this many instances of 0 into dataframe

#create uniq column to compare dfs by
P_pi$uniqseq <- paste(P_pi$CHROM, "-", P_pi$BIN_START, "-", P_pi$BIN_END)
all_pi$uniqseq <- paste(all_pi$CHROM, "-", all_pi$BIN_START, "-", all_pi$BIN_END)

#create dataframe of sequences to add back
contigs_needed_P <- as.data.frame(setdiff(all_pi$uniqseq, P_pi$uniqseq))
colnames(contigs_needed_P) <- c("uniqseq")

contigs_needed_P <- separate(contigs_needed_P, uniqseq, c("CHROM", "BIN_START", "BIN_END"), sep = "-", remove = FALSE)
contigs_needed_P$N_VARIANTS <- c(rep(0, times = 259))
contigs_needed_P$PI <- c(rep(0, times = 259))

#create full pi df
full_Ppi <- rbind(P_pi, contigs_needed_P)
full_Ppi <- full_Ppi[order(full_Ppi$uniqseq), ]
nrow(full_Ppi) #check equals 2253 (# total sequences)
full_Ppi$Pop <- c(rep("Philippines", times = 2253))
full_Ppi$NUM <- c(1:2253)

######## N pop calculations ########

#number of monomorphic sequences not included
nrow(all_pi) - nrow(N_pi) #469 --> need to add this many instances of 0 into dataframe

#create uniq column to compare dfs by
N_pi$uniqseq <- paste(N_pi$CHROM, "-", N_pi$BIN_START, "-", N_pi$BIN_END)
all_pi$uniqseq <- paste(all_pi$CHROM, "-", all_pi$BIN_START, "-", all_pi$BIN_END)

#create dataframe of sequences to add back
contigs_needed_N <- as.data.frame(setdiff(all_pi$uniqseq, N_pi$uniqseq))
colnames(contigs_needed_N) <- c("uniqseq")

contigs_needed_N <- separate(contigs_needed_N, uniqseq, c("CHROM", "BIN_START", "BIN_END"), sep = "-", remove = FALSE)
contigs_needed_N$N_VARIANTS <- c(rep(0, times = 469))
contigs_needed_N$PI <- c(rep(0, times = 469))

#create full pi df
full_Npi <- rbind(N_pi, contigs_needed_N)
full_Npi <- full_Npi[order(full_Npi$uniqseq), ]
nrow(full_Npi) #check equals 2253 (# total sequences)
full_Npi$Pop <- c(rep("Indonesia", times = 2253))
full_Npi$NUM <- c(1:2253)

######## Combine pop data ########

#add needed columns to all_pi df
all_pi <- all_pi[order(all_pi$uniqseq), ] #make sure ordered by uniqseq
all_pi$Pop <- c(rep("All", times = 2253))
all_pi$NUM <- c(1:2253)

#merge dataframes
pi_all <- rbind(full_Jpi, full_Ppi, full_Npi, all_pi)
lapply(pi_all, class) #check character class for columns

#write out combined data
write.csv(pi_all, "Data/pi_combined_full.csv")

#################################################################################################################################################

######## Bootstrap for 95% CIs ########

######## J pop calculations ########

J_pimean <- mean(full_Jpi$PI) #calculate mean

#mean function for bootstrapping
samp_mean <- function(x, i) {
  mean(x[i])
} #bc if use mean() in boot() throws trim error

#bootstrap for J pop pi
boot_J_pi <- boot(data = full_Jpi$PI, statistic = samp_mean, R = 1000) #1000 permutations of pi
J_95ci_pi <- boot.ci(boot_J_pi, conf = 0.95, type = "norm") #get 95% CI for pi
J_95ci_pi_normal <- J_95ci_pi$normal #pull out normal distribution  2.5 & 97.5 percentiles for pi

######## P pop calculations ########

P_pimean <- mean(full_Ppi$PI) #calculate mean

#bootstrap for P pop pi
boot_P_pi <- boot(data = full_Ppi$PI, statistic = samp_mean, R = 1000)
P_95ci_pi <- boot.ci(boot_P_pi, conf = 0.95, type = "norm")
P_95ci_pi_normal <- P_95ci_pi$normal

######## N pop calculations ########

N_pimean <- mean(full_Npi$PI) #calculate mean

#bootstrap for N pop pi
boot_N_pi <- boot(data = full_Npi$PI, statistic = samp_mean, R = 1000)
N_95ci_pi <- boot.ci(boot_N_pi, conf = 0.95, type = "norm")
N_95ci_pi_normal <- N_95ci_pi$normal

######## Create summary tables ########

pi_mean <- as.data.frame(c(J_pimean, P_pimean, N_pimean))
rownames(pi_mean) <- c("Japan", "Philippines", "Indonesia")

pi_norm_ci <- rbind(J_95ci_pi_normal, P_95ci_pi_normal, N_95ci_pi_normal) #combine df w/ci info for each pop into one dataframe
pi_sum <- cbind(pi_mean, pi_norm_ci)
colnames(pi_sum) <- c("mean", "ci", "2.5_per", "97.5_per")
pi_sum$Pop <- c("Japan", "Philippines", "Indonesia")

pi_sum$diff_lower <- pi_sum$mean - pi_sum$`2.5_per` #calculate diff btwn sample mean and 2.5 percentile for CI visualization
pi_sum$diff_upper <- pi_sum$`97.5_per` - pi_sum$mean # calculate diff btwn sample mean and 97.5 percentile for CI visualization

#write out data
write.csv(pi_sum, "Data/pi_cis.csv")

#################################################################################################################################################

######## Visualize data ########
#designed to be run separately from earlier sections

remove(list = ls())

#read in data
pi_all <- read.csv("Data/pi_combined_full.csv", header = TRUE, row.names = 1)
pi_sum <- read.csv("Data/pi_cis.csv", header = TRUE, row.names = 1)
outlierseq <- read.csv("Data/outlier_sequences_pi.csv", header = TRUE)

#add outlier status
setdiff(outlierseq$CONTIG, pi_all$uniqseq) #should = 0, as all should be present in pi_all
outliers <- pi_all[pi_all$uniqseq %in% outlierseq$CONTIG, ] #df with only outliers
nonoutliers <- pi_all[!(pi_all$uniqseq %in% outlierseq$CONTIG), ] #df with all other loci

#set seq status
outliers$Status <- c(rep("Outlier", times = 268))
nonoutliers$Status <- c(rep("Not_Outlier", times = 8744))
pi_all$Status <- c(rep("All", times = 9012))

#merge together
pi_all <- rbind(outliers, nonoutliers, pi_all)

######## Mean pi w/error bars plot ########

#plot of mean w/in pop pi w/95% CI error bars
mean_pi_plot <- ggplot(data = pi_sum, aes(x = Pop, y = mean)) + 
  geom_point(aes(size = 1), show.legend = FALSE) + 
  geom_errorbar(aes(ymin = mean - diff_lower, ymax = mean + diff_upper, width = 0.5, size = 0.5), show.legend = FALSE)
mean_pi_plot_annotated <- mean_pi_plot + ggtitle("Mean pi w/95% CI") + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 14, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 12, color = "black"))
mean_pi_plot_annotated

######## Boxplots ########

#ordering x-axis
pi_all$Pop2 <- factor(pi_all$Pop, levels = c("Japan", "Philippines", "Indonesia", "All")) #ordering X-axis

#boxplot
pi_boxplot <- ggplot(data = pi_all, aes(x = Pop2, y = PI, color = Status)) + 
  geom_boxplot()
pi_boxplot_annotated <- pi_boxplot + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(size = 1), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 14, color = "black"), 
        axis.title = element_text(size = 14, face = "bold"), legend.position = "top", 
        legend.text = element_text(size = 12), legend.title = element_text(size = 12))
pi_boxplot_annotated

#################################################################################################################################################

######## Site pi calculations ########
#calculated with mac = 2
#designed to be run separately from window pi calculations

remove(list = ls())

#load libraries
library(tidyverse)
library(boot)

#read in data
#site pi all calculated with mac1
allsites_pi <- read.csv("Data/allmac2_sitespi.csv", header = TRUE) #read in data from vcftools
Jsites_pi <- read.csv("Data/Jmac2_sitespi.csv", header = TRUE)
Nsites_pi <- read.csv("Data/Nmac2_sitespi.csv", header = TRUE)
Psites_pi <- read.csv("Data/Pmac2_sitespi.csv", header = TRUE)
contig_length <- read.csv("Data/Contig_length.csv", header = TRUE)
  colnames(contig_length) <- c("CHROM", "Contig.length")

######## Calculate mean pi by site ########

#read in total bp in transcriptome
total_bp_ref <- 45506763 #total bp in New_ref_N3.fa
total_bp_transcripts <- 1064586 #total bp in transcripts included in mac2 VCF

#### All pops ####
sum_allpi <- sum(allsites_pi$PI) #sum pi across all sites

#average pi across transcriptome
mean_allpi_total <- sum_allpi/total_bp_ref #0.0000235
mean_allpi_transcripts <- sum_allpi/total_bp_transcripts #0.000955

#add columns to write out together
allsites_pi$Pop <- c(rep("All", times = 4212))
allsites_pi$NUM <- c(1:4212)

#calculate pi by transcript
allsites_pi_aggregate <- aggregate(allsites_pi$PI~allsites_pi$CHROM, FUN = sum) #sum pi by transcript
  colnames(allsites_pi_aggregate) <- c("CHROM", "PI")

allsites_pi_aggregate <- merge(allsites_pi_aggregate, contig_length) #merge with contig length
allsites_pi_aggregate$pi_avg <- NA #create empty column to populate

for(i in 1:nrow(allsites_pi_aggregate)) {
  cat(paste(i, " ", sep = ''))
  {allsites_pi_aggregate$pi_avg[i] <- allsites_pi_aggregate$PI[i]/allsites_pi_aggregate$Contig.length[i]}
}

allsites_pi_aggregate$Pop <- c(rep("All", times = 999))
allsites_pi_aggregate$NUM <- c(1:999)

#### J pop ####
sum_Jpi <- sum(Jsites_pi$PI)
mean_Jpi_total <- sum_Jpi/total_bp_ref #0.0000198
mean_Jpi_transcripts <- sum_Jpi/total_bp_transcripts #0.000844

#add columns to write out together
Jsites_pi$Pop <- c(rep("Japan", times = 4212))
Jsites_pi$NUM <- c(1:4212)

#calculate pi by transcript
Jsites_pi_aggregate <- aggregate(Jsites_pi$PI~Jsites_pi$CHROM, FUN = sum) #sum pi by transcript
  colnames(Jsites_pi_aggregate) <- c("CHROM", "PI")

Jsites_pi_aggregate <- merge(Jsites_pi_aggregate, contig_length) #merge with contig length
Jsites_pi_aggregate$pi_avg <- NA #create empty column to populate

for(i in 1:nrow(Jsites_pi_aggregate)) {
  cat(paste(i, " ", sep = ''))
  {Jsites_pi_aggregate$pi_avg[i] <- Jsites_pi_aggregate$PI[i]/Jsites_pi_aggregate$Contig.length[i]}
}

Jsites_pi_aggregate$Pop <- c(rep("Japan", times = 999))
Jsites_pi_aggregate$NUM <- c(1:999)

#### P pop ####
sum_Ppi <- sum(Psites_pi$PI)
mean_Ppi_total <- sum_Ppi/total_bp_ref #0.0000225
mean_Ppi_transcripts <- sum_Ppi/total_bp_transcripts #0.000964

#add columns to write out together
Psites_pi$Pop <- c(rep("Philippines", times = 4212))
Psites_pi$NUM <- c(1:4212)

#calculate pi by transcript
Psites_pi_aggregate <- aggregate(Psites_pi$PI~Psites_pi$CHROM, FUN = sum) #sum pi by transcript
  colnames(Psites_pi_aggregate) <- c("CHROM", "PI")

Psites_pi_aggregate <- merge(Psites_pi_aggregate, contig_length) #merge with contig length
Psites_pi_aggregate$pi_avg <- NA #create empty column to populate

for(i in 1:nrow(Psites_pi_aggregate)) {
  cat(paste(i, " ", sep = ''))
  {Psites_pi_aggregate$pi_avg[i] <- Psites_pi_aggregate$PI[i]/Psites_pi_aggregate$Contig.length[i]}
}

Psites_pi_aggregate$Pop <- c(rep("Philippines", times = 999))
Psites_pi_aggregate$NUM <- c(1:999)

#### N pop ####
sum_Npi <- sum(Nsites_pi$PI)
mean_Npi_total <- sum_Npi/total_bp_ref #0.0000223
mean_Npi_transcripts <- sum_Npi/total_bp_transcripts #0.000953

#add columns to write out together
Nsites_pi$Pop <- c(rep("Indonesia", times = 4212))
Nsites_pi$NUM <- c(1:4212)

#calculate pi by transcript
Nsites_pi_aggregate <- aggregate(Nsites_pi$PI~Nsites_pi$CHROM, FUN = sum) #sum pi by transcript
  colnames(Nsites_pi_aggregate) <- c("CHROM", "PI")

Nsites_pi_aggregate <- merge(Nsites_pi_aggregate, contig_length) #merge with contig length
Nsites_pi_aggregate$pi_avg <- NA #create empty column to populate

for(i in 1:nrow(Nsites_pi_aggregate)) {
  cat(paste(i, " ", sep = ''))
  {Nsites_pi_aggregate$pi_avg[i] <- Nsites_pi_aggregate$PI[i]/Nsites_pi_aggregate$Contig.length[i]}
}

Nsites_pi_aggregate$Pop <- c(rep("Indonesia", times = 999))
Nsites_pi_aggregate$NUM <- c(1:999)

#merge dataframes
pi_site_all <- rbind(Jsites_pi, Psites_pi, Nsites_pi, allsites_pi)
lapply(pi_site_all, class) #check character class for columns

pi_site_all_aggregate <- rbind(Jsites_pi_aggregate, Psites_pi_aggregate, Nsites_pi_aggregate, allsites_pi_aggregate)
lapply(pi_site_all_aggregate, class)

#write out combined data
write.csv(pi_site_all, "Data/pi_combined_site_mac2_full.csv")
write.csv(pi_site_all_aggregate, "Data/pi_combined_site_full_mac2_aggregate.csv")

######## Bootstrap mean pi by site ########
#just working with transcript total bp for now

#sum function
samp_sum <- function(x, i) {
  sum(x[i])
} #bc if use sum() not summing properly

#### All pops ####
#bootstrap across loci to get sum of pi
boot_allsite_pi <- boot(data = allsites_pi$PI, statistic = samp_sum, R = 1000)
boot_allsite_pi_sum <- data.frame(boot_allsite_pi$t)
colnames(boot_allsite_pi_sum) <- c("boot_sum_pi")

#calculate avg pi by site
boot_allsite_pi_sum$mean_ref <- NA #create empty column to populate

#for loop to calculate avg pi across all VCF transcripts
for(i in 1:nrow(boot_allsite_pi_sum)) {
  cat(paste(i, " ", sep = ''))
  {boot_allsite_pi_sum$mean_ref[i] <- boot_allsite_pi_sum$boot_sum_pi[i]/total_bp_transcripts}
}

#get 95% CIs
allsite_95ci_pi <- quantile(boot_allsite_pi_sum$mean_ref, c(0.025, 0.975)) #(0.000938, 0.000974)
allsite_95ci_pi <- t(data.frame(allsite_95ci_pi))
rownames(allsite_95ci_pi) <- c("All")

#### J pop ####
#bootstrap across loci to get sum of pi
boot_Jsite_pi <- boot(data = Jsites_pi$PI, statistic = samp_sum, R = 1000)
boot_Jsite_pi_sum <- data.frame(boot_Jsite_pi$t)
colnames(boot_Jsite_pi_sum) <- c("boot_sum_pi")

#calculate avg pi by site
boot_Jsite_pi_sum$mean_ref <- NA #create empty column to populate

#for loop to calculate avg pi across all VCF transcripts
for(i in 1:nrow(boot_Jsite_pi_sum)) {
  cat(paste(i, " ", sep = ''))
  {boot_Jsite_pi_sum$mean_ref[i] <- boot_Jsite_pi_sum$boot_sum_pi[i]/total_bp_transcripts}
}

#get 95% CIs
Jsite_95ci_pi <- quantile(boot_Jsite_pi_sum$mean_ref, c(0.025, 0.975)) #(0.000823, 0.000867)
Jsite_95ci_pi <- t(data.frame(Jsite_95ci_pi))
rownames(Jsite_95ci_pi) <- c("Japan")

#### P pop ####
#bootstrap across loci to get sum of pi
boot_Psite_pi <- boot(data = Psites_pi$PI, statistic = samp_sum, R = 1000)
boot_Psite_pi_sum <- data.frame(boot_Psite_pi$t)
colnames(boot_Psite_pi_sum) <- c("boot_sum_pi")

#calculate avg pi by site
boot_Psite_pi_sum$mean_ref <- NA #create empty column to populate

#for loop to calculate avg pi across all VCF transcripts
for(i in 1:nrow(boot_Psite_pi_sum)) {
  cat(paste(i, " ", sep = ''))
  {boot_Psite_pi_sum$mean_ref[i] <- boot_Psite_pi_sum$boot_sum_pi[i]/total_bp_transcripts}
}

#get 95% CIs
Psite_95ci_pi <- quantile(boot_Psite_pi_sum$mean_ref, c(0.025, 0.975)) #(0.000945, 0.000984)
Psite_95ci_pi <- t(data.frame(Psite_95ci_pi))
rownames(Psite_95ci_pi) <- c("Philippines")

#### N pop ####
#bootstrap across loci to get sum of pi
boot_Nsite_pi <- boot(data = Nsites_pi$PI, statistic = samp_sum, R = 1000)
boot_Nsite_pi_sum <- data.frame(boot_Nsite_pi$t)
colnames(boot_Nsite_pi_sum) <- c("boot_sum_pi")

#calculate avg pi by site
boot_Nsite_pi_sum$mean_ref <- NA #create empty column to populate

#for loop to calculate avg pi across all VCF transcripts
for(i in 1:nrow(boot_Nsite_pi_sum)) {
  cat(paste(i, " ", sep = ''))
  {boot_Nsite_pi_sum$mean_ref[i] <- boot_Nsite_pi_sum$boot_sum_pi[i]/total_bp_transcripts}
}

#get 95% CIs
Nsite_95ci_pi <- quantile(boot_Nsite_pi_sum$mean_ref, c(0.025, 0.975)) #(0.000930, 0.000976)
Nsite_95ci_pi <- t(data.frame(Nsite_95ci_pi))
rownames(Nsite_95ci_pi) <- c("Indonesia")

######## Create summary tables for site pi ########
#just working with transcript total bp for now

#merge mean and CI data for all pops
pi_site_mean <- as.data.frame(c(mean_Jpi_transcripts, mean_Ppi_transcripts, mean_Npi_transcripts, mean_allpi_transcripts))
rownames(pi_site_mean) <- c("Japan", "Philippines", "Indonesia", "All")
pi_site_95ci <- rbind(Jsite_95ci_pi, Psite_95ci_pi, Nsite_95ci_pi, allsite_95ci_pi)
pi_site_sum <- cbind(pi_site_mean, pi_site_95ci)
colnames(pi_site_sum) <- c("mean", "2.5_per", "97.5_per")
pi_site_sum$Pop <- c("Japan", "Philippines", "Indonesia", "All")

pi_site_sum$diff_lower <- pi_site_sum$mean - pi_site_sum$`2.5_per` #calculate diff btwn sample mean and 2.5 percentile for CI visualization
pi_site_sum$diff_upper <- pi_site_sum$`97.5_per` - pi_site_sum$mean #calculate diff btwn sample mean and 97.5 percentile for CI visualization

#write out data
write.csv(pi_site_sum, "Data/pi_site_cis_mac2.csv")

######## Visualize data ########
#designed to be run separately from earlier sections
#need to sum across transcript

remove(list = ls())

#read in data
pi_site_all <- read.csv("Data/pi_combined_site_full_mac2_aggregate.csv", header = TRUE, row.names = 1)
  #pi_site_all$contig_bp <- paste(pi_site_all$Contig, pi_site_all$POS, sep = "_") #add column to enable addition of outlier status
  colnames(pi_site_all) <- c("CONTIG", "PI", "Contig.length", "pi_avg", "Pop", "NUM")  
pi_site_sum <- read.csv("Data/pi_site_cis_mac2.csv", header = TRUE, row.names = 1)
outliers <- read.csv("Data/outlier_contigs_mac2.csv", header = TRUE)

#add outlier status
setdiff(outliers$CONTIG, pi_site_all$CONTIG) #should = 0, as all should be present in pi_site_all
outlier_loci <- pi_site_all[pi_site_all$CONTIG %in% outliers$CONTIG, ] #df with only outliers
nonoutlier_loci <- pi_site_all[!(pi_site_all$CONTIG %in% outliers$CONTIG), ] #df with all other loci

#set loci status
outlier_loci$Status <- c(rep("Outlier", times = 164))
nonoutlier_loci$Status <- c(rep("Not_Outlier", times = 4336))
pi_site_all$Status <- c(rep("All", times = 4500))

#merge together
pi_site_all <- rbind(outlier_loci, nonoutlier_loci, pi_site_all)

#### Mean pi w/error bars plot ####

#plot of mean w/in pop pi w/95% CI error bars
mean_site_pi_plot <- ggplot(data = pi_site_sum, aes(x = Pop, y = mean)) + 
  geom_point(aes(size = 1), show.legend = FALSE) + 
  geom_errorbar(aes(ymin = mean - diff_lower, ymax = mean + diff_upper, width = 0.5, size = 0.5), show.legend = FALSE)
mean_site_pi_plot_annotated <- mean_site_pi_plot + ggtitle("Mean pi w/95% CI") + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 14, face = "bold"),
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 12, color = "black"))
mean_site_pi_plot_annotated


#### Boxplots ####
#doesn't work quite right bc not divided by sum of bp --> divide by transcript length?
#create new column with transcript length info, then divide pi by that to get per site and plot this
#Still need to do synonymous sites only too (box plot, but calc pi with only them to only look at those segregating sites (not all))

#ordering x-axis
pi_site_all$Pop2 <- factor(pi_site_all$Pop, levels = c("Japan", "Philippines", "Indonesia", "All"))

#boxplot
pi_site_boxplot <- ggplot(data = pi_site_all, aes(x = Pop2, y = pi_avg, color = Status)) + 
  geom_boxplot()
pi_site_boxplot_annotated <- pi_site_boxplot + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(size = 1), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 14, color = "black"), 
        axis.title = element_text(size = 14, face = "bold"), legend.position = "top", 
        legend.text = element_text(size = 12), legend.title = element_text(size = 12))
pi_site_boxplot_annotated

########################################################################################################

######## Site pi calculations (in HWE only) ########
#calculated with mac = 1 & only SNPs in HWE (putatively neutral)
#designed to be run separately from previous sections

remove(list = ls())

#load libraries
library(tidyverse)
library(boot)

#read in data
#site pi all calculated with mac1
inHWE_allsites_pi <- read.csv("Data/inHWE_allmac1_sitespi.csv", header = TRUE) #read in data from vcftools
inHWE_Jsites_pi <- read.csv("Data/inHWE_Jmac1_sitespi.csv", header = TRUE)
inHWE_Nsites_pi <- read.csv("Data/inHWE_Nmac1_sitespi.csv", header = TRUE)
inHWE_Psites_pi <- read.csv("Data/inHWE_Pmac1_sitespi.csv", header = TRUE)

######## Calculate mean pi by site ########

#read in total bp in transcriptome
total_bp_ref <- 45506763 #total bp in New_ref_N3.fa
total_bp_transcripts <- 1163361 #total bp in transcripts included in mac1 VCF

#### All pops ####
sum_inHWE_allpi <- sum(inHWE_allsites_pi$PI) #sum pi across all sites

#average pi across transcriptome
mean_inHWE_allpi_total <- sum_inHWE_allpi/total_bp_ref #0.0000201
mean_inHWE_allpi_transcripts <- sum_inHWE_allpi/total_bp_transcripts #0.000784

#add columns to write out together
inHWE_allsites_pi$Pop <- c(rep("All", times = 5098))
inHWE_allsites_pi$NUM <- c(1:5098)

#### J pop ####
sum_inHWE_Jpi <- sum(inHWE_Jsites_pi$PI)
mean_inHWE_Jpi_total <- sum_inHWE_Jpi/total_bp_ref #0.0000180
mean_inHWE_Jpi_transcripts <- sum_inHWE_Jpi/total_bp_transcripts #0.000705

#add columns to write out together
inHWE_Jsites_pi$Pop <- c(rep("Japan", times = 5098))
inHWE_Jsites_pi$NUM <- c(1:5098)

#### P pop ####
sum_inHWE_Ppi <- sum(inHWE_Psites_pi$PI)
mean_inHWE_Ppi_total <- sum_inHWE_Ppi/total_bp_ref #0.0000203
mean_inHWE_Ppi_transcripts <- sum_inHWE_Ppi/total_bp_transcripts #0.000795

#add columns to write out together
inHWE_Psites_pi$Pop <- c(rep("Philippines", times = 5098))
inHWE_Psites_pi$NUM <- c(1:5098)

#### N pop ####
sum_inHWE_Npi <- sum(inHWE_Nsites_pi$PI)
mean_inHWE_Npi_total <- sum_inHWE_Npi/total_bp_ref #0.0000199
mean_inHWE_Npi_transcripts <- sum_inHWE_Npi/total_bp_transcripts #0.000777

#add columns to write out together
inHWE_Nsites_pi$Pop <- c(rep("Indonesia", times = 5098))
inHWE_Nsites_pi$NUM <- c(1:5098)

#merge dataframes
pi_site_all_inHWE <- rbind(inHWE_Jsites_pi, inHWE_Psites_pi, inHWE_Nsites_pi, inHWE_allsites_pi)
lapply(pi_site_all_inHWE, class) #check character class for columns

#write out combined data
write.csv(pi_site_all_inHWE, "Data/pi_combined_site_full_inHWE_mac1.csv")

######## Bootstrap mean pi by site ########
#just working with transcript total bp for now

#sum function
samp_sum <- function(x, i) {
  sum(x[i])
} #bc if use sum() not summing properly

#### All pops ####
#bootstrap across loci to get sum of pi
boot_inHWE_allsite_pi <- boot(data = inHWE_allsites_pi$PI, statistic = samp_sum, R = 1000)
boot_inHWE_allsite_pi_sum <- data.frame(boot_inHWE_allsite_pi$t)
colnames(boot_inHWE_allsite_pi_sum) <- c("boot_sum_pi")

#calculate avg pi by site
boot_inHWE_allsite_pi_sum$mean_ref <- NA #create empty column to populate

#for loop to calculate avg pi across all VCF transcripts
for(i in 1:nrow(boot_inHWE_allsite_pi_sum)) {
  cat(paste(i, " ", sep = ''))
  {boot_inHWE_allsite_pi_sum$mean_ref[i] <- boot_inHWE_allsite_pi_sum$boot_sum_pi[i]/total_bp_transcripts}
}

#get 95% CIs
inHWE_allsite_95ci_pi <- quantile(boot_inHWE_allsite_pi_sum$mean_ref, c(0.025, 0.975)) #(0.000765, 0.000801)
inHWE_allsite_95ci_pi <- t(data.frame(inHWE_allsite_95ci_pi))
rownames(inHWE_allsite_95ci_pi) <- c("All")

#### J pop ####
#bootstrap across loci to get sum of pi
boot_inHWE_Jsite_pi <- boot(data = inHWE_Jsites_pi$PI, statistic = samp_sum, R = 1000)
boot_inHWE_Jsite_pi_sum <- data.frame(boot_inHWE_Jsite_pi$t)
colnames(boot_inHWE_Jsite_pi_sum) <- c("boot_sum_pi")

#calculate avg pi by site
boot_inHWE_Jsite_pi_sum$mean_ref <- NA #create empty column to populate

#for loop to calculate avg pi across all VCF transcripts
for(i in 1:nrow(boot_inHWE_Jsite_pi_sum)) {
  cat(paste(i, " ", sep = ''))
  {boot_inHWE_Jsite_pi_sum$mean_ref[i] <- boot_inHWE_Jsite_pi_sum$boot_sum_pi[i]/total_bp_transcripts}
}

#get 95% CIs
inHWE_Jsite_95ci_pi <- quantile(boot_inHWE_Jsite_pi_sum$mean_ref, c(0.025, 0.975)) #(0.000684, 0.000726)
inHWE_Jsite_95ci_pi <- t(data.frame(inHWE_Jsite_95ci_pi))
rownames(inHWE_Jsite_95ci_pi) <- c("Japan")

#### P pop ####
#bootstrap across loci to get sum of pi
boot_inHWE_Psite_pi <- boot(data = inHWE_Psites_pi$PI, statistic = samp_sum, R = 1000)
boot_inHWE_Psite_pi_sum <- data.frame(boot_inHWE_Psite_pi$t)
colnames(boot_inHWE_Psite_pi_sum) <- c("boot_sum_pi")

#calculate avg pi by site
boot_inHWE_Psite_pi_sum$mean_ref <- NA #create empty column to populate

#for loop to calculate avg pi across all VCF transcripts
for(i in 1:nrow(boot_inHWE_Psite_pi_sum)) {
  cat(paste(i, " ", sep = ''))
  {boot_inHWE_Psite_pi_sum$mean_ref[i] <- boot_inHWE_Psite_pi_sum$boot_sum_pi[i]/total_bp_transcripts}
}

#get 95% CIs
inHWE_Psite_95ci_pi <- quantile(boot_inHWE_Psite_pi_sum$mean_ref, c(0.025, 0.975)) #(0.000776, 0.000815)
inHWE_Psite_95ci_pi <- t(data.frame(inHWE_Psite_95ci_pi))
rownames(inHWE_Psite_95ci_pi) <- c("Philippines")

#### N pop ####
#bootstrap across loci to get sum of pi
boot_inHWE_Nsite_pi <- boot(data = inHWE_Nsites_pi$PI, statistic = samp_sum, R = 1000)
boot_inHWE_Nsite_pi_sum <- data.frame(boot_inHWE_Nsite_pi$t)
colnames(boot_inHWE_Nsite_pi_sum) <- c("boot_sum_pi")

#calculate avg pi by site
boot_inHWE_Nsite_pi_sum$mean_ref <- NA #create empty column to populate

#for loop to calculate avg pi across all VCF transcripts
for(i in 1:nrow(boot_inHWE_Nsite_pi_sum)) {
  cat(paste(i, " ", sep = ''))
  {boot_inHWE_Nsite_pi_sum$mean_ref[i] <- boot_inHWE_Nsite_pi_sum$boot_sum_pi[i]/total_bp_transcripts}
}

#get 95% CIs
inHWE_Nsite_95ci_pi <- quantile(boot_inHWE_Nsite_pi_sum$mean_ref, c(0.025, 0.975)) #(0.000755, 0.000798)
inHWE_Nsite_95ci_pi <- t(data.frame(inHWE_Nsite_95ci_pi))
rownames(inHWE_Nsite_95ci_pi) <- c("Indonesia")

######## Create summary tables for site pi ########
#just working with transcript total bp for now

#merge mean and CI data for all pops
inHWE_pi_site_mean <- as.data.frame(c(mean_inHWE_Jpi_transcripts, mean_inHWE_Ppi_transcripts, mean_inHWE_Npi_transcripts, mean_inHWE_allpi_transcripts))
rownames(inHWE_pi_site_mean) <- c("Japan", "Philippines", "Indonesia", "All")
inHWE_pi_site_95ci <- rbind(inHWE_Jsite_95ci_pi, inHWE_Psite_95ci_pi, inHWE_Nsite_95ci_pi, inHWE_allsite_95ci_pi)
inHWE_pi_site_sum <- cbind(inHWE_pi_site_mean, inHWE_pi_site_95ci)
colnames(inHWE_pi_site_sum) <- c("mean", "2.5_per", "97.5_per")
inHWE_pi_site_sum$Pop <- c("Japan", "Philippines", "Indonesia", "All")

inHWE_pi_site_sum$diff_lower <- inHWE_pi_site_sum$mean - inHWE_pi_site_sum$`2.5_per` #calculate diff btwn sample mean and 2.5 percentile for CI visualization
inHWE_pi_site_sum$diff_upper <- inHWE_pi_site_sum$`97.5_per` - inHWE_pi_site_sum$mean #calculate diff btwn sample mean and 97.5 percentile for CI visualization

#write out data
write.csv(inHWE_pi_site_sum, "Data/pi_site_cis_mac1_inHWE.csv")

######## Visualize data ########
#designed to be run separately from earlier sections

remove(list = ls())

#read in data
inHWE_pi_site_sum <- read.csv("Data/pi_site_cis_mac1_inHWE.csv", header = TRUE, row.names = 1)

#### Mean pi w/error bars plot ####

#plot of mean w/in pop pi w/95% CI error bars
mean_inHWE_site_pi_plot <- ggplot(data = inHWE_pi_site_sum, aes(x = Pop, y = mean)) + 
  geom_point(aes(size = 1), show.legend = FALSE) + 
  geom_errorbar(aes(ymin = mean - diff_lower, ymax = mean + diff_upper, width = 0.5, size = 0.5), show.legend = FALSE)
mean_inHWE_site_pi_plot_annotated <- mean_inHWE_site_pi_plot + ggtitle("Mean pi w/95% CI") + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 14, face = "bold"),
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 12, color = "black"))
mean_inHWE_site_pi_plot_annotated

########################################################################################################

######## Site pi calculations (synonymous sites only) ########
#calculated with mac = 1 & only synonymous sites
#designed to be run separately from previous sections

remove(list = ls())

#load libraries
library(tidyverse)
library(boot)

#read in data
#site pi all calculated with mac1
SYN_allsites_pi <- read.csv("Data/SYN_mac1_Allpi.csv", header = TRUE) #read in data from vcftools
SYN_Jsites_pi <- read.csv("Data/SYN_mac1_Jpi.csv", header = TRUE)
SYN_Nsites_pi <- read.csv("Data/SYN_mac1_Npi.csv", header = TRUE)
SYN_Psites_pi <- read.csv("Data/SYN_mac1_Ppi.csv", header = TRUE)
contig_length <- read.csv("Data/Contig_length.csv", header = TRUE)
  colnames(contig_length) <- c("CHROM", "Contig.length")

######## Calculate mean pi by site ########

#read in total bp in transcriptome
total_bp_ref <- 45506763 #total bp in New_ref_N3.fa
total_bp_transcripts <- 534055 #total bp in transcripts included in mac2 VCF

#### All pops ####
sum_allpi <- sum(SYN_allsites_pi$PI) #sum pi across all sites

#average pi across transcriptome
mean_allpi_total <- sum_allpi/total_bp_ref #0.0000235
mean_allpi_transcripts <- sum_allpi/total_bp_transcripts #0.000955

#add columns to write out together
SYN_allsites_pi$Pop <- c(rep("All", times = 1453))
SYN_allsites_pi$NUM <- c(1:1453)

#calculate pi by transcript
allsites_pi_aggregate <- aggregate(SYN_allsites_pi$PI~SYN_allsites_pi$CHROM, FUN = sum) #sum pi by transcript
colnames(allsites_pi_aggregate) <- c("CHROM", "PI")

allsites_pi_aggregate <- merge(allsites_pi_aggregate, contig_length) #merge with contig length
allsites_pi_aggregate$pi_avg <- NA #create empty column to populate

for(i in 1:nrow(allsites_pi_aggregate)) {
  cat(paste(i, " ", sep = ''))
  {allsites_pi_aggregate$pi_avg[i] <- allsites_pi_aggregate$PI[i]/allsites_pi_aggregate$Contig.length[i]}
}

allsites_pi_aggregate$Pop <- c(rep("All", times = 498))
allsites_pi_aggregate$NUM <- c(1:498)

#### J pop ####
sum_Jpi <- sum(SYN_Jsites_pi$PI)
mean_Jpi_total <- sum_Jpi/total_bp_ref #0.0000198
mean_Jpi_transcripts <- sum_Jpi/total_bp_transcripts #0.000844

#add columns to write out together
SYN_Jsites_pi$Pop <- c(rep("Japan", times = 1453))
SYN_Jsites_pi$NUM <- c(1:1453)

#calculate pi by transcript
Jsites_pi_aggregate <- aggregate(SYN_Jsites_pi$PI~SYN_Jsites_pi$CHROM, FUN = sum) #sum pi by transcript
colnames(Jsites_pi_aggregate) <- c("CHROM", "PI")

Jsites_pi_aggregate <- merge(Jsites_pi_aggregate, contig_length) #merge with contig length
Jsites_pi_aggregate$pi_avg <- NA #create empty column to populate

for(i in 1:nrow(Jsites_pi_aggregate)) {
  cat(paste(i, " ", sep = ''))
  {Jsites_pi_aggregate$pi_avg[i] <- Jsites_pi_aggregate$PI[i]/Jsites_pi_aggregate$Contig.length[i]}
}

Jsites_pi_aggregate$Pop <- c(rep("Japan", times = 498))
Jsites_pi_aggregate$NUM <- c(1:498)

#### P pop ####
sum_Ppi <- sum(SYN_Psites_pi$PI)
mean_Ppi_total <- sum_Ppi/total_bp_ref #0.0000225
mean_Ppi_transcripts <- sum_Ppi/total_bp_transcripts #0.000964

#add columns to write out together
SYN_Psites_pi$Pop <- c(rep("Philippines", times = 1453))
SYN_Psites_pi$NUM <- c(1:1453)

#calculate pi by transcript
Psites_pi_aggregate <- aggregate(SYN_Psites_pi$PI~SYN_Psites_pi$CHROM, FUN = sum) #sum pi by transcript
colnames(Psites_pi_aggregate) <- c("CHROM", "PI")

Psites_pi_aggregate <- merge(Psites_pi_aggregate, contig_length) #merge with contig length
Psites_pi_aggregate$pi_avg <- NA #create empty column to populate

for(i in 1:nrow(Psites_pi_aggregate)) {
  cat(paste(i, " ", sep = ''))
  {Psites_pi_aggregate$pi_avg[i] <- Psites_pi_aggregate$PI[i]/Psites_pi_aggregate$Contig.length[i]}
}

Psites_pi_aggregate$Pop <- c(rep("Philippines", times = 498))
Psites_pi_aggregate$NUM <- c(1:498)

#### N pop ####
sum_Npi <- sum(SYN_Nsites_pi$PI)
mean_Npi_total <- sum_Npi/total_bp_ref #0.0000223
mean_Npi_transcripts <- sum_Npi/total_bp_transcripts #0.000953

#add columns to write out together
SYN_Nsites_pi$Pop <- c(rep("Indonesia", times = 1453))
SYN_Nsites_pi$NUM <- c(1:1453)

#calculate pi by transcript
Nsites_pi_aggregate <- aggregate(SYN_Nsites_pi$PI~SYN_Nsites_pi$CHROM, FUN = sum) #sum pi by transcript
colnames(Nsites_pi_aggregate) <- c("CHROM", "PI")

Nsites_pi_aggregate <- merge(Nsites_pi_aggregate, contig_length) #merge with contig length
Nsites_pi_aggregate$pi_avg <- NA #create empty column to populate

for(i in 1:nrow(Nsites_pi_aggregate)) {
  cat(paste(i, " ", sep = ''))
  {Nsites_pi_aggregate$pi_avg[i] <- Nsites_pi_aggregate$PI[i]/Nsites_pi_aggregate$Contig.length[i]}
}

Nsites_pi_aggregate$Pop <- c(rep("Indonesia", times = 498))
Nsites_pi_aggregate$NUM <- c(1:498)

#merge dataframes
SYN_pi_site_all <- rbind(SYN_Jsites_pi, SYN_Psites_pi, SYN_Nsites_pi, SYN_allsites_pi)
lapply(SYN_pi_site_all, class) #check character class for columns

pi_site_all_aggregate <- rbind(Jsites_pi_aggregate, Psites_pi_aggregate, Nsites_pi_aggregate, allsites_pi_aggregate)
lapply(pi_site_all_aggregate, class)

#write out combined data
write.csv(SYN_pi_site_all, "Data/pi_combined_site_mac1_SYN_full.csv")
write.csv(pi_site_all_aggregate, "Data/pi_combined_site_full_mac1_SYN_aggregate.csv")

######## Bootstrap mean pi by site ########
#just working with transcript total bp for now

#sum function
samp_sum <- function(x, i) {
  sum(x[i])
} #bc if use sum() not summing properly

#### All pops ####
#bootstrap across loci to get sum of pi
boot_allsite_pi <- boot(data = allsites_pi$PI, statistic = samp_sum, R = 1000)
boot_allsite_pi_sum <- data.frame(boot_allsite_pi$t)
colnames(boot_allsite_pi_sum) <- c("boot_sum_pi")

#calculate avg pi by site
boot_allsite_pi_sum$mean_ref <- NA #create empty column to populate

#for loop to calculate avg pi across all VCF transcripts
for(i in 1:nrow(boot_allsite_pi_sum)) {
  cat(paste(i, " ", sep = ''))
  {boot_allsite_pi_sum$mean_ref[i] <- boot_allsite_pi_sum$boot_sum_pi[i]/total_bp_transcripts}
}

#get 95% CIs
allsite_95ci_pi <- quantile(boot_allsite_pi_sum$mean_ref, c(0.025, 0.975)) #(0.000938, 0.000974)
allsite_95ci_pi <- t(data.frame(allsite_95ci_pi))
rownames(allsite_95ci_pi) <- c("All")

#### J pop ####
#bootstrap across loci to get sum of pi
boot_Jsite_pi <- boot(data = Jsites_pi$PI, statistic = samp_sum, R = 1000)
boot_Jsite_pi_sum <- data.frame(boot_Jsite_pi$t)
colnames(boot_Jsite_pi_sum) <- c("boot_sum_pi")

#calculate avg pi by site
boot_Jsite_pi_sum$mean_ref <- NA #create empty column to populate

#for loop to calculate avg pi across all VCF transcripts
for(i in 1:nrow(boot_Jsite_pi_sum)) {
  cat(paste(i, " ", sep = ''))
  {boot_Jsite_pi_sum$mean_ref[i] <- boot_Jsite_pi_sum$boot_sum_pi[i]/total_bp_transcripts}
}

#get 95% CIs
Jsite_95ci_pi <- quantile(boot_Jsite_pi_sum$mean_ref, c(0.025, 0.975)) #(0.000823, 0.000867)
Jsite_95ci_pi <- t(data.frame(Jsite_95ci_pi))
rownames(Jsite_95ci_pi) <- c("Japan")

#### P pop ####
#bootstrap across loci to get sum of pi
boot_Psite_pi <- boot(data = Psites_pi$PI, statistic = samp_sum, R = 1000)
boot_Psite_pi_sum <- data.frame(boot_Psite_pi$t)
colnames(boot_Psite_pi_sum) <- c("boot_sum_pi")

#calculate avg pi by site
boot_Psite_pi_sum$mean_ref <- NA #create empty column to populate

#for loop to calculate avg pi across all VCF transcripts
for(i in 1:nrow(boot_Psite_pi_sum)) {
  cat(paste(i, " ", sep = ''))
  {boot_Psite_pi_sum$mean_ref[i] <- boot_Psite_pi_sum$boot_sum_pi[i]/total_bp_transcripts}
}

#get 95% CIs
Psite_95ci_pi <- quantile(boot_Psite_pi_sum$mean_ref, c(0.025, 0.975)) #(0.000945, 0.000984)
Psite_95ci_pi <- t(data.frame(Psite_95ci_pi))
rownames(Psite_95ci_pi) <- c("Philippines")

#### N pop ####
#bootstrap across loci to get sum of pi
boot_Nsite_pi <- boot(data = Nsites_pi$PI, statistic = samp_sum, R = 1000)
boot_Nsite_pi_sum <- data.frame(boot_Nsite_pi$t)
colnames(boot_Nsite_pi_sum) <- c("boot_sum_pi")

#calculate avg pi by site
boot_Nsite_pi_sum$mean_ref <- NA #create empty column to populate

#for loop to calculate avg pi across all VCF transcripts
for(i in 1:nrow(boot_Nsite_pi_sum)) {
  cat(paste(i, " ", sep = ''))
  {boot_Nsite_pi_sum$mean_ref[i] <- boot_Nsite_pi_sum$boot_sum_pi[i]/total_bp_transcripts}
}

#get 95% CIs
Nsite_95ci_pi <- quantile(boot_Nsite_pi_sum$mean_ref, c(0.025, 0.975)) #(0.000930, 0.000976)
Nsite_95ci_pi <- t(data.frame(Nsite_95ci_pi))
rownames(Nsite_95ci_pi) <- c("Indonesia")

######## Create summary tables for site pi ########
#just working with transcript total bp for now

#merge mean and CI data for all pops
pi_site_mean <- as.data.frame(c(mean_Jpi_transcripts, mean_Ppi_transcripts, mean_Npi_transcripts, mean_allpi_transcripts))
rownames(pi_site_mean) <- c("Japan", "Philippines", "Indonesia", "All")
pi_site_95ci <- rbind(Jsite_95ci_pi, Psite_95ci_pi, Nsite_95ci_pi, allsite_95ci_pi)
pi_site_sum <- cbind(pi_site_mean, pi_site_95ci)
colnames(pi_site_sum) <- c("mean", "2.5_per", "97.5_per")
pi_site_sum$Pop <- c("Japan", "Philippines", "Indonesia", "All")

pi_site_sum$diff_lower <- pi_site_sum$mean - pi_site_sum$`2.5_per` #calculate diff btwn sample mean and 2.5 percentile for CI visualization
pi_site_sum$diff_upper <- pi_site_sum$`97.5_per` - pi_site_sum$mean #calculate diff btwn sample mean and 97.5 percentile for CI visualization

#write out data
write.csv(pi_site_sum, "Data/pi_site_cis_mac2.csv")

######## Visualize data ########
#designed to be run separately from earlier sections
#need to sum across transcript

remove(list = ls())

#read in data
pi_site_all <- read.csv("Data/pi_combined_site_full_mac1_SYN_aggregate.csv", header = TRUE, row.names = 1)
#pi_site_all$contig_bp <- paste(pi_site_all$Contig, pi_site_all$POS, sep = "_") #add column to enable addition of outlier status
colnames(pi_site_all) <- c("CONTIG", "Contig.length", "pi_avg", "Pop", "NUM")  
#pi_site_sum <- read.csv("Data/pi_site_cis_mac2.csv", header = TRUE, row.names = 1)
outliers <- read.csv("Data/outlier_sequences_pi_mac1_SYN.csv", header = TRUE)

#add outlier status
setdiff(outliers$CONTIG, pi_site_all$CONTIG) #should = 0, as all should be present in pi_site_all
outlier_loci <- pi_site_all[pi_site_all$CONTIG %in% outliers$CONTIG, ] #df with only outliers
nonoutlier_loci <- pi_site_all[!(pi_site_all$CONTIG %in% outliers$CONTIG), ] #df with all other loci

#set loci status
outlier_loci$Status <- c(rep("Outlier", times = 112))
nonoutlier_loci$Status <- c(rep("Not_Outlier", times = 1880))
pi_site_all$Status <- c(rep("All", times = 1992))

#merge together
pi_site_all <- rbind(outlier_loci, nonoutlier_loci, pi_site_all)

#### Mean pi w/error bars plot ####

#plot of mean w/in pop pi w/95% CI error bars
mean_site_pi_plot <- ggplot(data = pi_site_sum, aes(x = Pop, y = mean)) + 
  geom_point(aes(size = 1), show.legend = FALSE) + 
  geom_errorbar(aes(ymin = mean - diff_lower, ymax = mean + diff_upper, width = 0.5, size = 0.5), show.legend = FALSE)
mean_site_pi_plot_annotated <- mean_site_pi_plot + ggtitle("Mean pi w/95% CI") + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 14, face = "bold"),
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 12, color = "black"))
mean_site_pi_plot_annotated

#### Boxplots ####

#ordering x-axis
pi_site_all$Pop2 <- factor(pi_site_all$Pop, levels = c("Japan", "Philippines", "Indonesia", "All"))

#boxplot
pi_site_boxplot <- ggplot(data = pi_site_all, aes(x = Pop2, y = pi_avg, color = Status)) + 
  geom_boxplot(outlier.size = 6, lwd = 3)
pi_site_boxplot_annotated <- pi_site_boxplot + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(size = 2), 
        axis.ticks = element_line(color = "black", size = 2), 
        axis.text = element_text(size = 34, color = "black"), 
        axis.title = element_text(size = 34, face = "bold"), legend.position = "top", 
        legend.text = element_text(size = 32), legend.title = element_text(size = 32))
pi_site_boxplot_annotated

fst <- read.csv("Data/Fst_perloc_mac1.csv")
all_pi <- read.csv("Data/SYN_mac1_Allpi.csv")
all_pi$X <- paste(all_pi$CHROM, all_pi$POS, sep = "_")

pi_fst_all <- merge(all_pi, fst, by = "X")

#boxplot
pi_fst_plot <- ggplot(data = pi_fst_all, aes(x = PI, y = Fst, color = status)) + 
  geom_point(size = 8) + geom_hline(yintercept = 0, color = "black", size = 3)
pi_fst_plot_annotated <- pi_fst_plot + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(size = 2), 
        axis.ticks = element_line(color = "black", size = 2), 
        axis.text = element_text(size = 34, color = "black"), 
        axis.title = element_text(size = 34, face = "bold"), legend.position = "top", 
        legend.text = element_text(size = 32), legend.title = element_text(size = 32))
pi_fst_plot_annotated