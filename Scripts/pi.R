################################################### Script for Pi  #######################################################

#Created for transcriptome project
#Uses site pi estimates from VCFtools and bootstraps to get 95% CIs 

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
allsites_pi <- read.csv("Data/allmac2_sitespi.csv", header = TRUE) #read in data from vcftools
Jsites_pi <- read.csv("Data/Jmac2_sitespi.csv", header = TRUE)
Nsites_pi <- read.csv("Data/Nmac2_sitespi.csv", header = TRUE)
Psites_pi <- read.csv("Data/Pmac2_sitespi.csv", header = TRUE)
contig_length <- read.csv("Data/Contig_length.csv", header = TRUE) #to calculate pi by transcript (need transcript bp length)
  colnames(contig_length) <- c("CHROM", "Contig.length")

#################################################################################################################################################

######## Site pi calculations ########
#calculated with mac = 2

######## Calculate mean pi by site ########

#read in total bp in transcriptome
total_bp_transcripts <- 1064586 #total bp in transcripts included in mac2 VCF

#### All pops ####
sum_allpi <- sum(allsites_pi$PI) #sum pi across all sites

#average pi across transcriptome
mean_allpi_transcripts <- sum_allpi/total_bp_transcripts #0.000955

#add columns to write out together
allsites_pi$Pop <- c(rep("All", times = 4212))
allsites_pi$NUM <- c(1:4212)

#sum pi by transcript
allsites_pi_aggregate <- aggregate(allsites_pi$PI~allsites_pi$CHROM, FUN = sum) #sum pi by transcript
  colnames(allsites_pi_aggregate) <- c("CHROM", "PI")

allsites_pi_aggregate <- merge(allsites_pi_aggregate, contig_length) #merge with contig length
allsites_pi_aggregate$pi_avg <- NA #create empty column to populate
setdiff(allsites_pi_aggregate$CHROM, contig_length$CHROM) #make sure all transcripts have length

#calculate pi by transcript
for(i in 1:nrow(allsites_pi_aggregate)) {
  cat(paste(i, " ", sep = ''))
  {allsites_pi_aggregate$pi_avg[i] <- allsites_pi_aggregate$PI[i]/allsites_pi_aggregate$Contig.length[i]}
}

allsites_pi_aggregate$Pop <- c(rep("All", times = 1002))
allsites_pi_aggregate$NUM <- c(1:1002)

#### J pop ####
sum_Jpi <- sum(Jsites_pi$PI)
mean_Jpi_transcripts <- sum_Jpi/total_bp_transcripts #0.000844

#add columns to write out together
Jsites_pi$Pop <- c(rep("Japan", times = 4212))
Jsites_pi$NUM <- c(1:4212)

#sum pi by transcript
Jsites_pi_aggregate <- aggregate(Jsites_pi$PI~Jsites_pi$CHROM, FUN = sum) #sum pi by transcript
  colnames(Jsites_pi_aggregate) <- c("CHROM", "PI")

Jsites_pi_aggregate <- merge(Jsites_pi_aggregate, contig_length) #merge with contig length
Jsites_pi_aggregate$pi_avg <- NA #create empty column to populate

#calculate pi by transcript
for(i in 1:nrow(Jsites_pi_aggregate)) {
  cat(paste(i, " ", sep = ''))
  {Jsites_pi_aggregate$pi_avg[i] <- Jsites_pi_aggregate$PI[i]/Jsites_pi_aggregate$Contig.length[i]}
}

Jsites_pi_aggregate$Pop <- c(rep("Japan", times = 1002))
Jsites_pi_aggregate$NUM <- c(1:1002)

#### P pop ####
sum_Ppi <- sum(Psites_pi$PI)
mean_Ppi_transcripts <- sum_Ppi/total_bp_transcripts #0.000964

#add columns to write out together
Psites_pi$Pop <- c(rep("Philippines", times = 4212))
Psites_pi$NUM <- c(1:4212)

#sum pi by transcript
Psites_pi_aggregate <- aggregate(Psites_pi$PI~Psites_pi$CHROM, FUN = sum) #sum pi by transcript
  colnames(Psites_pi_aggregate) <- c("CHROM", "PI")

Psites_pi_aggregate <- merge(Psites_pi_aggregate, contig_length) #merge with contig length
Psites_pi_aggregate$pi_avg <- NA #create empty column to populate

#calculate pi by transcript
for(i in 1:nrow(Psites_pi_aggregate)) {
  cat(paste(i, " ", sep = ''))
  {Psites_pi_aggregate$pi_avg[i] <- Psites_pi_aggregate$PI[i]/Psites_pi_aggregate$Contig.length[i]}
}

Psites_pi_aggregate$Pop <- c(rep("Philippines", times = 1002))
Psites_pi_aggregate$NUM <- c(1:1002)

#### N pop ####
sum_Npi <- sum(Nsites_pi$PI)
mean_Npi_transcripts <- sum_Npi/total_bp_transcripts #0.000953

#add columns to write out together
Nsites_pi$Pop <- c(rep("Indonesia", times = 4212))
Nsites_pi$NUM <- c(1:4212)

#sum pi by transcript
Nsites_pi_aggregate <- aggregate(Nsites_pi$PI~Nsites_pi$CHROM, FUN = sum) #sum pi by transcript
  colnames(Nsites_pi_aggregate) <- c("CHROM", "PI")

Nsites_pi_aggregate <- merge(Nsites_pi_aggregate, contig_length) #merge with contig length
Nsites_pi_aggregate$pi_avg <- NA #create empty column to populate

#calculate pi by transcript
for(i in 1:nrow(Nsites_pi_aggregate)) {
  cat(paste(i, " ", sep = ''))
  {Nsites_pi_aggregate$pi_avg[i] <- Nsites_pi_aggregate$PI[i]/Nsites_pi_aggregate$Contig.length[i]}
}

Nsites_pi_aggregate$Pop <- c(rep("Indonesia", times = 1002))
Nsites_pi_aggregate$NUM <- c(1:1002)

#### Merge dataframes ####

#per locus aggregate dataframe
pi_site_all <- rbind(Jsites_pi, Psites_pi, Nsites_pi, allsites_pi) 
  lapply(pi_site_all, class) #check character class for columns

#per transcript aggregate dataframe
pi_site_all_aggregate <- rbind(Jsites_pi_aggregate, Psites_pi_aggregate, Nsites_pi_aggregate, allsites_pi_aggregate) 
  lapply(pi_site_all_aggregate, class) #check character class for columns

#write out combined data
write.csv(pi_site_all, "Data/pi_combined_site_mac2_full.csv")
write.csv(pi_site_all_aggregate, "Data/pi_combined_site_full_mac2_aggregate.csv")

#################################################################################################################################################

######## Bootstrap mean pi by site ########

#sum function
samp_sum <- function(x, i) {
  sum(x[i])
} #bc if use sum() not summing properly

#### All pops ####
#bootstrap across loci to get sum of pi
boot_allsite_pi <- boot(data = allsites_pi$PI, statistic = samp_sum, R = 1000)
boot_allsite_pi_sum <- data.frame(boot_allsite_pi$t) #pull bootstrapped sum pi statistic
  colnames(boot_allsite_pi_sum) <- c("boot_sum_pi")

#calculate avg pi by site
boot_allsite_pi_sum$mean_ref <- NA #create empty column to populate

#for loop to calculate avg pi across all VCF transcripts
for(i in 1:nrow(boot_allsite_pi_sum)) {
  cat(paste(i, " ", sep = ''))
  {boot_allsite_pi_sum$mean_ref[i] <- boot_allsite_pi_sum$boot_sum_pi[i]/total_bp_transcripts}
}

#get 95% CIs
allsite_95ci_pi <- quantile(boot_allsite_pi_sum$mean_ref, c(0.025, 0.975)) #(0.000937, 0.000972)
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
Jsite_95ci_pi <- quantile(boot_Jsite_pi_sum$mean_ref, c(0.025, 0.975)) #(0.000822, 0.000867)
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
Nsite_95ci_pi <- quantile(boot_Nsite_pi_sum$mean_ref, c(0.025, 0.975)) #(0.000930, 0.000973)
Nsite_95ci_pi <- t(data.frame(Nsite_95ci_pi))
  rownames(Nsite_95ci_pi) <- c("Indonesia")

######## Create summary tables for site pi ########

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

#################################################################################################################################################

######## Visualize data ########
#designed to be run separately from earlier sections

remove(list = ls())

#read in data
pi_site_all <- read.csv("Data/pi_combined_site_full_mac2_aggregate.csv", header = TRUE, row.names = 1)
    colnames(pi_site_all) <- c("CONTIG", "PI", "Contig.length", "pi_avg", "Pop", "NUM")  
pi_site_sum <- read.csv("Data/pi_site_cis_mac2.csv", header = TRUE, row.names = 1)
outliers <- read.csv("Data/outlier_contigs_mac2.csv", header = TRUE)

#add outlier status
setdiff(outliers$CONTIG, pi_site_all$CONTIG) #should = 0, as all should be present in pi_site_all
outlier_loci <- pi_site_all[pi_site_all$CONTIG %in% outliers$CONTIG, ] #df with only outliers
nonoutlier_loci <- pi_site_all[!(pi_site_all$CONTIG %in% outliers$CONTIG), ] #df with all other loci

#set loci status
outlier_loci$Status <- c(rep("Outlier", times = 172))
nonoutlier_loci$Status <- c(rep("Not_Outlier", times = 3836))
pi_site_all$Status <- c(rep("All", times = 4008))

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
#not per-transcript calculations, just per-site

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