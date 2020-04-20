################################################### Script for Pi  #######################################################

#Created for transcriptome project
#Uses window pi estimates from VCFtools and bootstraps to get 95% CIs 

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
outliers$Status <- c(rep("Outlier", times = 292))
nonoutliers$Status <- c(rep("Not_Outlier", times = 8720))
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

TD_boxplot <- ggplot(data = TD_only_all, aes(x = Pop2, y = TajimaD, color = Status)) + 
  geom_boxplot()
TD_boxplot_annotated <- TD_boxplot + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(size = 1), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 14, face = "bold"), legend.position = "top",
        legend.text = element_text(size = 12), legend.title = element_text(size = 12))
TD_boxplot_annotated
