################################################### Script for Fst  #######################################################

#Created for transcriptome project
#Calculates Fst (pairwise and per locus) using the hierfstat package

#################################################################################################################################################

######## Set-up ########

#set working directory
setwd("C:/Users/rclar/Dropbox/Pinsky_Lab/Transcriptome_Proj/R_scripts/A_clarkii_transcriptomics/")
getwd()

remove(list = ls())

#load libraries
library(tidyverse)
library(vcfR)
library(adegenet)
library(hierfstat)

#read in data
mac2_vcf <- read.vcfR("../../VCFs_and_PLINK/output.hicov2.snps.only.mac2.vcf", verbose = FALSE)
  mac2_genind <- vcfR2genind(mac2_vcf) #convert to genind object for analyses
locnames <- read_table2("Data/Loc_Names_mac2.txt", col_names = TRUE) #read in contig & bp for loci passing mac2 filter
outlierloc <- read.csv("Data/outlier_loci_mac2.csv", header = TRUE)
mac1_vcf <- read.vcfR("../../VCFs_and_PLINK/output.hicov2.snps.only.mac1.vcf", verbose = FALSE)
  mac1_genind <- vcfR2genind(mac1_vcf)
locnames_mac1 <- read_table2("Data/Loc_names_mac1.txt", col_names = FALSE)
  colnames(locnames_mac1) <- c("Contig_bp")
outlierloc_mac1 <- read.csv("Data/outlier_loci_mac1.csv", header = TRUE)

#################################################################################################################################################

######## Calculate Fst per locus --> mac>2 ########

#add pop levels to genind objects
Pop <- c(rep(1, times = 8), rep(2, times = 7), rep(3, times = 10)) #1 = Japan, 2 = Indonesia, 3 = Philippines
pop(mac2_genind) <- Pop
mac2_genind #check to make sure three pops 

#calculate Fst per locus
mac2_stats <- basic.stats(mac2_genind)
stats_perloc <- data.frame(mac2_stats$perloc)

#can bootstrap with Het datasets if want

######## Calculate pairwise-Fst ########

mac2_hierf <- genind2hierfstat(mac2_genind) #convert to hierfstat db for pairwise analyses
pairwise_fst <- genet.dist(mac2_hierf, method = "WC84") #calculates Weir & Cockerham's Fst

#bootstrap pairwise_fst for 95% CI

#need to convert pop character to numeric for bootstrap to work
mac2_hierf$pop <- as.numeric(mac2_hierf$pop)
class(mac2_hierf$pop) #check to make sure numeric

#bootstrap
pairwise_boot <- boot.ppfst(dat = mac2_hierf, nboot = 1000, quant = c(0.025, 0.975), diploid = TRUE)

#get 95% CI limits
ci_upper <- pairwise_boot$ul
ci_lower <- pairwise_boot$ll

######## Add outlier status to stats_perloc df ########

#reorder Fst stats df
stats_perloc <- cbind(stats_perloc, locnames)
stats_perloc <- stats_perloc[order(stats_perloc$Contig_bp), ]
stats_perloc$NUM <- c(1:4212)

#pull outlier loci data
setdiff(outlierloc$contig_bp, stats_perloc$Contig_bp) #should = 0, as all should be present in stats_perloc
outlier_stats <- stats_perloc[stats_perloc$Contig_bp %in% outlierloc$contig_bp, ] #df with only outliers
nonoutlier_stats <- stats_perloc[!(stats_perloc$Contig_bp %in% outlierloc$contig_bp), ] #df with all other loci

#set outlier status
outlier_stats$status <- c(rep("Outlier", times = 82))
nonoutlier_stats$status <- c(rep("Not_Outlier", times = 4130))

#merge together
stats_perloc <- rbind(outlier_stats, nonoutlier_stats)
stats_perloc <-stats_perloc[order(stats_perloc$NUM), ]

#write out
write.csv(stats_perloc, "Data/Fst_perloc.csv")

#################################################################################################################################################

######## Calculate Fst per locus --> mac>1 ########

#add pop levels to genind objects
Pop <- c(rep(1, times = 8), rep(2, times = 7), rep(3, times = 10)) #1 = Japan, 2 = Indonesia, 3 = Philippines
pop(mac1_genind) <- Pop
mac1_genind #check to make sure three pops 

#calculate Fst per locus
mac1_stats <- basic.stats(mac1_genind)
stats_perloc_mac1 <- data.frame(mac1_stats$perloc)

#can bootstrap with Het datasets if want

######## Calculate pairwise-Fst ########

mac1_hierf <- genind2hierfstat(mac1_genind) #convert to hierfstat db for pairwise analyses
pairwise_fst <- genet.dist(mac1_hierf, method = "WC84") #calculates Weir & Cockerham's Fst

#bootstrap pairwise_fst for 95% CI

#need to convert pop character to numeric for bootstrap to work
mac1_hierf$pop <- as.numeric(mac1_hierf$pop)
class(mac1_hierf$pop) #check to make sure numeric

#bootstrap
pairwise_boot_mac1 <- boot.ppfst(dat = mac1_hierf, nboot = 1000, quant = c(0.025, 0.975), diploid = TRUE)

#get 95% CI limits
ci_upper_mac1 <- pairwise_boot_mac1$ul
ci_lower_mac1 <- pairwise_boot_mac1$ll

######## Add outlier status to stats_perloc_mac1 df ########

#reorder Fst stats df
stats_perloc_mac1 <- cbind(stats_perloc_mac1, locnames_mac1)
stats_perloc_mac1 <- stats_perloc_mac1[order(stats_perloc_mac1$Contig_bp), ]
stats_perloc_mac1$NUM <- c(1:5718)

#pull outlier loci data
setdiff(outlierloc_mac1$contig_bp, stats_perloc_mac1$Contig_bp) #should = 0, as all should be present in stats_perloc
outlier_stats_mac1 <- stats_perloc_mac1[stats_perloc_mac1$Contig_bp %in% outlierloc_mac1$contig_bp, ] #df with only outliers
nonoutlier_stats_mac1 <- stats_perloc_mac1[!(stats_perloc_mac1$Contig_bp %in% outlierloc_mac1$contig_bp), ] #df with all other loci

#set outlier status
outlier_stats_mac1$status <- c(rep("Outlier", times = 80))
nonoutlier_stats_mac1$status <- c(rep("Not_Outlier", times = 5638))

#merge together
stats_perloc_mac1 <- rbind(outlier_stats_mac1, nonoutlier_stats_mac1)
stats_perloc_mac1 <-stats_perloc_mac1[order(stats_perloc_mac1$NUM), ]

#write out
write.csv(stats_perloc_mac1, "Data/Fst_perloc_mac1.csv")

#################################################################################################################################################

######## Visualize data ########
#designed to be run separately from earlier sections

remove(list = ls())

#read in data
stats_perloc <- read.csv("Data/Fst_perloc.csv", header = TRUE, row.names = 1)
stats_perloc_mac1 <- read.csv("Data/Fst_perloc_mac1.csv", header = TRUE, row.names = 1)

#add rows with loci for all boxplot
stats_perloc_all <- subset(stats_perloc, select = -c(status)) #duplicate df without status column
stats_perloc_all$status <- c(rep("All", times = 4212)) #add status back in (set to "all")

stats_perloc_all_mac1 <- subset(stats_perloc_mac1, select = -c(status)) #duplicate df without status column
stats_perloc_all_mac1$status <- c(rep("All", times = 5718)) #add status back in (set to "all")

#merge together
stats_perloc_boxplot <- rbind(stats_perloc, stats_perloc_all) #special df bc only want for boxplots (not pseudo-Manhattan plots)
stats_perloc_boxplot_mac1 <- rbind(stats_perloc_mac1, stats_perloc_all_mac1) #same but for mac1

######## mac > 2 plots ########

######## Scatter plot ########

#plot with outliers highlighted
fst_plot <- ggplot(data = stats_perloc, aes(x = NUM, y = Fst, color = status)) + 
  geom_point() + geom_hline(yintercept = 0.15, color = "black", size = 1, linetype = "dashed") + 
  scale_color_manual(values = c("gray74", "dodgerblue4"), labels = c("Not Outlier", "Outlier"))
fst_plot_annotated <- fst_plot + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 14, face = "bold"), legend.position = "top", 
        legend.text = element_text(size = 12), legend.title = element_text(size = 12))
fst_plot_annotated

######## Boxplots ########

#ordering x-axis
stats_perloc_boxplot$status2 <- factor(stats_perloc_boxplot$status, levels = c("Outlier", "Not_Outlier", "All")) #ordering X-axis

#boxplot
Fst_boxplot <- ggplot(data = stats_perloc_boxplot, aes(x = status2, y = Fst)) + geom_boxplot()
Fst_boxplot_annotated <- Fst_boxplot + xlab("Status") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(size = 1), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 14, color = "black"), 
        axis.title = element_text(size = 14, face = "bold"), legend.position = "top", 
        legend.text = element_text(size = 12), legend.title = element_text(size = 12))
Fst_boxplot_annotated

######## mac > 1 plots ########

######## Scatter plot ########

#plot with outliers highlighted
fst_mac1_plot <- ggplot(data = stats_perloc_mac1, aes(x = NUM, y = Fst, color = status)) + 
  geom_point() + geom_hline(yintercept = 0.15, color = "black", size = 1, linetype = "dashed") + 
  scale_color_manual(values = c("gray74", "dodgerblue4"), labels = c("Not Outlier", "Outlier"))
fst_mac1_plot_annotated <- fst_mac1_plot + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 14, face = "bold"), legend.position = "top", 
        legend.text = element_text(size = 12), legend.title = element_text(size = 12))
fst_mac1_plot_annotated

######## Boxplots ########

#ordering x-axis
stats_perloc_boxplot_mac1$status2 <- factor(stats_perloc_boxplot_mac1$status, levels = c("Outlier", "Not_Outlier", "All")) #ordering X-axis

#boxplot
Fst_mac1_boxplot <- ggplot(data = stats_perloc_boxplot_mac1, aes(x = status2, y = Fst)) + geom_boxplot()
Fst_mac1_boxplot_annotated <- Fst_mac1_boxplot + xlab("Status") + theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(size = 1), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 14, color = "black"), 
        axis.title = element_text(size = 14, face = "bold"), legend.position = "top", 
        legend.text = element_text(size = 12), legend.title = element_text(size = 12))
Fst_mac1_boxplot_annotated


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

allsites_pi$CHROM_BP <- paste(allsites_pi$CHROM, "_", allsites_pi$POS)

CHROM_BP <- allsites_pi$CHROM_BP
All_pi <- allsites_pi$PI
J_pi <- Jsites_pi$PI
N_pi <- Nsites_pi$PI
P_pi <- Psites_pi$PI
pi_bypop <- data.frame(cbind(CHROM_BP, All_pi, J_pi, P_pi, N_pi))
Contig_bp <- as.character(pi_bypop$Contig_bp)
J_pi <- as.numeric(as.character(pi_bypop$J))
N_pi <- as.numeric(as.character(pi_bypop$N))
P_pi <- as.numeric(as.character(pi_bypop$P))
All_pi <- as.numeric(as.character(pi_bypop$All))

colnames(pi_bypop) <- c("Contig_bp", "All", "J", "P", "N")
pi_bypop <- data.frame(pi_bypop[order(pi_bypop$Contig_bp),])

Fst_data <- read.csv("Data/Fst_perloc.csv", header = TRUE)
Fst <- Fst_data$Fst
Outlier_status <- as.character(Fst_data$status)
pi_fst <- data.frame(cbind())
fst_tot <- data.frame(cbind(Contig_bp, Outlier_status, All_pi, J_pi, P_pi, N_pi, Fst))

Allpi_fst <- ggplot(data = fst_tot, aes(x = All_pi, y = Fst, color = Outlier_status)) + 
  geom_point(size = 4) + theme_bw() + ggtitle("All pi") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = c(1, 0), axis.line = element_line(size = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_text(size = 26), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26, face = "bold")) + guides(color = guide_legend(override.aes = list(size = 8)))
Allpi_fst

Jpi_fst <- ggplot(data = fst_tot, aes(x = J_pi, y = Fst, color = Outlier_status)) + 
  geom_point(size = 4) + theme_bw() + ggtitle("J pi") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = c(1, 0), axis.line = element_line(size = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_text(size = 26), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26, face = "bold")) + guides(color = guide_legend(override.aes = list(size = 8)))
Jpi_fst

Ppi_fst <- ggplot(data = fst_tot, aes(x = P_pi, y = Fst, color = Outlier_status)) + 
  geom_point(size = 4) + theme_bw() + ggtitle("P pi") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = c(1, 0), axis.line = element_line(size = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_text(size = 26), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26, face = "bold")) + guides(color = guide_legend(override.aes = list(size = 8)))
Ppi_fst

Npi_fst <- ggplot(data = fst_tot, aes(x = N_pi, y = Fst, color = Outlier_status)) + 
  geom_point(size = 4) + theme_bw() + ggtitle("N pi") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.justification = c(1, 0), axis.line = element_line(size = 1), plot.title = element_blank(), 
        legend.position = "top", legend.text = element_text(size = 26), legend.title = element_text(size = 26), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 26, color = "black"), 
        axis.title = element_text(size = 26, face = "bold")) + guides(color = guide_legend(override.aes = list(size = 8)))
Npi_fst
