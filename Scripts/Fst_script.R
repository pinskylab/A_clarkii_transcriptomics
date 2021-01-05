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
outlierloc <- read.csv("Data/outlier_loci_mac2.csv", header = TRUE) #to add outlier status

#################################################################################################################################################

######## Calculate Fst per locus --> mac>2 ########

#add pop levels to genind objects
Pop <- c(rep(1, times = 8), rep(2, times = 7), rep(3, times = 10)) #1 = Japan, 2 = Indonesia, 3 = Philippines
pop(mac2_genind) <- Pop
mac2_genind #check to make sure three pops 

#calculate Fst per locus
mac2_stats <- basic.stats(mac2_genind)
stats_perloc <- data.frame(mac2_stats$perloc)

######## Calculate pairwise-Fst ########

mac2_hierf <- genind2hierfstat(mac2_genind) #convert to hierfstat db for pairwise analyses
pairwise_fst <- genet.dist(mac2_hierf, method = "WC84") #calculates Weir & Cockerham's Fst

#bootstrap pairwise_fst for 95% CI

#need to convert pop character to numeric for bootstrap to work
mac2_hierf$pop <- as.numeric(mac2_hierf$pop)
class(mac2_hierf$pop) #check to make sure numeric

#bootstrap pairwise estimates
pairwise_boot <- boot.ppfst(dat = mac2_hierf, nboot = 1000, quant = c(0.025, 0.975), diploid = TRUE)

#get 95% CI limits
ci_upper <- pairwise_boot$ul
ci_lower <- pairwise_boot$ll

######## Add outlier status to stats_perloc df ########

#reorder Fst stats df
stats_perloc <- cbind(stats_perloc, locnames) #add locus names to identify by
stats_perloc <- stats_perloc[order(stats_perloc$Contig_bp), ] #order by contig_bp
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
stats_perloc <-stats_perloc[order(stats_perloc$NUM), ] #reorder

#write out
write.csv(stats_perloc, "Data/Fst_perloc.csv")

#################################################################################################################################################

######## Visualize data ########
#designed to be run separately from earlier sections

remove(list = ls())

#load libraries
library(tidyverse)

#read in data
stats_perloc <- read.csv("Data/Fst_perloc.csv", header = TRUE, row.names = 1)

#add rows with loci for all boxplot
stats_perloc_all <- subset(stats_perloc, select = -c(status)) #duplicate df without status column
stats_perloc_all$status <- c(rep("All", times = 4212)) #add status back in (set to "all")

#merge together
stats_perloc_boxplot <- rbind(stats_perloc, stats_perloc_all) #special df bc only want for boxplots (not pseudo-Manhattan scatterplots)

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