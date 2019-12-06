#################################### Script to write simulation pops for Het metrics ##########################################

#Script to write simulation populations (bootstrap across individuals) to get 95% CIs for heterozygosity & genetic diversity estimates
#Run 10x for 1000 bootstrapped samples for each population (unless run on HPC then increase # repeats in resampling step)
#Written for transcriptome project

#################################################################################################################################################

######## Set-up ########

#set working directory
setwd("C:/Users/rclar/Dropbox/Pinsky_Lab/Transcriptome_Proj/R_scripts/A_clarkii_transcriptomics/")
getwd()

remove(list = ls())

#load libraries
library(readr)
library(tidyverse)
library(adegenet)
library(hierfstat)
library(purrr)

#read in data
output_hicov2_snps_only_mac2 <- read_csv("../../VCFs_and_PLINK/output.hicov2.snps.only.mac2.csv", col_names = FALSE)
#quantiles <- read.csv("diversity_bootstrapped_cis.csv", header = TRUE, row.names = 1) #only bc if run all the way through hit errors where Rstudio can't support that many files open at once so won't create plots
#diversity_boot <- read.csv("diversity_bootstrapped.csv", header = TRUE, row.names = 1)

#################################################################################################################################################

######## Clean vcf for bootstrapping ########

#pull loci names (contig & bp for each SNP)
loci_df <- output_hicov2_snps_only_mac2[, 1:2] #pull list of contigs & SNP bp
names(loci_df) <- c("CHROM", "pos") #change headers
loci_df$contig_bp <- paste(loci_df$CHROM, loci_df$pos, sep = "_") #create row of contig & bp combined
contig_bp <- loci_df$contig_bp[2:4213] #pull contig-bp column into a vector

#pull out first 9 columns containing information on SNPs & genotyping
output_info <- output_hicov2_snps_only_mac2[, 1:9]

#subset vcf into 3 populations for bootstrapping w/in a given population
output_J <- output_hicov2_snps_only_mac2[, 10:17]
output_I <- output_hicov2_snps_only_mac2[, 18:24]
output_P <- output_hicov2_snps_only_mac2[, 25:34]

######## Bootstrap pops ########
sample_list_J <- replicate(100, as.data.frame(t(apply(output_J, MARGIN = 1, FUN = sample, replace = TRUE))), simplify = FALSE) #creates list of 1000 matrices that were created by sampling with replacement across rows for each pop (randomly select N genotypes at each locus independently)
sample_list_I <- replicate(100, as.data.frame(t(apply(output_I, MARGIN = 1, FUN = sample, replace = TRUE))), simplify = FALSE)
sample_list_P <- replicate(100, as.data.frame(t(apply(output_P, MARGIN = 1, FUN = sample, replace = TRUE))), simplify = FALSE)

######## Combine bootstrapped samples into one large dataframe for analysis ########

#unlist matrices into separate df for each boostrapped sample
for (i in 1:100) {
  df <- data.frame(sample_list_J[[i]])
  assign(paste('J', i , sep = ''), df)
}

for (i in 1:100) {
  df <- data.frame(sample_list_I[[i]])
  assign(paste('I', i , sep = ''), df)
}

for (i in 1:100) {
  df <- data.frame(sample_list_P[[i]])
  assign(paste('P', i , sep = ''), df)
}

#combine all dataframes into one large dataframe per pop
tot_J <- do.call("cbind", mget(ls(pattern = "^J"))) #search in env for any object starting with J and combine by columns
tot_J <- tot_J[2:4213, ] #remove header row
colnames(tot_J) <- paste("J", 1:800, sep = '') #create header row with unique individual ID for each sample

tot_I <- do.call("cbind", mget(ls(pattern = "^I"))) #search in env for any object starting with I and combine by columns
tot_I <- tot_I[2:4213, ]
colnames(tot_I) <- paste("I", 1:700, sep = '')

tot_P <- do.call("cbind", mget(ls(pattern = "^P"))) #search in env for any object starting with P and combine by columns
tot_P <- tot_P[2:4213, ]
colnames(tot_P) <- paste("I", 1:1000, sep = '')

#combine all pops into one large dataframe with all simulations for all populations
output_boot <- cbind(tot_J, tot_I, tot_P)

######## Tidy data for genetic diversity calculations ########

#write & apply function to clean genotype data
clean <- function(x) {
  substitute <- gsub(pattern = ":.*", "", x) #removes everything after : and replaces with nothing
}

output_boot_clean <- apply(output_boot, MARGIN = 1, FUN = clean) #apply cleaning function to remove extra genotype information
output_boot <- gsub("/", "", output_boot_clean) #remove / separator in each genotype
output_boot <- as.data.frame(output_boot) #turn matrix from apply into datframe
colnames(output_boot) <- contig_bp #create header row w/contig_bp notation

#write out dataframe as table
write.table(output_boot, file = "bootstrap_het.txt") #rename with subscript numbers so don't overwrite as run 10x

#################################################################################################################################################

######## Calculate Ho & He for each boostrapped sample ########

#read in data from bootstrapping runs earlier
bs_run_1 <- read.table(file = "bootstrap_het1.txt", colClasses = "character", header = TRUE, row.names = 1) #colClasses as character to keep leading zeros
bs_run_2 <- read.table(file = "bootstrap_het2.txt", colClasses = "character", header = TRUE, row.names = 1)
bs_run_3 <- read.table(file = "bootstrap_het3.txt", colClasses = "character", header = TRUE, row.names = 1)
bs_run_4 <- read.table(file = "bootstrap_het4.txt", colClasses = "character", header = TRUE, row.names = 1)
bs_run_5 <- read.table(file = "bootstrap_het5.txt", colClasses = "character", header = TRUE, row.names = 1)
bs_run_6 <- read.table(file = "bootstrap_het6.txt", colClasses = "character", header = TRUE, row.names = 1)
bs_run_7 <- read.table(file = "bootstrap_het7.txt", colClasses = "character", header = TRUE, row.names = 1)
bs_run_8 <- read.table(file = "bootstrap_het8.txt", colClasses = "character", header = TRUE, row.names = 1)
bs_run_9 <- read.table(file = "bootstrap_het9.txt", colClasses = "character", header = TRUE, row.names = 1)
bs_run_10 <- read.table(file = "bootstrap_het10.txt", colClasses = "character", header = TRUE, row.names = 1)

#set population vector to assign individuals to discrete bootstrap sample
pop <- c(rep("J1", 8), rep("J2", 8), rep("J3", 8), rep("J4", 8), rep("J5", 8), rep("J6", 8), rep("J7", 8), rep("J8", 8), 
         rep("J9", 8), rep("J10", 8), rep("J11", 8), rep("J12", 8), rep("J13", 8), rep("J14", 8), rep("J15", 8), 
         rep("J16", 8), rep("J17", 8), rep("J18", 8), rep("J19", 8), rep("J20", 8), rep("J21", 8), rep("J22", 8), 
         rep("J23", 8), rep("J24", 8), rep("J25", 8), rep("J26", 8), rep("J27", 8), rep("J28", 8), rep("J29", 8), 
         rep("J30", 8), rep("J31", 8), rep("J32", 8), rep("J33", 8), rep("J34", 8), rep("J35", 8), rep("J36", 8), 
         rep("J37", 8), rep("J38", 8), rep("J39", 8), rep("J40", 8), rep("J41", 8), rep("J42", 8), rep("J43", 8), 
         rep("J44", 8), rep("J45", 8), rep("J46", 8), rep("J47", 8), rep("J48", 8), rep("J49", 8), rep("50", 8), 
         rep("J51", 8), rep("J52", 8), rep("J53", 8), rep("J54", 8), rep("J55", 8), rep("J56", 8), rep("J57", 8), 
         rep("J58", 8), rep("J59", 8), rep("J60", 8), rep("J61", 8), rep("J62", 8), rep("J63", 8), rep("J64", 8), 
         rep("J65", 8), rep("J66", 8), rep("J67", 8), rep("J68", 8), rep("J69", 8), rep("J70", 8), rep("J71", 8), 
         rep("J72", 8), rep("J73", 8), rep("J74", 8), rep("J75", 8), rep("J76", 8), rep("J77", 8), rep("J78", 8), 
         rep("J79", 8), rep("J80", 8), rep("J81", 8), rep("J82", 8), rep("J83", 8), rep("J84", 8), rep("J85", 8), 
         rep("J86", 8), rep("J87", 8), rep("J88", 8), rep("J89", 8), rep("J90", 8), rep("J91", 8), rep("J92", 8), 
         rep("J93", 8), rep("J94", 8), rep("J95", 8), rep("J96", 8), rep("J97", 8), rep("J98", 8), rep("J99", 8), 
         rep("J100", 8), rep("I1", 7), rep("I2", 7), rep("I3", 7), rep("I4", 7), rep("I5", 7), rep("II6", 7), 
         rep("I7", 7), rep("I8", 7), rep("I9", 7), rep("I10", 7), rep("I11", 7), rep("I12", 7), rep("I13", 7), 
         rep("I14", 7), rep("I15", 7), rep("I16", 7), rep("I17", 7), rep("I18", 7), rep("I19", 7), rep("I20", 7), 
         rep("I21", 7), rep("I22", 7), rep("I23", 7), rep("I24", 7), rep("I25", 7), rep("I26", 7), rep("I27", 7), 
         rep("I28", 7), rep("I29", 7), rep("I30", 7), rep("I31", 7), rep("I32", 7), rep("I33", 7), rep("I34", 7), 
         rep("I35", 7), rep("I36", 7), rep("I37", 7), rep("I38", 7), rep("I39", 7), rep("I40", 7), rep("I41", 7), 
         rep("I42", 7), rep("I43", 7), rep("I44", 7), rep("I45", 7), rep("I46", 7), rep("I47", 7), rep("I48", 7), 
         rep("I49", 7), rep("I50", 7), rep("I51", 7), rep("I52", 7), rep("I53", 7), rep("I54", 7), rep("I55", 7), 
         rep("I56", 7), rep("I57", 7), rep("I58", 7), rep("I59", 7), rep("I60", 7), rep("I61", 7), rep("I62", 7), 
         rep("I63", 7), rep("I64", 7), rep("I65", 7), rep("I66", 7), rep("I67", 7), rep("I68", 7), rep("I69", 7), 
         rep("I70", 7), rep("I71", 7), rep("I72", 7), rep("I73", 7), rep("I74", 7), rep("I75", 7), rep("I76", 7), 
         rep("I77", 7), rep("I78", 7), rep("I79", 7), rep("I80", 7), rep("I81", 7), rep("I82", 7), rep("I83", 7), 
         rep("I84", 7), rep("I85", 7), rep("I86", 7), rep("I87", 7), rep("I88", 7), rep("I89", 7), rep("I90", 7), 
         rep("I91", 7), rep("I92", 7), rep("I93", 7), rep("I94", 7), rep("I95", 7), rep("I96", 7), rep("I97", 7), 
         rep("I98", 7), rep("I99", 7), rep("I100", 7),rep("P1", 10), rep("P2", 10), rep("P3", 10), rep("P4", 10), 
         rep("P5", 10), rep("PI6", 10), rep("P7", 10), rep("P8", 10), rep("P9", 10), rep("P10", 10), rep("P11", 10), 
         rep("P12", 10), rep("P13", 10), rep("P14", 10), rep("P15", 10), rep("P16", 10), rep("P17", 10), rep("P18", 10), 
         rep("P19", 10), rep("P20", 10), rep("P21", 10), rep("P22", 10), rep("P23", 10), rep("P24", 10), rep("P25", 10), 
         rep("P26", 10), rep("P27", 10), rep("P28", 10), rep("P29", 10), rep("P30", 10), rep("P31", 10), rep("P32", 10), 
         rep("P33", 10), rep("P34", 10), rep("P35", 10), rep("P36", 10), rep("P37", 10), rep("P38", 10), rep("P39", 10), 
         rep("P40", 10), rep("P41", 10), rep("P42", 10), rep("P43", 10), rep("P44", 10), rep("P45", 10), rep("P46", 10), 
         rep("P47", 10), rep("P48", 10), rep("P49", 10), rep("P50", 10), rep("P51", 10), rep("P52", 10), rep("P53", 10), 
         rep("P54", 10), rep("P55", 10), rep("P56", 10), rep("P57", 10), rep("P58", 10), rep("P59", 10), rep("P60", 10), 
         rep("P61", 10), rep("P62", 10), rep("P63", 10), rep("P64", 10), rep("P65", 10), rep("P66", 10), rep("P67", 10), 
         rep("P68", 10), rep("P69", 10), rep("P70", 10), rep("P71", 10), rep("P72", 10), rep("P73", 10), rep("P74", 10), 
         rep("P75", 10), rep("P76", 10), rep("P77", 10), rep("P78", 10), rep("P79", 10), rep("P80", 10), rep("P81", 10), 
         rep("P82", 10), rep("P83", 10), rep("P84", 10), rep("P85", 10), rep("P86", 10), rep("P87", 10), rep("P88", 10), 
         rep("P89", 10), rep("P90", 10), rep("P91", 10), rep("P92", 10), rep("P93", 10), rep("P94", 10), rep("P95", 10), 
         rep("P96", 10), rep("P97", 10), rep("P98", 10), rep("P99", 10), rep("P100", 10))

#turn dataframe into genind object to calculate het stats with
bs_genind_1 <- df2genind(bs_run_1, ploidy = 2, ncode = 1, pop = pop)
bs_genind_2 <- df2genind(bs_run_2, ploidy = 2, ncode = 1, pop = pop)
bs_genind_3 <- df2genind(bs_run_3, ploidy = 2, ncode = 1, pop = pop)
bs_genind_4 <- df2genind(bs_run_4, ploidy = 2, ncode = 1, pop = pop)
bs_genind_5 <- df2genind(bs_run_5, ploidy = 2, ncode = 1, pop = pop)
bs_genind_6 <- df2genind(bs_run_6, ploidy = 2, ncode = 1, pop = pop)
bs_genind_7 <- df2genind(bs_run_7, ploidy = 2, ncode = 1, pop = pop)
bs_genind_8 <- df2genind(bs_run_8, ploidy = 2, ncode = 1, pop = pop)
bs_genind_9 <- df2genind(bs_run_9, ploidy = 2, ncode = 1, pop = pop)
bs_genind_10 <- df2genind(bs_run_10, ploidy = 2, ncode = 1, pop = pop)

#calculate het stats (per locus)
sum_stats_1 <- basic.stats(bs_genind_1)
sum_stats_2 <- basic.stats(bs_genind_2)
sum_stats_3 <- basic.stats(bs_genind_3)
sum_stats_4 <- basic.stats(bs_genind_4)
sum_stats_5 <- basic.stats(bs_genind_5)
sum_stats_6 <- basic.stats(bs_genind_6)
sum_stats_7 <- basic.stats(bs_genind_7)
sum_stats_8 <- basic.stats(bs_genind_8)
sum_stats_9 <- basic.stats(bs_genind_9)
sum_stats_10 <- basic.stats(bs_genind_10)

#pull out Ho for each run
Ho_J_1 <- colMeans(sum_stats_1$Ho[,1:100])
Ho_I_1 <- colMeans(sum_stats_1$Ho[,101:200])
Ho_P_1 <- colMeans(sum_stats_1$Ho[,201:300])
Means_Ho_1 <- data.frame(Ho_J_1, Ho_I_1, Ho_P_1)
rownames(Means_Ho_1) <- paste("S", 1:100, sep = '')
colnames(Means_Ho_1) <- c("J", "I", "P")

Ho_J_2 <- colMeans(sum_stats_2$Ho[,1:100])
Ho_I_2 <- colMeans(sum_stats_2$Ho[,101:200])
Ho_P_2 <- colMeans(sum_stats_2$Ho[,201:300])
Means_Ho_2 <- data.frame(Ho_J_2, Ho_I_2, Ho_P_2)
rownames(Means_Ho_2) <- paste("S", 101:200, sep = '')
colnames(Means_Ho_2) <- c("J", "I", "P")

Ho_J_3 <- colMeans(sum_stats_3$Ho[,1:100])
Ho_I_3 <- colMeans(sum_stats_3$Ho[,101:200])
Ho_P_3 <- colMeans(sum_stats_3$Ho[,201:300])
Means_Ho_3 <- data.frame(Ho_J_3, Ho_I_3, Ho_P_3)
rownames(Means_Ho_3) <- paste("S", 201:300, sep = '')
colnames(Means_Ho_3) <- c("J", "I", "P")

Ho_J_4 <- colMeans(sum_stats_4$Ho[,1:100])
Ho_I_4 <- colMeans(sum_stats_4$Ho[,101:200])
Ho_P_4 <- colMeans(sum_stats_4$Ho[,201:300])
Means_Ho_4 <- data.frame(Ho_J_4, Ho_I_4, Ho_P_4)
rownames(Means_Ho_4) <- paste("S", 301:400, sep = '')
colnames(Means_Ho_4) <- c("J", "I", "P")

Ho_J_5 <- colMeans(sum_stats_5$Ho[,1:100])
Ho_I_5 <- colMeans(sum_stats_5$Ho[,101:200])
Ho_P_5 <- colMeans(sum_stats_5$Ho[,201:300])
Means_Ho_5 <- data.frame(Ho_J_5, Ho_I_5, Ho_P_5)
rownames(Means_Ho_5) <- paste("S", 401:500, sep = '')
colnames(Means_Ho_5) <- c("J", "I", "P")

Ho_J_6 <- colMeans(sum_stats_6$Ho[,1:100])
Ho_I_6 <- colMeans(sum_stats_6$Ho[,101:200])
Ho_P_6 <- colMeans(sum_stats_6$Ho[,201:300])
Means_Ho_6 <- data.frame(Ho_J_6, Ho_I_6, Ho_P_6)
rownames(Means_Ho_6) <- paste("S", 501:600, sep = '')
colnames(Means_Ho_6) <- c("J", "I", "P")

Ho_J_7 <- colMeans(sum_stats_7$Ho[,1:100])
Ho_I_7 <- colMeans(sum_stats_7$Ho[,101:200])
Ho_P_7 <- colMeans(sum_stats_7$Ho[,201:300])
Means_Ho_7 <- data.frame(Ho_J_7, Ho_I_7, Ho_P_7)
rownames(Means_Ho_7) <- paste("S", 601:700, sep = '')
colnames(Means_Ho_7) <- c("J", "I", "P")

Ho_J_8 <- colMeans(sum_stats_8$Ho[,1:100])
Ho_I_8 <- colMeans(sum_stats_8$Ho[,101:200])
Ho_P_8 <- colMeans(sum_stats_8$Ho[,201:300])
Means_Ho_8 <- data.frame(Ho_J_8, Ho_I_8, Ho_P_8)
rownames(Means_Ho_8) <- paste("S", 701:800, sep = '')
colnames(Means_Ho_8) <- c("J", "I", "P")

Ho_J_9 <- colMeans(sum_stats_9$Ho[,1:100])
Ho_I_9 <- colMeans(sum_stats_9$Ho[,101:200])
Ho_P_9 <- colMeans(sum_stats_9$Ho[,201:300])
Means_Ho_9 <- data.frame(Ho_J_9, Ho_I_9, Ho_P_9)
rownames(Means_Ho_9) <- paste("S", 801:900, sep = '')
colnames(Means_Ho_9) <- c("J", "I", "P")

Ho_J_10 <- colMeans(sum_stats_10$Ho[,1:100])
Ho_I_10 <- colMeans(sum_stats_10$Ho[,101:200])
Ho_P_10 <- colMeans(sum_stats_10$Ho[,201:300])
Means_Ho_10 <- data.frame(Ho_J_10, Ho_I_10, Ho_P_10)
rownames(Means_Ho_10) <- paste("S", 901:1000, sep = '')
colnames(Means_Ho_10) <- c("J", "I", "P")

#merge all Ho dataframes together
Ho_means_all <- rbind(Means_Ho_1, Means_Ho_2, Means_Ho_3, Means_Ho_4, Means_Ho_5, Means_Ho_6, Means_Ho_7, Means_Ho_8,
                      Means_Ho_9, Means_Ho_10)
colnames(Ho_means_all) <- c("J_Ho", "I_Ho", "P_Ho")
#write.csv(Ho_means_all, "Ho_bootstrapped_all.csv") #to write out along the way

#pull out He for each run
He_J_1 <- colMeans(sum_stats_1$Hs[,1:100])
He_I_1 <- colMeans(sum_stats_1$Hs[,101:200])
He_P_1 <- colMeans(sum_stats_1$Hs[,201:300])
Means_He_1 <- data.frame(He_J_1, He_I_1, He_P_1)
rownames(Means_He_1) <- paste("S", 1:100, sep = '')
colnames(Means_He_1) <- c("J", "I", "P")

He_J_2 <- colMeans(sum_stats_2$Hs[,1:100])
He_I_2 <- colMeans(sum_stats_2$Hs[,101:200])
He_P_2 <- colMeans(sum_stats_2$Hs[,201:300])
Means_He_2 <- data.frame(He_J_2, He_I_2, He_P_2)
rownames(Means_He_2) <- paste("S", 101:200, sep = '')
colnames(Means_He_2) <- c("J", "I", "P")

He_J_3 <- colMeans(sum_stats_3$Hs[,1:100])
He_I_3 <- colMeans(sum_stats_3$Hs[,101:200])
He_P_3 <- colMeans(sum_stats_3$Hs[,201:300])
Means_He_3 <- data.frame(He_J_3, He_I_3, He_P_3)
rownames(Means_He_3) <- paste("S", 201:300, sep = '')
colnames(Means_He_3) <- c("J", "I", "P")

He_J_4 <- colMeans(sum_stats_4$Hs[,1:100])
He_I_4 <- colMeans(sum_stats_4$Hs[,101:200])
He_P_4 <- colMeans(sum_stats_4$Hs[,201:300])
Means_He_4 <- data.frame(He_J_4, He_I_4, He_P_4)
rownames(Means_He_4) <- paste("S", 301:400, sep = '')
colnames(Means_He_4) <- c("J", "I", "P")

He_J_5 <- colMeans(sum_stats_5$Hs[,1:100])
He_I_5 <- colMeans(sum_stats_5$Hs[,101:200])
He_P_5 <- colMeans(sum_stats_5$Hs[,201:300])
Means_He_5 <- data.frame(He_J_5, He_I_5, He_P_5)
rownames(Means_He_5) <- paste("S", 401:500, sep = '')
colnames(Means_He_5) <- c("J", "I", "P")

He_J_6 <- colMeans(sum_stats_6$Hs[,1:100])
He_I_6 <- colMeans(sum_stats_6$Hs[,101:200])
He_P_6 <- colMeans(sum_stats_6$Hs[,201:300])
Means_He_6 <- data.frame(He_J_6, He_I_6, He_P_6)
rownames(Means_He_6) <- paste("S", 501:600, sep = '')
colnames(Means_He_6) <- c("J", "I", "P")

He_J_7 <- colMeans(sum_stats_7$Hs[,1:100])
He_I_7 <- colMeans(sum_stats_7$Hs[,101:200])
He_P_7 <- colMeans(sum_stats_7$Hs[,201:300])
Means_He_7 <- data.frame(He_J_7, He_I_7, He_P_7)
rownames(Means_He_7) <- paste("S", 601:700, sep = '')
colnames(Means_He_7) <- c("J", "I", "P")

He_J_8 <- colMeans(sum_stats_8$Hs[,1:100])
He_I_8 <- colMeans(sum_stats_8$Hs[,101:200])
He_P_8 <- colMeans(sum_stats_8$Hs[,201:300])
Means_He_8 <- data.frame(He_J_8, He_I_8, He_P_8)
rownames(Means_He_8) <- paste("S", 701:800, sep = '')
colnames(Means_He_8) <- c("J", "I", "P")

He_J_9 <- colMeans(sum_stats_9$Hs[,1:100])
He_I_9 <- colMeans(sum_stats_9$Hs[,101:200])
He_P_9 <- colMeans(sum_stats_9$Hs[,201:300])
Means_He_9 <- data.frame(He_J_9, He_I_9, He_P_9)
rownames(Means_He_9) <- paste("S", 801:900, sep = '')
colnames(Means_He_9) <- c("J", "I", "P")

He_J_10 <- colMeans(sum_stats_10$Hs[,1:100])
He_I_10 <- colMeans(sum_stats_10$Hs[,101:200])
He_P_10 <- colMeans(sum_stats_10$Hs[,201:300])
Means_He_10 <- data.frame(He_J_10, He_I_10, He_P_10)
rownames(Means_He_10) <- paste("S", 901:1000, sep = '')
colnames(Means_He_10) <- c("J", "I", "P")

#merge all Ho dataframes together
He_means_all <- rbind(Means_He_1, Means_He_2, Means_He_3, Means_He_4, Means_He_5, Means_He_6, Means_He_7, Means_He_8,
                      Means_He_9, Means_He_10)
colnames(He_means_all) <- c("J_He", "I_He", "P_He")
#write.csv(He_means_all, "He_bootstrapped_all.csv") #to write out along the way

#combine He & Ho into one df
diversity_boot <- cbind(Ho_means_all, He_means_all)

######## Calculate Fis for bootstrapped samples########

#calculate Fis as (1- Ho/He)
diversity_boot$J_Fis <- 1 - (diversity_boot$J_Ho/diversity_boot$J_He)
diversity_boot$I_Fis <- 1 - (diversity_boot$I_Ho/diversity_boot$I_He)
diversity_boot$P_Fis <- 1 - (diversity_boot$P_Ho/diversity_boot$P_He)

######## Calculate quantiles for all metrics ########

#calculate 2.5 & 97.5 quantiles
quantiles <- data.frame(t(apply(diversity_boot, MARGIN = 2, FUN = quantile, c(0.025, 0.50, 0.975)))) #calculates quantiles for each column in diversity_boot and transform df so in tidy format
colnames(quantiles) <- c("2.5_per", "median", "97.5_per")
quantiles$pop <- c(rep(c("Japan", "Indonesia", "Philippines"), 3)) #add pop column for CI visualization
quantiles$pop <- factor(quantiles$pop, levels = c("Japan", "Philippines", "Indonesia")) #add levels for CI visualization
quantiles$metric <- c(rep("Ho", 3), rep("He", 3), rep("Fis", 3)) #add metric column for CI visualization
quantiles$diff_lower <- quantiles$median - quantiles$`2.5_per` #calculate diff btwn median and 2.5 percentile for CI visualization
quantiles$diff_upper <- quantiles$`97.5_per` - quantiles$median # calculate diff btwn median and 97.5 percentile for CI visualization

#write out data
write.csv(diversity_boot, "diversity_bootstrapped.csv")
write.csv(quantiles, "diversity_bootstrapped_cis.csv")

################################################################################################################################################

######## Visualize results ########

#plot of median Ho w/95% CI error bars
Ho_plot <- ggplot(data = quantiles[which(quantiles$metric == "Ho"), ], aes(x = pop, y = median)) + 
  geom_point(aes(size = 1), show.legend = FALSE) + 
  geom_errorbar(aes(ymin = median - diff_lower, ymax = median + diff_upper, width = 0.5, size = 0.5), show.legend = FALSE) + 
  ylim(0.19, 0.26) + ggtitle("Ho w/95% CI") + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 14, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 12, color = "black"))
Ho_plot

#plot of median He w/95% CI error bars
He_plot <- ggplot(data = quantiles[which(quantiles$metric == "He"), ], aes(x = pop, y = median)) + 
  geom_point(aes(size = 1), show.legend = FALSE) + 
  geom_errorbar(aes(ymin = median - diff_lower, ymax = median + diff_upper, width = 0.5, size = 0.5), show.legend = FALSE) + 
  ylim(0.19, 0.26) + ggtitle("He w/95% CI") + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 14, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 12, color = "black"))
He_plot

#plot of median Fis w/95% CI error bars
Fis_plot <- ggplot(data = quantiles[which(quantiles$metric == "Fis"), ], aes(x = pop, y = median)) + 
  geom_point(aes(size = 1), show.legend = FALSE) + 
  geom_errorbar(aes(ymin = median - diff_lower, ymax = median + diff_upper, width = 0.5, size = 0.5), show.legend = FALSE) + 
  ggtitle("Fis w/95% CI") + theme_bw() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 14, face = "bold"), 
        axis.ticks = element_line(color = "black", size = 1), axis.text = element_text(size = 12, color = "black"))
Fis_plot


######## Optional: turn genotypes into separate columns (1 per allele) ########
#ideally split this way, randomly sample from alleles and then merge back together two columns at a time
#dt <- data.frame(id = 1:2,
 #                AIN5997 = c("01/02", "01/02"),
  #               AIN7452 = c("02/02", NA),
   #              AIN8674 = c("02/02", "02/02"), stringsAsFactors = F)

#function for splitting genotypes into different columns
#f <- function(x) {
 # dt %>% select("id", x) %>% separate_(x, paste0(x, c(".1", ".2")))
#}

#dt_new <- names(dt)[2:4] %>%
 # map(f) %>%
#reduce(left_join, by = "id")