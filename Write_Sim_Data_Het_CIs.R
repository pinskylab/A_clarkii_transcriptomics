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
library(purrr)

#read in data
output_hicov2_snps_only <- read_csv("output.hicov2.snps.only.csv", col_names = FALSE)

#################################################################################################################################################

######## Clean vcf for bootstrapping ########

#redo for mac2 vcf

#pull loci names (contig & bp for each SNP)
loci_df <- output_hicov2_snps_only[, 1:2] #pull list of contigs & SNP bp
names(loci_df) <- c("CHROM", "pos") #change headers
loci_df$contig_bp <- paste(loci_df$CHROM, loci_df$pos, sep = "_") #create row of contig & bp combined
contig_bp <- loci_df$contig_bp[2:5730] #pull contig-bp column into a vector

#pull out first 9 columns containing information on SNPs & genotyping
output_info <- output_hicov2_snps_only[, 1:9]

#subset vcf into 3 populations for bootstrapping w/in a given population
output_J <- output_hicov2_snps_only[, 10:17]
output_I <- output_hicov2_snps_only[, 18:24]
output_P <- output_hicov2_snps_only[, 25:34]

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
tot_J <- tot_J[2:5730, ] #remove header row
colnames(tot_J) <- paste("J", 1:80, sep = '') #create header row with unique individual ID for each sample

tot_I <- do.call("cbind", mget(ls(pattern = "^I"))) #search in env for any object starting with I and combine by columns
tot_I <- tot_I[2:5730, ]
colnames(tot_I) <- paste("I", 1:70, sep = '')

tot_P <- do.call("cbind", mget(ls(pattern = "^P"))) #search in env for any object starting with P and combine by columns
tot_P <- tot_P[2:5730, ]
colnames(tot_P) <- paste("I", 1:100, sep = '')

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
         rep("J100", 8) )

#write out dataframe as table
write.table(output_boot, file = "bootstrap_het.txt")

csv <- read.table(file = "bootstrap_het.txt", colClasses = "character", header = TRUE, row.names = 1) #colClasses as character to keep leading zeros

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
         rep("J100", 8) )


#turn dataframe into genind object to calculate het stats with
obj <- df2genind(output_boot, ploidy = 2, ncode = 1, pop = pop)

#alternatively write out as csv then read in as dataframe and run in diversity script

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