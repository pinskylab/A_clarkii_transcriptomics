#################################### Script to write simulation pops for Het metrics ##########################################

#Script to write simulation populations (bootstrap across individuals) to get 95% CIs for heterozygosity & genetic diversity estimates
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

#read in data
output_hicov2_snps_only <- read_csv("output.hicov2.snps.only.csv", col_names = FALSE)

#################################################################################################################################################

######## Create simulated VCFs ########

#redo for mac2 vcf

#pull loci names (contig & bp for each SNP)
loci_df <- output_hicov2_snps_only[, 1:2] #pull list of contigs & SNP bp
names(loci_df) <- c("CHROM", "pos") #change headers
loci_df$contig-bp <- paste(loci_df$CHROM, loci_df$pos, sep = "-") #create row of contig & bp combined
contig_bp <- loci_df$contig_bp[2:5730] #pull contig-bp column into a vector

#pull out first 9 columns containing information on SNPs & genotyping
output_info <- output_hicov2_snps_only[, 1:9]

#subset vcf into 3 populations for bootstrapping w/in a given population
output_J <- output_hicov2_snps_only[, 10:17]
output_I <- output_hicov2_snps_only[, 18:24]
output_P <- output_hicov2_snps_only[, 25:34]

#bootstrap pops
output_sample_J <- as.data.frame(t(apply(output_J, MARGIN = 1, FUN = sample, replace = TRUE))) #sample with replacement across rows for each pop (randomly select N genotypes at each locus independently)
output_sample_I <- as.data.frame(t(apply(output_I, MARGIN = 1, FUN = sample, replace = TRUE)))
output_sample_P <- as.data.frame(t(apply(output_P, MARGIN = 1, FUN = sample, replace = TRUE)))

######## Tidy data for genetic diversity calculations ########

output_boot <- cbind(output_J, output_sample_I, output_sample_P) #combine all bootstrapped pops into one dataframe
output_boot <- output_boot[2:5730, ] #remove row with header information
colnames(output_boot) <- c("J1", "J2", "J3", "J4", "J5", "J6", "J7", "J8", "N1", "N2", "N3", "N4", "N5", "N6", "N7", 
                           "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10") #create header row

#write & apply function to clean genotype data
clean <- function(x) {
  substitute <- gsub(pattern = ":.*", "", x) #removes everything after : and replaces with nothing
}

output_boot_clean <- apply(output_boot, MARGIN = 1, FUN = clean) #apply cleaning function to remove extra genotype information
output_boot <- gsub("/", "", output_boot_clean) #remove / separator in each genotype
output_boot <- as.data.frame(output_boot) #turn matrix from apply into datframe
colnames(output_boot) <- contig_bp #create header row w/contig_bp notation

pop <- c("J", "J", "J", "J", "J", "J", "J", "J", "N", "N", "N", "N", "N", "N", "N", 
         "P", "P", "P", "P", "P", "P", "P", "P", "P", "P") #vector w/population information for each individual

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