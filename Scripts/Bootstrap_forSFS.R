#################################### Script to create bootstrapped vcfs for easySFS ##########################################

#Script to write bootstrapped populations (across loci) to read into easySFS and create bootstrapped SFS for fastsimcoal2 CI estimations
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

#read in data
snps_data <- read_csv("Data/output.hicov2.snps.only.mac1.nooutliersfinal.contigsonly.reordered.csv", 
                                                   col_names = TRUE) #mac1 so keep informative rare alleles, no outliers to be neutral

#################################################################################################################################################

######## Clean vcf for bootstrapping ########

#pull out first 9 columns containing information on SNPs & genotyping
output_info <- snps_data[, 1:9]

#transform data so sample across loci
t_snps_data <- as.data.frame(t(snps_data))

######## Bootstrap across loci ########
bootstrap_list <- replicate(100, as.data.frame(t(sample(t_snps_data, replace = TRUE))), simplify = FALSE) #bootstrap w/replacement across loci, transform back into original VCF format
list <- as.data.frame(t(sample(t_snps_data, replace = TRUE)))

#unlist matrices into separate df for each bootstrapped sample
for (i in 1:100) {
  df <- data.frame(bootstrap_list[[i]])
  assign(paste('all', i , sep = ''), df)
}

#col for unique SNP id --> otherwise won't take all SNPs (if same name appears twice)
chrom <- c(rep("Contig", 5662))
num <- c(1:5662)
df <- as.data.frame(cbind(chrom, num))
df$Contig <- gsub(" ", "", paste(df$chrom, df$num))

#write out each df (repeat 10x, for 1-100 df)
all1 <- cbind(df$Contig, all1)
all1 <- all1[,-2]
colnames(all1)[1] <- "#CHROM"
write.table(all1, file = "Data/bootstrapped_vcfs_forSFS/mac1_nooutliers_reordered_boot1.txt", 
            sep = "\t",row.names = FALSE, col.names = TRUE, quote = FALSE)

#after write out all 100, add in VCF header info (in txt editor) and save as ".vcf"