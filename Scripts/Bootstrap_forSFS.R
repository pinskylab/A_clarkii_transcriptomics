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

#write out each df
all91 <- cbind(df$Contig, all91)
all91 <- all91[,-2]
colnames(all91)[1] <- "#CHROM"
write.table(all91, file = "Data/bootstrapped_vcfs_forSFS/mac1_nooutliers_reordered_boot91.txt", 
            sep = "\t",row.names = FALSE, col.names = TRUE, quote = FALSE)
all92 <- cbind(df$Contig, all92)
all92 <- all92[,-2]
colnames(all92)[1] <- "#CHROM"
write.table(all92, file = "Data/bootstrapped_vcfs_forSFS/mac1_nooutliers_reordered_boot92.txt", 
            sep = "\t",row.names = FALSE, col.names = TRUE, quote = FALSE)
all93 <- cbind(df$Contig, all93)
all93 <- all93[,-2]
colnames(all93)[1] <- "#CHROM"
write.table(all93, file = "Data/bootstrapped_vcfs_forSFS/mac1_nooutliers_reordered_boot93.txt", 
            sep = "\t",row.names = FALSE, col.names = TRUE, quote = FALSE)
all94 <- cbind(df$Contig, all94)
all94 <- all94[,-2]
colnames(all94)[1] <- "#CHROM"
write.table(all94, file = "Data/bootstrapped_vcfs_forSFS/mac1_nooutliers_reordered_boot94.txt", 
            sep = "\t",row.names = FALSE, col.names = TRUE, quote = FALSE)
all95 <- cbind(df$Contig, all95)
all95 <- all95[,-2]
colnames(all95)[1] <- "#CHROM"
write.table(all95, file = "Data/bootstrapped_vcfs_forSFS/mac1_nooutliers_reordered_boot95.txt", 
            sep = "\t",row.names = FALSE, col.names = TRUE, quote = FALSE)
all96 <- cbind(df$Contig, all96)
all96 <- all96[,-2]
colnames(all96)[1] <- "#CHROM"
write.table(all96, file = "Data/bootstrapped_vcfs_forSFS/mac1_nooutliers_reordered_boot96.txt", 
            sep = "\t",row.names = FALSE, col.names = TRUE, quote = FALSE)
all97 <- cbind(df$Contig, all97)
all97 <- all97[,-2]
colnames(all97)[1] <- "#CHROM"
write.table(all97, file = "Data/bootstrapped_vcfs_forSFS/mac1_nooutliers_reordered_boot97.txt", 
            sep = "\t",row.names = FALSE, col.names = TRUE, quote = FALSE)
all98 <- cbind(df$Contig, all98)
all98 <- all98[,-2]
colnames(all98)[1] <- "#CHROM"
write.table(all98, file = "Data/bootstrapped_vcfs_forSFS/mac1_nooutliers_reordered_boot98.txt", 
            sep = "\t",row.names = FALSE, col.names = TRUE, quote = FALSE)
all99 <- cbind(df$Contig, all99)
all99 <- all99[,-2]
colnames(all99)[1] <- "#CHROM"
write.table(all99, file = "Data/bootstrapped_vcfs_forSFS/mac1_nooutliers_reordered_boot99.txt", 
            sep = "\t",row.names = FALSE, col.names = TRUE, quote = FALSE)
all100 <- cbind(df$Contig, all100)
all100 <- all100[,-2]
colnames(all100)[1] <- "#CHROM"
write.table(all100, file = "Data/bootstrapped_vcfs_forSFS/mac1_nooutliers_reordered_boot100.txt", 
            sep = "\t",row.names = FALSE, col.names = TRUE, quote = FALSE)