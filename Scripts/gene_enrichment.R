################################################## Script for Gene Enrichment  ######################################################

#Created for transcriptome project
#Uses snp2go package among others (snp2go pulled from snp2go.R function on GitHub -> davidszkiba/snp2go)
#Run on Amphiprion (correct version of GenomicRanges, etc.)

#################################################################################################################################################

######## Set-up ########

#set working directory
setwd("C:/Users/Rene/Dropbox/Pinsky_Lab/Transcriptome_Proj/R_scripts/A_clarkii_transcriptomics/")
getwd()

remove(list = ls())

#load libraries
library(hash)
library(goProfiles)
library(GenomicRanges)

#load functions
source("Scripts/snp2go.R") #pulled from davidszkiba/snp2go/snp2go.R on GitHub (as all other functions are)
source("Scripts/ProcessGFF.R")
source("Scripts/ProcessMartExport.R")
source("Scripts/GetProgressBarIncrementer.R")
source("Scripts/ComputeLevels.R")
source("Scripts/RemoveFromUpperLevels.R")
source("Scripts/ComputeSNPsByLevel.R")
source("Scripts/GetAllRegionsAndSNPs.R")
source("Scripts/GetCandidateTerms.R")
source("Scripts/GetStatistic.R")
source("Scripts/FindLevel.R")
source("Scripts/SamplingTest.R")
source("Scripts/AddOffspringInformation.R")

#read in data
cand_snps <- read.delim("../../VCFs_and_PLINK/output.hicov2.snps.only.outliersonly.strict.vcf", header = FALSE, comment.char = "#") #reading in VCF file with cand SNPs only, ignoring header comments
noncand_snps <- read.delim("../../VCFs_and_PLINK/output.hicov2.snps.only.mac2.nooutliers.strict.vcf", header = FALSE, comment.char = "#") #reading in VCF file with noncand SNPs only, ignoring header comments
#noncand_snps <- read.delim("../Desktop/Transcriptome_Proj/SNP_VCFs/output.hicov2.snps.only.nooutliers.strict.vcf", header = FALSE, comment.char = "#") #reading in VCF file with noncand SNPs only, ignoring header comments

################################################################################################################################################

######## Create GenomicRanges object from VCFs ########
cand_snps[, 2] <- as.numeric(cand_snps[, 2]) #making sure SNP bp column is numeric
cand_snps_GR <- GRanges(seqnames = cand_snps[, 1], ranges = IRanges(cand_snps[, 2], cand_snps[, 2])) #converting cand_snps to GenomicRanges object

noncand_snps[, 2] <- as.numeric(noncand_snps[, 2])
noncand_snps_GR <- GRanges(seqnames = noncand_snps[, 1], ranges = IRanges(noncand_snps[, 2], noncand_snps[, 2]))

######## Gene enrichment analysis ########

#run on Amphiprion bc have older version of bioconductor that is compatable with lapply()
snp2go_analysis <- snp2go(gff = "snp2go_annotation_MFBP.gff", candidateSNPs = cand_snps_GR, noncandidateSNPs = noncand_snps_GR)

gff.significant.terms <- snp2go_analysis$enriched$GO
gff.significant.regions <- snp2go_analysis$regions
enriched <- snp2go_analysis$enriched

#write out
write.table(file = "../../Annotation_Results/snp2go_enriched_gff.tsv", enriched, sep = "\t", row.names = FALSE)