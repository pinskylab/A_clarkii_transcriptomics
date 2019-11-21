############################################## Calibrate Xtx significance cut-off for BayPass ##################################################

#follows Gautier (2015) instructions for calibrating and determining Xtx significance cut-offs
#creates and assesses quality of PODs (pseudo-observed datasets) and generates cut-off for Xtx based on dataset
#sources functions included in R script that comes with BayPass download

################################################################################################################################################

######## Set-up ########

remove(list = ls())

#Set working directory
setwd("C:/Users/Rene/Dropbox/Pinsky_Lab/Transcriptome_Proj/R_scripts/A_clarkii_transcriptomics/")
getwd()

#load libraries
library(corrplot)
library(mvtnorm)

#source functions
source("baypass_utils.R") #script that comes with BayPass download with functions used below (Gautier 2015)

#read in data
#read in output from BayPass runs with real dataset (4212 loci)
omega_realdata <- as.matrix(read.table("../../BayPass_Output/mac2_SNPs/clownfish_mac2_mat.cov")) #read in posterior distribution of covariance matrix
pi.beta.coef_realdata <- read.table("../../BayPass_Output/mac2_SNPs/aux1_mac2_summary_beta_params.out", h = T)$Mean #read in posterior distributions of mean a(pi) and mean b(pi)
XtX_realdata <- read.table("../../BayPass_Output/mac2_SNPs/aux1_mac2_summary_pi_xtx.out", h = T)$M_XtX #read in posterior distribution of Xtx values for each loci
geno <- geno2YN("../../BayPass_Output/mac2_SNPs/clownfish_mac2.geno") #read in allele counts in BayPass format to convert to total allele counts (per SNP)

#read in output from POD runs with pseudo-observed dataset
omega_POD <- as.matrix(read.table(".../../BayPass_Output/mac2_SNPs/core1_PODs/POD_4212_mat_omega.out"))
#pi.beta.coef_realdata <- read.table("../../BayPass_Output/mac2_SNPs/core1_PODs/core1_summary_beta_params.out", h = T)$Mean
pi.beta.coef_POD <- read.table("../../BayPass_Output/mac2_SNPs/core1_PODs/POD_4212_summary_beta_params.out", h = T)$Mean
XtX_POD <- read.table("../../BayPass_Output/mac2_SNPs/core1_PODs/POD_4212_summary_pi_xtx.out", h = T)$M_XtX

################################################################################################################################################

######## Create POD ########

#create POD
simulate.baypass(omega.mat = omega_realdata, nsnp = 4212, 
                 sample.size = geno$NN, beta.pi = pi.beta.coef_realdata, pi.maf = 0, 
                 suffix = "aclarkiipods") #simulate 4212 loci POD

######## Compare POD and real data BayPass output ########

#format omega_realdata dataframe
pop.names = c("J", "P", "N")
dimnames(omega_realdata) <- list(pop.names, pop.names)

#compute and visualize correlation matrix
cor.mat_realdata <- cov2cor(omega_realdata)
corrplot(cor.mat_realdata, method = "color", mar = c(2, 1, 2, 2) + 0.1, 
         main = expression("Correlation map based on" ~hat(Omega))) #heat map of allele freq correlation by population

#compare estimate of omega from POD & real data
omega_compar <- plot(omega_POD, omega_realdata) + abline(a = 0, b = 1) #scatterplot of omegas with 1:1 line
omega_compar

#compare estimate of a_pi & b_pi hyperparameters from POD & real data
a_b_pi_compar <- plot(pi.beta.coef_POD, pi.beta.coef_realdata) + abline(a = 0, b = 1) #scatter plot of hyperparameters with 1:1 line
a_b_pi_compar

######## Calculate XtX significance cut-off ########

#XtX calibration
POD_0.01_thresh <- quantile(XtX_POD, probs = 0.99) #calculate XtX 99% cut-off (value at which only 1% of POD posterior distribution of XtX is greater than)
XtX_plot <- plot(XtX_realdata) + abline(h = POD_0.01_thresh, lty = 2) #plot of XtXs from real dataset with cut-off line
XtX_plot
