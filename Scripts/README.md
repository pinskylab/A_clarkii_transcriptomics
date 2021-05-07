**Scripts:**
 * **Allele Freq Line Plots.R:** Creates polarized allele frequency plots (outlier and non-outlier).
 * **Pull_BFs.R:** Reads in raw "summary_betai.txt" files from BayPass and reorganizes into one .csv file for easier downstream analysis. Also pulls out candidate SNPs.
 * **XtX_Calibration.R:** Modified code from Gautier (2015) to create pseudo-observed datasets (PODs) to generate 95% significance thresholds for XtX values. Once PODs are run in BayPass, contains modified code to calculate XtX significance cut-off.
 * **BF Correlation Plots.R:** Reads in .csv file of Bayes Factors and XtX values outputted by BayPass and creates correlation plots by environmental variables.
 * **Write Simulation Pop Data for BayPass.R:** Creates permuted datasets to run in BayPass and create distributions of Bayes Factors under null hypothesis (no association between allele frequencies and environmental factors).
 * **ECDFS for Sim v Real BFs.R:** Contains code for Mann-Whitney U-tests and plots of empirical cumulative distribution functions for Bayes Factors from raw data and permuted data.
 * **RDA.R:** Runs RDA, creates bi-plots for visualization, and identifies outliers from RDA analysis using two different methods.
 * **Fst_script.R:** Calculates per-SNP Fst and site pairwise Fst values from VCF.
 * **PCAs.R:** Reads in eigenvector information from plink and creates PCA plots.
 * **STRUCTURE_script.R:** Reads in STRUCTURE output files, runs CLUMPP and creates output plots for visualization of STRUCTURE results. Also creates ML and Evanno method plots to identify the "best" value of K.
 * **Bootstrap_forSFS.R:** Creates bootstrapped VCF files. These files can then be read into easySFS to create bootstrapped SFS for reading into *fastsimcoal2* and calculating 95% CIs.
 * **fsc_CIs.R:** Reads in the parameter estimates the best maximum likelihood run for each bootstrapped SFS from *fastsimcoal2* and calculates 95% CIs.
 * **TajimaD_script.R:** Reads in .csv files containing Tajima's D output from VCFtools and calculates the mean (+ SE) Tajima's D in each site and with all sites pooled. Also contains code for Mann-Whitney U-tests and plots of empirical cumulative distribution functions (outlier vs. all transcripts).
 * **TD_v_pi.R:** Creates plots of Tajima's D v. pi for each sampling site and for all sampling sites pooled.
 * **Diversity_Script.R:** Calculates Ho, He and Fis from allele frequencies.
 * **Write_Sim_Data_Het_CIs.R:** Writes simulation populations and bootstraps to create CIs for Ho, He and Fis.
 * **pi.R:** Reads in .csv files containing site pi output from VCFtools and calculates the mean (+ bootstrapped CI) for pi in each site and with all sites pooled.
 * **relatedness.R:** Calculates pairwise relatedness (point estimates for all possible pairs, mean within-site relatedness, and CI for within-site relatedness).
 
 *If input files for a script were created by running code and/or calling programs on a remote workstation, details on the code used can be found in /Scripts/Upstream_Analyses.*
