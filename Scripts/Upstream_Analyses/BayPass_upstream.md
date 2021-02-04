BayPass\_upstream
================

Code for BayPass analyses run on an HPC cluster (Amarel at Rutgers) before transferring output files to local computer for downstream analyses and visualization in R.

Downloaded BAYPASS v.2.1 to Amarel and unzipped/installed in `~/Programs`.

Created BAYPASS input file with PGDSpider v.2.1.1.5 & modifications on local computer.

-   Took `output.hicov2.snps.only.mac2.vcf` file and converted to `clownfish_mac2.snps`.

Used following code to run PGDSpider v.2.1.1.5 on Amarel:

``` bash
#grabbed an interactive node with following line of code
srun --partition=main --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=2000 --time=00:30:00 --export=ALL --pty bash -i

#called #PGDSpider in ~/Programs/PGDSpider_2.1.1.5
java -Xmx1024m -Xms512m -jar PGDSpider2-cli.jar -inputfile output.hicov2.snps.only.mac2.vcf -inputformat VCF -outputfile clownfish_mac2.snps -outputformat BAYENV -spid VCF_BAYENV.spid

#NOTE: if don't have spid file made, or know proper format, running the previous line of code without the -spid argument will generate a template spid file to modify as needed
```

Spid format for PGDSpider:

``` bash
# VCF Parser questions
PARSER_FORMAT=VCF

# Only output SNPs with a phred-scaled quality of at least:
VCF_PARSER_QUAL_QUESTION=
# Select population definition file:
VCF_PARSER_POP_FILE_QUESTION= Pop_assignments.txt #txt file with indv names in one column and pop assignment (Pop_1, etc.) in next (no headers, tab-delimited)
# What is the ploidy of the data?
VCF_PARSER_PLOIDY_QUESTION=  DIPLOID
# Do you want to include a file with population definitions?
VCF_PARSER_POP_QUESTION= TRUE
# Output genotypes as missing if the phred-scale genotype quality is below:
VCF_PARSER_GTQUAL_QUESTION=
# Do you want to include non-polymorphic SNPs?
VCF_PARSER_MONOMORPHIC_QUESTION= TRUE
# Only output following individuals (ind1, ind2, ind4, ...):
VCF_PARSER_IND_QUESTION=
# Only input following regions (refSeqName:start:end, multiple regions: whitespace separated):
VCF_PARSER_REGION_QUESTION=
# Output genotypes as missing if the read depth of a position for the sample is below:
VCF_PARSER_READ_QUESTION=
# Take most likely genotype if "PL" or "GL" is given in the genotype field?
VCF_PARSER_PL_QUESTION= FALSE
# Do you want to exclude loci with only missing data?
VCF_PARSER_EXC_MISSING_LOCI_QUESTION= FALSE

# BAYENV Writer questions
WRITER_FORMAT=BAYENV

# Save sample file
BAYENV_WRITER_SAMPLE_FILE_QUESTION=
# Assign half missing genotypes (one allele missing) as complete missing?
BAYENV_WRITER_HALF_MISSING_QUESTION=
# Do you want to save two additional files with used sample and loci names?
BAYENV_WRITER_WRITE_INFO_FILE_QUESTION=
# Do you want to save an additional sample file with sample sizes?
BAYENV_WRITER_WRITE_SAMPLE_FILE_QUESTION=
# Save sample/loci names file
BAYENV_WRITER_INFO_FILE_QUESTION=
```

Currently, there is no option to convert a VCF file (or any other file format) to the input file format BAYPASS needs. Best option is to create input file for BayEnv (similar program to BAYPASS) and then convert to BAYPASS using either regex or on local computer.

-   To convert to BAYPASS input, need to turn 3 columns (each pop has 1 column) into 6 columns (each pop has 2 columns). BayEnv input file has every allele as one row (e.g., the first row is the count for the first allele at the first SNP in all three pops, the second row is the cound for the second allele at the first SNP, etc.). The BayPass input file has every SNP as one row (e.g., the first row is the count for the first and second allele at the first SNP in all three pops - thus the 2 columns for each pop).

Transformed BayEnv input file to BAYPASS input file titled `clownfish_mac2.geno`.

Next, standardized env variables of interest (SST mean, SST min, SST max, lat, SSS mean). All env variables (except for latitude) were pulled from the MARSPEC GIS database for each sampling location.

To standardize, ran following code in R:

``` r
#example for SSS mean

#create vector with raw data
SSS_mean_raw <- c(34.38, 34.09, 33.66)
SSS_mean_st <- scale(SSS_mean_raw)

#done for each env of interest
#NOTE: SST min & SST max were put in one vector and standardized together because they use same original scale (degrees celsius)
```

Once standardized, values put into a separate `*.env` text file for each env variable (3 columns, 1 row - each pop in own column in same order as in `*.geno` file).

-   Ex: `clownfish_data_ssmean.env`

Created `BayPass.sbatch` file to run BAYPASS:

``` bash
#!/bin/bash

#SBATCH --partition=p_mlp195 #p_mlp195 for pinsky lab partition, otherwise use main (or EOAS or E&E)
#SBATCH --job-name=baypass
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8000
#SBATCH --time=2:00:00

#SBATCH --requeue
#SBATCH --export=ALL
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=rene.clark@rutgers.edu

#script for running BayPass on Amarel

cd ~/Aclarkii_Transcriptome/BayPass

#uncomment if need to run core model to create covariance matrix and standardize env variables
g_baypass -npop 3 -gfile clownfish_mac2.geno -outprefix mac2core1

#uncomment if need to run auxiliary model (make sure env variables are standardized!)
g_baypass -npop 3 -gfile clownfish_mac2.geno -efile clownfish_data_sssmean.env -auxmodel -omegafile clownfish_mat.cov -outprefix mac2aux1_ssmean
```

Ran with default settings: burn-in of 5,000 iterations followed by an additional 25,000 MCMC steps with thinning every 25 iterations (for more info/options, call `g_baypass -help`).

First, ran just core model to (1) create covariance matrix for the standard covariate model (aux model) and (2) generate XtX statistics. BAYPASS generates many output files, but the useful ones here are:

1.  `mac2core1_mat_omega.out`: This is the covariance matrix needed to run aux model (file ending needs to be changed to `*.cov`.)
2.  `mac2core1_summary_pi_xtx.out`: This has the XtX variables for each SNP.

Then, ran the axuiliary model to generate Bayes Factors for each env variable/SNP combination. Ran model for each env variable separately, to avoid issues with correlation. Again, got a lot of output files, but the `*_summary_betai.out` files contain the Bayes Factors.

-   Ran 2X for each env variable to make sure model converged (got same Bayes Factors each time).

Copied `*_summary_betai.out` files to local computer and read into R for downstream analyses (`Scripts/Pull_BFs.R`).

**To create XtX threshold value**

After the original core model was run (with raw data), copied `clownfish_mac2_mat.cov`, `mac2core1_summary_beta_params.out`, and `clownfish_mac2.geno` to local computer and read into R (`Scripts/XtX_Calibration.R`) to generate pseudo-observed data sets (PODS) under a null model of no selection.

-   Used functions from `baypass_utils.R` script that comes with BAYPASS download/installation (found in `baypass_2.1/utils`).

Took `*.geno` file that was generated using the `simulate.baypass()` function and re-ran core model in BAYPASS with that as input. Then, copied the `*_summary_pi_xtx.out` file created from that BAYPASS run, and copied it to local computer. Read into R, and used the XtX values in the file to generate an empirical distribution of XtX under no selection. Used the 99% quantile of this distribution as the selection/neutrality threshold. SNPs with XtX values above this cut-off (from raw data) are considered adaptively differentiated.

-   Code for this in `Scripts/XtX_Calibration.R`.

**To assess SNP-environment associations**

To determine if we had more significant SNP-environmental associations than expected by chance, we created ten permuted datasets to represent the null hypothesis of no associatio nbetween allele frequencies and environmental conditions. These datasets by randomly reassigning individuals to one of the three sampling locations (`Scripts/Write Simulation Pop Data for BayPass.R` --&gt; creates BayEnv-formatted file which needs to be changed to BAYPASS input format).

Once input files were created, ran the permuted datasets through BAYPASS in the same manner as the raw data (first ran core model to create covariance matrix, then ran auxiliary model separately for each env variable).

Copied the `*_summary_betai.out` files to local computer, and calculated the mean Bayes Factor for each unique SSNP-environmental covariable combination across all ten permuted datasets to create a null distribution of BFs to compare wit the resutls from the raw data. These values were combined with the Bayes Factors from the raw data in `Data/M-W_test_BF_data.csv`, which was read into `Scripts/ECDFS for Sim v Real BFs.R` to compare the two distributions.
