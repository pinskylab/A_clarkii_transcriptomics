Contains files either read into R scripts or uploaded to Amarel for analyses. Also contains metadata and general output files. *Files generated in R scripts not included in directory, with a few exceptions.*

Metadata and data unrelated to specific scripts:
 * **Aclarkii_metadata.csv:** Metadata for individuals (sampling location, time, etc.).
 * **annotations_combined.xlsx:** Functional and structural annotations for every outlier SNP. Also contains BFs for each environmental variable and XtX values.

Data for R scripts or upstream analyses on Amarel:
* **output.hicov2.snps.only.mac2.vcf:** VCF file filtered to contain only SNPs with minor allele count > 2. Read into `Fst_script.R`. Also used for most upstream analyses on Amarel.
* **output.hicov2.snps.only.mac2.csv:** Body of `output.hicov2.snps.only.mac2.vcf` (no header info). Read into `Write Simulation Pop Data for BayPass.R` & `Write_Sim_Data_Het_CIs.R`.
* **Loc_Names_mac2.txt:** Contig_bp names for every SNP. Read into `Diversity_Script.R`, `Fst_script.R` & `Pull_BFs.R`.
* **polarized_allele_freqs.csv:** Polarized allele frequencies for every SNP at each sampling site. Read into `Allele Freq Line Plots.R`.
* **mac2aux1_summary_betai.txt:** Output file from running BayPass auxiliary model on Amarel (experimental data). BFs for each SNP. One file for every environmental variable. Read into `Pull_BFs.R`.
* **BF_total_mac2.csv:** BFs for each SNP/environmental variable combination, as well as XtX values for each SNP. Data comes from BayPass runs on Amarel. Created with `Pull_BFs.R` (XtX values added afterwards). Read into `BF Correlation Plots.R`.
* **clownfish_mac2.geno:** Input file with genotype information for BayPass (experimental data). Read into `XtX_Calibration.R` and used for running BayPass core and auxiliary models on Amarel.
* **clownfish_data.env** & **clownfish_data_std.env:** Raw and standardized environmental variable data for each sample site. Input files for BayPass core and auxiliary models on Amarel (once separated into files for each environmental variable).
* **clownfish_mac2_mat.cov:** Input file with covariance matrix for BayPass (experimental data). Created by running BayPass core model on Amarel. Read into `XtX_Calibration.R`.
* **mac2core1_pod_mat_omega.out:** Output file from running BayPass core model on Amarel (pseudo-observed datasets (PODs)). Contains covariance matrix. Read into `XtX_Calibration.R`.
* **mac2core1_summary_beta_params.out:**Output file from running BayPass core model on Amarel (either experimental or POD data). Contains estimates for beta parameters. Read into `XtX_Calibration.R`.
* **mac2core1_summary_pi_xtx.out:** Output file from running BayPass core model on Amarel (either experimental or POD data). Contains XtX estimates for each SNP. Read into `XtX_Calibration.R`.
* **Individual_Env_Data.csv:** Environmental variable data for each individual. Read into `Write Simulation Pop Data for BayPass.R`.
