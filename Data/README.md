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
* **M-W_test_BF_data.csv:** BFs for experimental dataset (raw) and permuted datasets (mean) for each SNP/environmental variable combination. Read into `ECDFs for Sim v Real BFs.R`.
* **outlier_loci_mac2.csv:** List of outlier SNPs (contig_bp). Read into `Fst_script.R`.
* **output.hicov2.snps.only.mac2_eigenvec.csv:** Eigenvectors from creating PCAs with Plink1.9 on Amarel (output imported into Excel and saved as .csv). Read into `PCAs.R`.
* **STRUCTURE_Output:** Directories with STRUCTURE & CLUMPP output files generated with different sets of SNPs (all SNPs, SNPs in HWE, outlier SNPs only, no outlier SNPs). Created by running STRUCTURE on Amarel and `STRUCTURE_script.R`. Read into `STRUCTURE_script.R`.
* **SYN_mac1_TajimasD.csv:** Tajima's D for every contig (that was able to be mapped to *A. frenatus* and contained at least one synonymous SNP). Calculated within sampling sites and with all sampling sites pooled. Created by running VCFtools on Amarel (output imported into Excel and saved as .csv). Read into `TajimaD_script.R`.
* **outlier_sequences_TD_mac1_SYN.csv:** List of contigs that contain at least one outlier SNP. Formatted to be read into `TD_v_pi.R` & `TajimaD_script.R`.
* **SYN_mac1_TajimasD_combined_full.csv** & **SYN_mac1_TajimasD_outlier_combined_full.csv:** Tajima's D for every contig in each samplign site and with all sampling sites pooled (all contigs or only contigs with at least one outlier SNP). Formatted to combine with pi data. Created with `TajimaD_script.R`. Read into `TD_v_pi.R`.
* **pi_combined_site_full_mac1_SYN_aggregate.csv:** Pi for every contig in each sampling site and with all sampling sites pooled (all contigs or only contigs with at least one outlier SNP). Formatted to combine with Tajima's D data. Created with `pi.R`. Read into `TD_v_pi.R`.
* **clownfish_mac2_geno.txt:** Allele counts for each SNP in each sampling site. Read into `Diversity_Script.R`.
* **output.hicov2.snps.only.mac2_genetix.gtx:** gtx file with genotype information. Read into `Diversity_SCript.R`.
* **Contig_length.csv:** Contig lengths (bp) for every contig. Read into `pi.R`.
* **outlier_contigs_mac2.csv:** List of contigs that contain at least one outlier SNP. Read into `pi.R`.
* **mac2_sitespi.csv:** Per-site pi for every SNP. Calculated within sampling sites and with all sampling sites pooled. Created by running VCFtools on Amarel (output imported into Excel and saved as .csv). Read into `pi.R`.
* **SYN_mac1_pi.csv:** Pi for every contig (that was able to be mapped to A. frenatus and contained at least one synonymous SNP). Calculated within sampling sites and with all sampling sites pooled. Created by running VCFtools on Amarel (output imported into Excel and saved as .csv). Read into `pi.R`.
* **relatedness_input_mac2.txt:** Genotype information for every SNP. Read into `relatedness.R`.
