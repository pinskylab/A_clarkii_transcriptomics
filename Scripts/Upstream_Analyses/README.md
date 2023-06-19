Contains code for upstream analyses that create input files for some of the R scripts:
 * **BayPass_upstream.md:** Contains code for running BayPass on Amarel. Creates input files for `Pull_BFs.R` & `XtX_Calibration.R` and takes input files from `Write Simulation Pop Data for BayPass.R`.
 * **PCA_upstream.md:** Contains code for running plink on Amarel to generate eigenvalue and eigenvector inputs for PCAs. Creates input files for `PCAs.R`.
 * **STRUCTURE_upstream.md:** Contains code for running STRUCTURE on Amarel. Creates input files for `STRUCTURE_script.R`.
 * **Fastsimcoal2_upstream.md:** Contains code for running *fastsimcoal2* on Amarel. Also contains information on how to create SFS with easySFS.py. Creates best point estimates for demographic parameters and input files for `fsc_CIs.R`.
 * **StairwayPlot_upstream.md:** Contains fode for running Stairway Plot v.2 on Amarel. Creates figures with estimates of demographic histories for each sampling site.
 * **TajimasD_upstream.md:** Contains code for running SnpEff and calculating Tajima's D using VCFtools on Amarel. Also contains information on how mapping to *Amphiprion frenatus* was done. Creates input files for `TajimaD_script.R` & `TD_v_pi.R`.
 * **pi_upstream.md:** Contains code for calculating site-pi using VCFtools on Amarel. Creates input files for `pi.R` & `TD_v_pi.R`.

***Bioinformatics_Notes.txt** has code and notes from the bioinformatics pipeline.*
