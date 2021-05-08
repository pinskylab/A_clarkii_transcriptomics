Manuscript Plots:
* **Aclarkii_RangeMarginSelection_Manuscript_Fig1.tiff:** Map of sampling sites and range extent of *A. clarkii*. Red shading indicates the relative probability of occurrence of *A. clarkii* based on Aquamaps data. Created with arcGIS. Figure 1 in manuscript.
* **Aclarkii_RangeMarginSelection_Manuscript_Fig2.tiff:** Scatterplots of XtX v. BFs for each environmental variable. Black triangles represent the 56 candidate SNPs while grey circles represent the remaining SNPs. Created with `BF Correlation Plots.R`. Figure 2 in manuscript.
* **Aclarkii_RangeMarginSelection_Manuscript_Fig3.tiff:** PCA and STRUCTURE plots generated with all loci and outlier loci only. Created with `PCAs.R` and `STRUCTURE_script.R`.

Supplemental Plots:
* **TajimasD_boxplot_mac1_SYN_wRDA.png:** Boxplot of Tajima's D for each sampling site and for all sampling sites pooled. Color-coded by outlier status. Created with `TajimaD_script.R`.
* **BF_correlation_plot_mac2.png:** Scatterplots of BFs for every combination of environmental variables. Created with `BF Correlation Plots.R`.
* **BF_permuted_distribution_plots_mac2.png:** ECDFs of BFs for experimental data and permuted data for each environmental variable. Created with `ECDFs for Sim v Real BFs.R`.
* **RDA_biplot_edited.png:** RDA biplot with the distribution of SNPs (red plus signs) and genotypes of each individual (open grey circles). Blue vectors indicate correlations between the constraining axes and the environmental (predictor) variables. Created with `RDA.R` and edited in PowerPoint.
* **AF_arranged_plot_wRDA.png:** Plot of polarized allele frequencies across sampling sites (outleir loci and non-outlier loci). Created with `Allele Freqs Line Plots.R`.
* **SYN_TajimaD_ECDFs_wRDA.png:** ECDFs of Tajima's D for outlier v. non-outlier contigs for each sampling site and for all sampling sites pooled. Created with `TajimaD_script.R`.
* **PCA_STRUCTURE_inHWE_nooutliers_mac2_wRDA.png:** PCA and STRUCTURE plots generated with only loci in HWE and without outlier loci. Created with `PCAs.R` and `STRUCTURE_script.R`.
* **PCA_STRUCTURE_noN4_all_outliersonly_mac2.png:** PCA and STRUCTURE plots generated with all loci and outlier loci only. Individual N4 is not included. Created with `PCAs.R` and `STRUCTURE_script.R`.
* **Japan_fold.final.summary.png:** Effective population size (Ne) through time for Japan (also one for the Philippines and Indonesia). The solid red line is the estimate of median Ne, and the light and dark grey shaded lines represent 95% & 75% CIs, respectively. Created with Stairway Plot v.2.
* **fastsimcoal_model.png:** Simple demographic model for *fastsimcoal2*.

Other plots:
 * **Fst_boxplot_mac2.png:** Boxplot of per-locus Fsts, separated by outlier status. Created with `Fst_script.R`.
 * **Fst_plot_mac2.png:** "Manhattan plot" of per-locus Fsts. Color-coded by outlier status. X-axis ordered by contig # (not necessarily place in genome). Created with `Fst_script.R.`.
 * **TajimasD_all_mac1_SYN.png** & **TajimasD_outliersonly_mac1_SYN.png:** Scatterplots of Tajima's D for either all contigs or outlier only contigs (outlier = contig with at least one outlier locus). Color-coded by sampling site. Created with `TajimaD_script.R`.
  * **TD_pi_plot.png:** Scatterplots of Tajima's D v. pi for each sampling site, for all sampling sites pooled, and for outlier contigs only.Color-coded by outlier status (outlier contigs only is color-coded by sampling site). Created by `TD_v_pi.R`.
 * **pi_boxplot_mac2.png:** Boxplot of pi for each sampling site and for all sampling sites pooled. Color-coded by outlier status. Created by `pi.R`.
 * **pi_plot_mac2.png:** Average pi in each sampling site with 95% CIs. Created by `pi.R`.
 * **wang_relatedness_95CI_mac2.png:** Average within-site pair-wise relatedness in each sampling site with 95% CIs. created by `relatedness.R`.
