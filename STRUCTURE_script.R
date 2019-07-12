#################################################### Script for Creating STRUCTURE Plots  ########################################################

#using pophelper library as described in Francis 2016 - Molecular Ecology Resources
#manual found at: royfrancis.github.io/pophelper/#1_introduction

#################################################################################################################################################

######## Set-up ########

#set working directory
setwd("C:/Users/Rene/Dropbox/Pinsky_Lab/Transcriptome_Proj/R_scripts/A_clarkii_transcriptomics/")
getwd()

remove(list = ls())

#load libraries
library(devtools)
library(tidyverse)
library(gridExtra)
library(gtable)
#devtools::install_github('royfrancis/pophelper')
library(pophelper)

#read in data
all_SNPs_sfiles <- list.files(path = "../../STRUCTURE_Output/allSNPs_results/allSNPs/allSNPs_StructureResults/", full.names = TRUE) #create list of file names of STRUCTURE runs with all SNPs
all_SNPs_slist <- readQ(files = all_SNPs_sfiles, filetype = "structure") #read in all STRUCTURE files with all SNPs included (creates set of dataframes)
HWE_SNPs_sfiles <- list.files(path = "../../STRUCTURE_Output/inhweSNPs_results/inhwe_StructureResults/", full.names = TRUE)
HWE_SNPs_slist <- readQ(files = HWE_SNPs_sfiles, filetype = "structure")
outlier_SNPs_sfiles <- list.files(path = "../../STRUCTURE_Output/outliersonly_results/outliersonly_strict/outliersonly_strict_StructureResults/", full.names = TRUE)
outlier_SNPs_slist <- readQ(files = outlier_SNPs_sfiles, filetype = "structure")
nooutlier_SNPs_sfiles <- list.files(path = "../../STRUCTURE_Output/nooutliers_results/nooutliers_strict/nooutliers_strict_StructureResults/", full.names = TRUE)
nooutlier_SNPs_slist <- readQ(files = nooutlier_SNPs_sfiles, filetype = "structure")

################################################################################################################################################

######## All SNPs Set-up ########

#summarize all_SNPs results
all_SNPs_table <- tabulateQ(qlist = all_SNPs_slist) #pulls parameters (including mean ln likelihood - mvli) from all_SNPs runs and puts into dataframe
all_SNPs_table_sum <- summariseQ(all_SNPs_table) #summarizes all_SNPs_table by K

#evanno method
all_SNPs_em <- evannoMethodStructure(data = all_SNPs_table_sum) #evanno method to pick best value of K
all_SNPs_em_plot <- evannoMethodStructure(data = all_SNPs_table_sum, exportplot = TRUE)

#clumpp
clumppExport(qlist = all_SNPs_slist, useexe = TRUE) #run clumpp to order clusters properly

######## Create STRUCTURE Plot ########

#read in CLUMPP data
all_SNPs_aligned_K2 <- readQ("../../STRUCTURE_Output/allSNPs_results/allSNPs/allSNPs_CLUMPP/K2/pop_K2-combined-aligned.txt")
all_SNPs_aligned_K3 <- readQ("../../STRUCTURE_Output/allSNPs_results/allSNPs/allSNPs_CLUMPP/K3/pop_K3-combined-aligned.txt")
all_SNPs_aligned_K4 <- readQ("../../STRUCTURE_Output/allSNPs_results/allSNPs/allSNPs_CLUMPP/K4/pop_K4-combined-aligned.txt")
all_SNPs_aligned_K5 <- readQ("../../STRUCTURE_Output/allSNPs_results/allSNPs/allSNPs_CLUMPP/K5/pop_K5-combined-aligned.txt")

#create group labels
grplab <- c(rep("Japan", 8), rep("Indonesia", 7), rep("Philippines", 10))#create group labels
meta.data <- data.frame(loc = grplab)
meta.data$loc <- as.character(meta.data$loc)

#create plots
all_SNPs_K2 <- plotQ(all_SNPs_aligned_K2[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, 
                     clustercol = c("#2121D9", "#9999FF"),
                     showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
                     showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
                     grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
                     showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "all SNPs included")

all_SNPs_K3 <- plotQ(all_SNPs_aligned_K3[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,  
                     clustercol = c("#2121D9", "#9999FF", "#FF9329"),
                     showsp = TRUE, spbgcol = "white", splab = "K = 3", splabsize = 6,
                     showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
                     grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1,
                     showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "all SNPs included")

all_SNPs_K4 <- plotQ(all_SNPs_aligned_K4[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,  
                     clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23"),
                     showsp = TRUE, spbgcol = "white", splab = "K = 4", splabsize = 6,
                     showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
                     grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1,
                     showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "all SNPs included")

all_SNPs_K5 <- plotQ(all_SNPs_aligned_K5[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,  
                     clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23", "#610B5E"),
                     showsp = TRUE, spbgcol = "white", splab = "K = 5", splabsize = 6,
                     showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
                     grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1,
                     showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "all SNPs included")

################################################################################################################################################

######## SNPs in HWE Set-up ########

#summarize inHWE_SNPs results
HWE_SNPs_table <- tabulateQ(qlist = HWE_SNPs_slist) #pulls parameters (including mean ln likelihood - mvli) from HWE_SNPs runs and puts into dataframe
HWE_SNPs_table_sum <- summariseQ(HWE_SNPs_table) #summarizes HWE_SNPs_table by K

#evanno method
HWE_SNPs_em <- evannoMethodStructure(data = HWE_SNPs_table_sum) #evanno method to pick best value of K
HWE_SNPs_em_plot <- evannoMethodStructure(data = HWE_SNPs_table_sum, exportplot = TRUE)

#clumpp
clumppExport(qlist = HWE_SNPs_slist, useexe = TRUE) #run clumpp to order clusters properly

######## Create STRUCTURE Plot ########

#read in CLUMPP data
HWE_SNPs_aligned_K2 <- readQ("../../STRUCTURE_Output/inhweSNPs_results/inhwe_CLUMPP/K2/pop_K2-combined-aligned.txt")
HWE_SNPs_aligned_K3 <- readQ("../../STRUCTURE_Output/inhweSNPs_results/inhwe_CLUMPP/K3/pop_K3-combined-aligned.txt")
HWE_SNPs_aligned_K4 <- readQ("../../STRUCTURE_Output/inhweSNPs_results/inhwe_CLUMPP/K4/pop_K4-combined-aligned.txt")
HWE_SNPs_aligned_K5 <- readQ("../../STRUCTURE_Output/inhweSNPs_results/inhwe_CLUMPP/K5/pop_K5-combined-aligned.txt")

#create plots
HWE_SNPs_K2 <- plotQ(HWE_SNPs_aligned_K2[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, 
                     clustercol = c("#2121D9", "#9999FF"),
                     showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
                     showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
                     grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
                     showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "SNPs in HWE included")

HWE_SNPs_K3 <- plotQ(HWE_SNPs_aligned_K3[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,  
                     clustercol = c("#2121D9", "#9999FF", "#FF9329"),
                     showsp = TRUE, spbgcol = "white", splab = "K = 3", splabsize = 6,
                     showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
                     grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1,
                     showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "SNPs in HWE included")

HWE_SNPs_K4 <- plotQ(HWE_SNPs_aligned_K4[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,  
                     clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23"),
                     showsp = TRUE, spbgcol = "white", splab = "K = 4", splabsize = 6,
                     showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
                     grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1,
                     showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "SNPS in HWE included")

HWE_SNPs_K5 <- plotQ(HWE_SNPs_aligned_K5[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,  
                     clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23", "#610B5E"),
                     showsp = TRUE, spbgcol = "white", splab = "K = 5", splabsize = 6,
                     showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
                     grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1,
                     showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "SNPs in HWE included")

################################################################################################################################################

######## Outlier SNPs Only Set-up ########

#summarize outlieronly_SNPs results
outlier_SNPs_table <- tabulateQ(qlist = outlier_SNPs_slist) #pulls parameters (including mean ln likelihood - mvli) from outlier_SNPs runs and puts into dataframe
outlier_SNPs_table_sum <- summariseQ(outlier_SNPs_table) #summarizes outlier_SNPs_table by K

#evanno method
outlier_SNPs_em <- evannoMethodStructure(data = outlier_SNPs_table_sum) #evanno method to pick best value of K
outlier_SNPs_em_plot <- evannoMethodStructure(data = outlier_SNPs_table_sum, exportplot = TRUE)

#clumpp
clumppExport(qlist = outlier_SNPs_slist, useexe = TRUE) #run clumpp to order clusters properly

######## Create STRUCTURE Plot ########

#read in CLUMPP data
outlier_SNPs_aligned_K2 <- readQ("../../STRUCTURE_Output/outliersonly_results/outliersonly_strict/outliersonly_strict_CLUMPP/K2/pop_K2-combined-aligned.txt")
outlier_SNPs_aligned_K3 <- readQ("../../STRUCTURE_Output/outliersonly_results/outliersonly_strict/outliersonly_strict_CLUMPP/K3/pop_K3-combined-aligned.txt")
outlier_SNPs_aligned_K4 <- readQ("../../STRUCTURE_Output/outliersonly_results/outliersonly_strict/outliersonly_strict_CLUMPP/K4/pop_K4-combined-aligned.txt")
outlier_SNPs_aligned_K5 <- readQ("../../STRUCTURE_Output/outliersonly_results/outliersonly_strict/outliersonly_strict_CLUMPP/K5/pop_K5-combined-aligned.txt")

#create plots
outlier_SNPs_K2 <- plotQ(outlier_SNPs_aligned_K2[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, 
                     clustercol = c("#2121D9", "#9999FF"),
                     showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
                     showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
                     grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
                     showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "outlier SNPs included")

outlier_SNPs_K3 <- plotQ(outlier_SNPs_aligned_K3[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,  
                     clustercol = c("#2121D9", "#9999FF", "#FF9329"),
                     showsp = TRUE, spbgcol = "white", splab = "K = 3", splabsize = 6,
                     showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
                     grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1,
                     showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "outlier SNPs included")

outlier_SNPs_K4 <- plotQ(outlier_SNPs_aligned_K4[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,  
                     clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23"),
                     showsp = TRUE, spbgcol = "white", splab = "K = 4", splabsize = 6,
                     showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
                     grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1,
                     showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "outlier SNPs included")

outlier_SNPs_K5 <- plotQ(outlier_SNPs_aligned_K5[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,  
                     clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23", "#610B5E"),
                     showsp = TRUE, spbgcol = "white", splab = "K = 5", splabsize = 6,
                     showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
                     grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1,
                     showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "outlier SNPs included")

################################################################################################################################################

######## No Outlier SNPs Only Set-up ########

#summarize nooutlier_SNPs results
nooutlier_SNPs_table <- tabulateQ(qlist = nooutlier_SNPs_slist) #pulls parameters (including mean ln likelihood - mvli) from nooutlier_SNPs runs and puts into dataframe
nooutlier_SNPs_table_sum <- summariseQ(nooutlier_SNPs_table) #summarizes nooutlier_SNPs_table by K

#evanno method
nooutlier_SNPs_em <- evannoMethodStructure(data = nooutlier_SNPs_table_sum) #evanno method to pick best value of K
nooutlier_SNPs_em_plot <- evannoMethodStructure(data = nooutlier_SNPs_table_sum, exportplot = TRUE)

#clumpp
clumppExport(qlist = nooutlier_SNPs_slist, useexe = TRUE) #run clumpp to order clusters properly

######## Create STRUCTURE Plot ########

#read in CLUMPP data
nooutlier_SNPs_aligned_K2 <- readQ("../../STRUCTURE_Output/nooutliers_results/nooutliers_strict/nooutliers_strict_CLUMPP/K2/pop_K2-combined-aligned.txt")
nooutlier_SNPs_aligned_K3 <- readQ("../../STRUCTURE_Output/nooutliers_results/nooutliers_strict/nooutliers_strict_CLUMPP/K3/pop_K3-combined-aligned.txt")
nooutlier_SNPs_aligned_K4 <- readQ("../../STRUCTURE_Output/nooutliers_results/nooutliers_strict/nooutliers_strict_CLUMPP/K4/pop_K4-combined-aligned.txt")
nooutlier_SNPs_aligned_K5 <- readQ("../../STRUCTURE_Output/nooutliers_results/nooutliers_strict/nooutliers_strict_CLUMPP/K5/pop_K5-combined-aligned.txt")
#create plots
nooutlier_SNPs_K2 <- plotQ(nooutlier_SNPs_aligned_K2[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, 
                         clustercol = c("#2121D9", "#9999FF"),
                         showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
                         showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
                         grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
                         showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "no outlier SNPs included")

nooutlier_SNPs_K3 <- plotQ(nooutlier_SNPs_aligned_K3[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,  
                         clustercol = c("#2121D9", "#9999FF", "#FF9329"),
                         showsp = TRUE, spbgcol = "white", splab = "K = 3", splabsize = 6,
                         showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
                         grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1,
                         showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "no outlier SNPs included")

nooutlier_SNPs_K4 <- plotQ(nooutlier_SNPs_aligned_K4[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,  
                         clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23"),
                         showsp = TRUE, spbgcol = "white", splab = "K = 4", splabsize = 6,
                         showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
                         grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1,
                         showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "no outlier SNPs included")

nooutlier_SNPs_K5 <- plotQ(nooutlier_SNPs_aligned_K5[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,  
                         clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23", "#610B5E"),
                         showsp = TRUE, spbgcol = "white", splab = "K = 5", splabsize = 6,
                         showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
                         grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1,
                         showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "no outlier SNPs included")