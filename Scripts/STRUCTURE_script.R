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
HWE_SNPs_sfiles <- list.files(path = "../../STRUCTURE_Output/inhweSNPs_results/inhwe/inhwe_StructureResults/", full.names = TRUE)
  HWE_SNPs_slist <- readQ(files = HWE_SNPs_sfiles, filetype = "structure")
outlier_SNPs_sfiles <- list.files(path = "../../STRUCTURE_Output/outliersonly_results/outliersonly_strict/outliersonly_strict_StructureResults/", full.names = TRUE)
  outlier_SNPs_slist <- readQ(files = outlier_SNPs_sfiles, filetype = "structure")
nooutlier_SNPs_sfiles <- list.files(path = "../../STRUCTURE_Output/nooutliers_results/nooutliers_strict/nooutliers_strict_StructureResults/", full.names = TRUE)
  nooutlier_SNPs_slist <- readQ(files = nooutlier_SNPs_sfiles, filetype = "structure")
mac2_SNPs_sfiles <- list.files(path = "../../STRUCTURE_Output/allSNPs_results/mac2/mac2SNPs_StructureResults/", full.names = TRUE)
  mac2_SNPs_slist <- readQ(files = mac2_SNPs_sfiles, filetype = "structure")
mac2_HWE_SNPs_sfiles <- list.files(path = "../../STRUCTURE_Output/inhweSNPs_results/inhwe_mac2/inhwe_mac2_StructureResults/", full.names = TRUE)
  mac2_HWE_SNPs_slist <- readQ(files = mac2_HWE_SNPs_sfiles, filetype = "structure")
nooutlier_mac2_sfiles <- list.files(path = "../../STRUCTURE_Output/nooutliers_results/nooutliers_mac2_strict/nooutliers_mac2_StructureResults/", full.names = TRUE)
  nooutlier_mac2_slist <- readQ(files = nooutlier_mac2_sfiles, filetype = "structure")
mac1_SNPs_sfiles <- list.files(path = "../../STRUCTURE_Output/allSNPs_results/mac1/mac1SNPS_StructureResults/", full.names= TRUE)
  mac1_SNPs_slist <- readQ(files = mac1_SNPs_sfiles, filetype = "structure")
mac1_HWE_SNPs_sfiles <- list.files(path = "../../STRUCTURE_Output/inhweSNPs_results/inhwe_mac1/inhwe_mac1_StructureResults/", full.names = TRUE)
  mac1_HWE_SNPs_slist <- readQ(files = mac1_HWE_SNPs_sfiles, filetype = "structure")
outlier_mac1_sfiles <- list.files(path = "../../STRUCTURE_Output/outliersonly_results/outliersonly_mac1/outliersonly_mac1_StructureResults/", full.names = TRUE)
  outlier_mac1_slist <- readQ(files = outlier_mac1_sfiles, filetype = "structure")
nooutlier_mac1_sfiles <- list.files(path = "../../STRUCTURE_Output/nooutliers_results/nooutliers_mac1/nooutliers_mac1_StructureResults/", full.names = TRUE)
  nooutlier_mac1_slist <- readQ(files = nooutlier_mac1_sfiles, filetype = "structure")
  
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

################################################################################################################################################

######## SNPs w/mac > 2 Only Set-up ########

#summarize mac2_SNPs results
mac2_SNPs_table <- tabulateQ(qlist = mac2_SNPs_slist) #pulls parameters (including mean ln likelihood - mvli) from mac2_SNPs runs and puts into dataframe
mac2_SNPs_table_sum <- summariseQ(mac2_SNPs_table) #summarizes mac2_SNPs_table by K

#evanno method
mac2_SNPs_em <- evannoMethodStructure(data = mac2_SNPs_table_sum) #evanno method to pick best value of K
mac2_SNPs_em_plot <- evannoMethodStructure(data = mac2_SNPs_table_sum, exportplot = TRUE)

#clumpp
clumppExport(qlist = mac2_SNPs_slist, useexe = TRUE) #run clumpp to order clusters properly

######## Create STRUCTURE Plot ########

#read in CLUMPP data
mac2_SNPs_aligned_K2 <- readQ("../../STRUCTURE_Output/allSNPs_results/mac2/mac2SNPs_CLUMPP/K2/pop_K2-combined-aligned.txt")
mac2_SNPs_aligned_K3 <- readQ("../../STRUCTURE_Output/allSNPs_results/mac2/mac2SNPs_CLUMPP/K3/pop_K3-combined-aligned.txt")
mac2_SNPs_aligned_K4 <- readQ("../../STRUCTURE_Output/allSNPs_results/mac2/mac2SNPs_CLUMPP/K4/pop_K4-combined-aligned.txt")
mac2_SNPs_aligned_K5 <- readQ("../../STRUCTURE_Output/allSNPs_results/mac2/mac2SNPs_CLUMPP/K5/pop_K5-combined-aligned.txt")

#create plots
mac2_SNPs_K2 <- plotQ(mac2_SNPs_aligned_K2[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, 
                           clustercol = c("#2121D9", "#FF9329"),
                           showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
                           showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
                           grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
                           showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "SNPs w/mac > 2 included")

mac2_SNPs_K3 <- plotQ(mac2_SNPs_aligned_K3[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,  
                           clustercol = c("#9999FF", "#2121D9", "#FF9329"),
                           showsp = TRUE, spbgcol = "white", splab = "K = 3", splabsize = 6,
                           showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
                           grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1,
                           showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "SNPs w/mac > 2 included")

mac2_SNPs_K4 <- plotQ(mac2_SNPs_aligned_K4[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,  
                           clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23"),
                           showsp = TRUE, spbgcol = "white", splab = "K = 4", splabsize = 6,
                           showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
                           grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1,
                           showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "SNPs w/mac > 2 included")

mac2_SNPs_K5 <- plotQ(mac2_SNPs_aligned_K5[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,  
                           clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23", "#610B5E"),
                           showsp = TRUE, spbgcol = "white", splab = "K = 5", splabsize = 6,
                           showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
                           grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1,
                           showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "SNPs w/mac > 2 included")

################################################################################################################################################

######## HWE SNPs w/mac > 2 Only Set-up ########

#summarize mac2_HWE_SNPs results
mac2_HWE_SNPs_table <- tabulateQ(qlist = mac2_HWE_SNPs_slist) #pulls parameters (including mean ln likelihood - mvli) from mac2_HWE_SNPs runs and puts into dataframe
mac2_HWE_SNPs_table_sum <- summariseQ(mac2_HWE_SNPs_table) #summarizes mac2_HWE_SNPs_table by K

#evanno method
mac2_HWE_SNPs_em <- evannoMethodStructure(data = mac2_HWE_SNPs_table_sum) #evanno method to pick best value of K
mac2_HWE_SNPs_em_plot <- evannoMethodStructure(data = mac2_HWE_SNPs_table_sum, exportplot = TRUE)

#clumpp
clumppExport(qlist = mac2_HWE_SNPs_slist, useexe = TRUE) #run clumpp to order clusters properly

######## Create STRUCTURE Plot ########

#read in CLUMPP data
mac2_HWE_SNPs_aligned_K2 <- readQ("../../STRUCTURE_Output/inhweSNPs_results/inhwe_mac2/inhwe_mac2_CLUMPP/K2/pop_K2-combined-aligned.txt")
mac2_HWE_SNPs_aligned_K3 <- readQ("../../STRUCTURE_Output/inhweSNPs_results/inhwe_mac2/inhwe_mac2_CLUMPP/K3/pop_K3-combined-aligned.txt")
mac2_HWE_SNPs_aligned_K4 <- readQ("../../STRUCTURE_Output/inhweSNPs_results/inhwe_mac2/inhwe_mac2_CLUMPP/K4/pop_K4-combined-aligned.txt")
mac2_HWE_SNPs_aligned_K5 <- readQ("../../STRUCTURE_Output/inhweSNPs_results/inhwe_mac2/inhwe_mac2_CLUMPP/K5/pop_K5-combined-aligned.txt")

#create plots
mac2_HWE_SNPs_K2 <- plotQ(mac2_HWE_SNPs_aligned_K2[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, 
                      clustercol = c("#2121D9", "#9999FF"),
                      showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
                      showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
                      grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
                      showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "HWE SNPs w/mac > 2 included")

mac2_HWE_SNPs_K3 <- plotQ(mac2_HWE_SNPs_aligned_K3[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,  
                      clustercol = c("#FF9329", "#9999FF", "#2121D9"),
                      showsp = TRUE, spbgcol = "white", splab = "K = 3", splabsize = 6,
                      showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
                      grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1,
                      showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "HWE SNPs w/mac > 2 included")

mac2_HWE_SNPs_K4 <- plotQ(mac2_HWE_SNPs_aligned_K4[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,  
                      clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23"),
                      showsp = TRUE, spbgcol = "white", splab = "K = 4", splabsize = 6,
                      showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
                      grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1,
                      showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "HWE SNPs w/mac > 2 included")

mac2_HWE_SNPs_K5 <- plotQ(mac2_HWE_SNPs_aligned_K5[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,  
                      clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23", "#610B5E"),
                      showsp = TRUE, spbgcol = "white", splab = "K = 5", splabsize = 6,
                      showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
                      grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1,
                      showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "HWE SNPs w/mac > 2 included")

################################################################################################################################################

######## No Outlier mac > 2 SNPs Only Set-up ########

#summarize nooutlier_mac2_SNPs results
nooutlier_mac2_SNPs_table <- tabulateQ(qlist = nooutlier_mac2_slist) #pulls parameters (including mean ln likelihood -mvli) from nooutlier_mac2_SNPs runs and puts into dataframe
nooutlier_mac2_SNPs_table_sum <- summariseQ(nooutlier_mac2_SNPs_table) #summarizes nooutlier_mac2_SNPs_table_sum by K

#evanno method
nooutlier_mac2_SNPs_em <- evannoMethodStructure(data = nooutlier_mac2_SNPs_table_sum) #evanno method to pick best value of K
nooutlier_mac2_SNPs_em_plot <- evannoMethodStructure(data = nooutlier_mac2_SNPs_table_sum, exportplot = TRUE)

#clumpp
clumppExport(qlist = nooutlier_mac2_slist, useexe = TRUE) #run clumpp to order clusters properly

######## Create STRUCTURE Plot ########

#read in CLUMPP data
nooutlier_mac2_SNPs_aligned_K2 <- readQ("../../STRUCTURE_Output/nooutliers_results/nooutliers_mac2_strict/nooutliers_mac2_CLUMPP/K2/pop_K2-combined-aligned.txt")
nooutlier_mac2_SNPs_aligned_K3 <- readQ("../../STRUCTURE_Output/nooutliers_results/nooutliers_mac2_strict/nooutliers_mac2_CLUMPP/K3/pop_K3-combined-aligned.txt")
nooutlier_mac2_SNPs_aligned_K4 <- readQ("../../STRUCTURE_Output/nooutliers_results/nooutliers_mac2_strict/nooutliers_mac2_CLUMPP/K4/pop_K4-combined-aligned.txt")
nooutlier_mac2_SNPs_aligned_K5 <- readQ("../../STRUCTURE_Output/nooutliers_results/nooutliers_mac2_strict/nooutliers_mac2_CLUMPP/K5/pop_K5-combined-aligned.txt")

#create plots
nooutlier_mac2_SNPs_K2 <- plotQ(nooutlier_mac2_SNPs_aligned_K2[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, 
                           clustercol = c("#2121D9", "#9999FF"),
                           showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
                           showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
                           grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
                           showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "no outlier SNPs included")

nooutlier_mac2_SNPs_K3 <- plotQ(nooutlier_mac2_SNPs_aligned_K3[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,  
                           clustercol = c("#FF9329", "#9999FF", "#2121D9"),
                           showsp = TRUE, spbgcol = "white", splab = "K = 3", splabsize = 6,
                           showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
                           grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1,
                           showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "no outlier SNPs included")

nooutlier_mac2_SNPs_K4 <- plotQ(nooutlier_mac2_SNPs_aligned_K4[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,  
                           clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23"),
                           showsp = TRUE, spbgcol = "white", splab = "K = 4", splabsize = 6,
                           showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
                           grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1,
                           showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "no outlier SNPs included")

nooutlier_mac2_SNPs_K5 <- plotQ(nooutlier_mac2_SNPs_aligned_K5[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,  
                           clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23", "#610B5E"),
                           showsp = TRUE, spbgcol = "white", splab = "K = 5", splabsize = 6,
                           showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
                           grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1,
                           showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "no outlier SNPs included")

################################################################################################################################################

######## SNPs w/mac > 1 Only Set-up ########

#summarize mac1_SNPs results
mac1_SNPs_table <- tabulateQ(qlist = mac1_SNPs_slist) #pulls parameters (including mean ln likelihood -mvli) from mac1_SNPs runs and puts into dataframe
mac1_SNPs_table_sum <- summariseQ(mac1_SNPs_table) #summarizes mac1_SNPs_table by K

#evanno method
mac1_SNPs_em <- evannoMethodStructure(data = mac1_SNPs_table_sum) #evanno method to pick best value of K
mac1_SNPs_em_plot <- evannoMethodStructure(data = mac1_SNPs_table_sum, exportplot = TRUE)

#clumpp
clumppExport(qlist = mac1_SNPs_slist, useexe = TRUE) #run clumpp to order clusters properly

######## Create STRUCTURE Plot ########

#read in CLUMPP data
mac1_SNPs_aligned_K2 <- readQ("../../STRUCTURE_Output/allSNPs_results/mac1/mac1SNPS_CLUMPP/K2/pop_K2-combined-aligned.txt")
mac1_SNPs_aligned_K3 <- readQ("../../STRUCTURE_Output/allSNPs_results/mac1/mac1SNPS_CLUMPP/K3/pop_K3-combined-aligned.txt")
mac1_SNPs_aligned_K4 <- readQ("../../STRUCTURE_Output/allSNPs_results/mac1/mac1SNPS_CLUMPP/K4/pop_K4-combined-aligned.txt")
mac1_SNPs_aligned_K5 <- readQ("../../STRUCTURE_Output/allSNPs_results/mac1/mac1SNPS_CLUMPP/K5/pop_K5-combined-aligned.txt")

#create plots
mac1_SNPs_K2 <- plotQ(mac1_SNPs_aligned_K2[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,
                      clustercol = c("#2121D9", "#9999FF"), 
                      showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6,
                      showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
                      grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
                      showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "SNPs w/mac > 1")

mac1_SNPs_K3 <- plotQ(mac1_SNPs_aligned_K3[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,
                      clustercol = c("#FF9329", "#9999FF", "#2121D9"), 
                      showsp = TRUE, spbgcol = "white", splab = "K = 3", splabsize = 6,
                      showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
                      grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
                      showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "SNPs w/mac > 1")

mac1_SNPs_K4 <- plotQ(mac1_SNPs_aligned_K4[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,
                      clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23"), 
                      showsp = TRUE, spbgcol = "white", splab = "K = 4", splabsize = 6,
                      showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
                      grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
                      showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "SNPs w/mac > 1")

mac1_SNPs_K5 <- plotQ(mac1_SNPs_aligned_K5[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,
                      clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23", "#610B5E"), 
                      showsp = TRUE, spbgcol = "white", splab = "K = 5", splabsize = 6,
                      showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
                      grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
                      showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "SNPs w/mac > 1")

################################################################################################################################################

######## HWE SNPs w/mac > 1 Only Set-up ########

#summarize mac1_HWE_SNPs results
mac1_HWE_SNPs_table <- tabulateQ(qlist = mac1_HWE_SNPs_slist) #pulls parameters (including mean ln likelihood -mvli) from mac1_HWE_SNPs runs and puts into dataframe
mac1_HWE_SNPs_table_sum <- summariseQ(mac1_HWE_SNPs_table) #summarizes mac1_HWE_SNPs_table by K

#evanno method
mac1_HWE_SNPs_em <- evannoMethodStructure(data = mac1_HWE_SNPs_table_sum) #evanno method to pick best value of K
mac1_HWE_SNPs_em_plot <- evannoMethodStructure(data = mac1_HWE_SNPs_table_sum, exportplot = TRUE)

#clumpp
clumppExport(qlist = mac1_HWE_SNPs_slist, useexe = TRUE) #run clumpp to order clusters properly

######## Create STRUCTURE Plot ########

#read in CLUMPP data
inhwe_mac1_SNPs_aligned_K2 <- readQ("../../STRUCTURE_Output/inhweSNPs_results/inhwe_mac1/inhwe_mac1_CLUMPP/K2/pop_K2-combined-aligned.txt")
inhwe_mac1_SNPs_aligned_K3 <- readQ("../../STRUCTURE_Output/inhweSNPs_results/inhwe_mac1/inhwe_mac1_CLUMPP/K3/pop_K3-combined-aligned.txt")
inhwe_mac1_SNPs_aligned_K4 <- readQ("../../STRUCTURE_Output/inhweSNPs_results/inhwe_mac1/inhwe_mac1_CLUMPP/K4/pop_K4-combined-aligned.txt")
inhwe_mac1_SNPs_aligned_K5 <- readQ("../../STRUCTURE_Output/inhweSNPs_results/inhwe_mac1/inhwe_mac1_CLUMPP/K5/pop_K5-combined-aligned.txt")

#create plots
inhwe_mac1_SNPs_K2 <- plotQ(inhwe_mac1_SNPs_aligned_K2[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,
                      clustercol = c("#9999FF", "#2121D9"), 
                      showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6,
                      showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
                      grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
                      showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "HWE SNPs w/mac > 1")

inhwe_mac1_SNPs_K3 <- plotQ(inhwe_mac1_SNPs_aligned_K3[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,
                      clustercol = c("#FF9329", "#9999FF", "#2121D9"), 
                      showsp = TRUE, spbgcol = "white", splab = "K = 3", splabsize = 6,
                      showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
                      grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
                      showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "HWE SNPs w/mac > 1")

inhwe_mac1_SNPs_K4 <- plotQ(inhwe_mac1_SNPs_aligned_K4[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,
                      clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23"), 
                      showsp = TRUE, spbgcol = "white", splab = "K = 4", splabsize = 6,
                      showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
                      grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
                      showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "HWE SNPs w/mac > 1")

inhwe_mac1_SNPs_K5 <- plotQ(inhwe_mac1_SNPs_aligned_K5[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,
                      clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23", "#610B5E"), 
                      showsp = TRUE, spbgcol = "white", splab = "K = 5", splabsize = 6,
                      showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
                      grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
                      showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "HWE SNPs w/mac > 1")

################################################################################################################################################

######## No Outlier mac > 1 SNPs Only Set-up ########

#summarize nooutlier_mac1_SNPs results
nooutlier_mac1_table <- tabulateQ(qlist = nooutlier_mac1_slist) #pulls parameters (including mean ln likelihood -mvli) from nooutlier_mac1 runs and puts into dataframe
nooutlier_mac1_table_sum <- summariseQ(nooutlier_mac1_table) #summarizes nooutlier_mac1_table by K

#evanno method
nooutlier_mac1_em <- evannoMethodStructure(data = nooutlier_mac1_table_sum) #evanno method to pick best value of K
nooutlier_mac1_em_plot <- evannoMethodStructure(data = nooutlier_mac1_table_sum, exportplot = TRUE)

#clumpp
clumppExport(qlist = nooutlier_mac1_slist, useexe = TRUE) #run clumpp to order clusters properly

######## Create STRUCTURE Plot ########

#read in CLUMPP data
nooutlier_mac1_SNPs_aligned_K2 <- readQ("../../STRUCTURE_Output/nooutliers_results/nooutliers_mac1/nooutliers_mac1_CLUMPP/K2/pop_K2-combined-aligned.txt")
nooutlier_mac1_SNPs_aligned_K3 <- readQ("../../STRUCTURE_Output/nooutliers_results/nooutliers_mac1/nooutliers_mac1_CLUMPP/K3/pop_K3-combined-aligned.txt")
nooutlier_mac1_SNPs_aligned_K4 <- readQ("../../STRUCTURE_Output/nooutliers_results/nooutliers_mac1/nooutliers_mac1_CLUMPP/K4/pop_K4-combined-aligned.txt")
nooutlier_mac1_SNPs_aligned_K5 <- readQ("../../STRUCTURE_Output/nooutliers_results/nooutliers_mac1/nooutliers_mac1_CLUMPP/K5/pop_K5-combined-aligned.txt")

#create plots
nooutlier_mac1_SNPs_K2 <- plotQ(nooutlier_mac1_SNPs_aligned_K2[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,
                            clustercol = c("#2121D9", "#9999FF"), 
                            showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6,
                            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
                            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
                            showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "No Outlier SNPs w/mac > 1")

nooutlier_mac1_SNPs_K3 <- plotQ(nooutlier_mac1_SNPs_aligned_K3[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,
                            clustercol = c("#FF9329", "#2121D9", "#9999FF"), 
                            showsp = TRUE, spbgcol = "white", splab = "K = 3", splabsize = 6,
                            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
                            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
                            showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "No Outlier SNPs w/mac > 1")

nooutlier_mac1_SNPs_K4 <- plotQ(nooutlier_mac1_SNPs_aligned_K4[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,
                            clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23"), 
                            showsp = TRUE, spbgcol = "white", splab = "K = 4", splabsize = 6,
                            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
                            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
                            showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "No Outlier SNPs w/mac > 1")

nooutlier_mac1_SNPs_K5 <- plotQ(nooutlier_mac1_SNPs_aligned_K5[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,
                            clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23", "#610B5E"), 
                            showsp = TRUE, spbgcol = "white", splab = "K = 5", splabsize = 6,
                            showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
                            grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
                            showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "No Outlier SNPs w/mac > 1")

################################################################################################################################################

######## Outlier Only mac > 1 SNPs Only Set-up ########

#summarize outlieronly_mac1_SNPs results
outlier_mac1_table <- tabulateQ(qlist = outlier_mac1_slist) #pulls parameters (including mean ln likelihood -mvli) from outlier_mac1 runs and puts into dataframe
outlier_mac1_table_sum <- summariseQ(outlier_mac1_table) #summarizes outlier_mac1_table by K

#evanno method
outlier_mac1_em <- evannoMethodStructure(data = outlier_mac1_table_sum) #evanno method to pick best value of K
outlier_mac1_em_plot <- evannoMethodStructure(data = outlier_mac1_table_sum, exportplot = TRUE)

#clumpp
clumppExport(qlist = outlier_mac1_slist, useexe = TRUE) #run clumpp to order clusters properly

######## Create STRUCTURE Plot ########

#read in CLUMPP data
outlier_mac1_SNPs_aligned_K2 <- readQ("../../STRUCTURE_Output/outliersonly_results/outliersonly_mac1/outliersonly_mac1_CLUMPP/K2/pop_K2-combined-aligned.txt")
outlier_mac1_SNPs_aligned_K3 <- readQ("../../STRUCTURE_Output/outliersonly_results/outliersonly_mac1/outliersonly_mac1_CLUMPP/K3/pop_K3-combined-aligned.txt")
outlier_mac1_SNPs_aligned_K4 <- readQ("../../STRUCTURE_Output/outliersonly_results/outliersonly_mac1/outliersonly_mac1_CLUMPP/K4/pop_K4-combined-aligned.txt")
outlier_mac1_SNPs_aligned_K5 <- readQ("../../STRUCTURE_Output/outliersonly_results/outliersonly_mac1/outliersonly_mac1_CLUMPP/K5/pop_K5-combined-aligned.txt")

#create plots
outlier_mac1_SNPs_K2 <- plotQ(outlier_mac1_SNPs_aligned_K2[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,
                                clustercol = c("#2121D9", "#9999FF"), 
                                showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6,
                                showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
                                grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
                                showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "Outlier Only SNPs w/mac > 1")

outlier_mac1_SNPs_K3 <- plotQ(outlier_mac1_SNPs_aligned_K3[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,
                                clustercol = c("#FF9329", "#9999FF", "#2121D9"), 
                                showsp = TRUE, spbgcol = "white", splab = "K = 3", splabsize = 6,
                                showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
                                grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
                                showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "Outlier Only SNPs w/mac > 1")

outlier_mac1_SNPs_K4 <- plotQ(outlier_mac1_SNPs_aligned_K4[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,
                                clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23"), 
                                showsp = TRUE, spbgcol = "white", splab = "K = 4", splabsize = 6,
                                showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
                                grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
                                showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "Outlier Only SNPs w/mac > 1")

outlier_mac1_SNPs_K5 <- plotQ(outlier_mac1_SNPs_aligned_K5[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,
                                clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23", "#610B5E"), 
                                showsp = TRUE, spbgcol = "white", splab = "K = 5", splabsize = 6,
                                showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
                                grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
                                showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "Outlier Only SNPs w/mac > 1")