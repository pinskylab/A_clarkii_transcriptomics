#################################################### Script for Creating STRUCTURE Plots  ########################################################

#using pophelper library as described in Francis 2016 - Molecular Ecology Resources
#manual found at: royfrancis.github.io/pophelper/#1_introduction
#STRUCTURE plots with all loci, loci only in HWE, outlier loci excluded, and only outlier loci

#################################################################################################################################################

######## Set-up ########

#set working directory
setwd("C:/Users/Rene/Dropbox/Pinsky_Lab/Transcriptome_Proj/R_scripts/A_clarkii_transcriptomics/")
getwd()

remove(list = ls())

#load libraries
library(devtools) #v.2.4.5
library(tidyverse) #v.2.0.0
library(gridExtra) #v.2.3
library(gtable) #v.0.3.1
#devtools::install_github('royfrancis/pophelper')
library(pophelper) #v.2.3.1

#read in data
mac2_SNPs_sfiles <- list.files(path = "Data/STRUCTURE_Output/allSNPs_results/mac2SNPs_StructureResults/", full.names = TRUE)
  mac2_SNPs_slist <- readQ(files = mac2_SNPs_sfiles, filetype = "structure") 
mac2_HWE_SNPs_sfiles <- list.files(path = "Data/STRUCTURE_Output/inhweSNPs_results/inhwe_mac2_StructureResults/", full.names = TRUE)
  mac2_HWE_SNPs_slist <- readQ(files = mac2_HWE_SNPs_sfiles, filetype = "structure")
nooutlier_mac2_sfiles <- list.files(path = "Data/STRUCTURE_Output/nooutlierswRDA_results/nooutlierswRDA_mac2_StructureResults/", full.names = TRUE)
  nooutlier_mac2_slist <- readQ(files = nooutlier_mac2_sfiles, filetype = "structure")
outlier_SNPs_sfiles <- list.files(path = "Data/STRUCTURE_Output/outliersonlywRDA_results/outliersonlywRDA_mac2_StructureResults/", full.names = TRUE)
  outlier_SNPs_slist <- readQ(files = outlier_SNPs_sfiles, filetype = "structure")
  
#create group labels 
grplab <- c(rep("Japan", 8), rep("Indonesia", 7), rep("Philippines", 10))#create group labels
meta.data <- data.frame(loc = grplab)
meta.data$loc <- as.character(meta.data$loc)
  
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
mac2_SNPs_aligned_K2 <- readQ("Data/STRUCTURE_Output/allSNPs_results/mac2_CLUMPP/pop_K2-combined-aligned.txt")
mac2_SNPs_aligned_K3 <- readQ("Data/STRUCTURE_Output/allSNPs_results/mac2_CLUMPP/pop_K3-combined-aligned.txt")
mac2_SNPs_aligned_K4 <- readQ("Data/STRUCTURE_Output/allSNPs_results/mac2_CLUMPP/pop_K4-combined-aligned.txt")
mac2_SNPs_aligned_K5 <- readQ("Data/STRUCTURE_Output/allSNPs_results/mac2_CLUMPP/pop_K5-combined-aligned.txt")

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
mac2_HWE_SNPs_aligned_K2 <- readQ("Data/STRUCTURE_Output/inhweSNPs_results/inhwe_mac2_CLUMPP/pop_K2-combined-aligned.txt")
mac2_HWE_SNPs_aligned_K3 <- readQ("Data/STRUCTURE_Output/inhweSNPs_results/inhwe_mac2_CLUMPP/pop_K3-combined-aligned.txt")
mac2_HWE_SNPs_aligned_K4 <- readQ("Data/STRUCTURE_Output/inhweSNPs_results/inhwe_mac2_CLUMPP/pop_K4-combined-aligned.txt")
mac2_HWE_SNPs_aligned_K5 <- readQ("Data/STRUCTURE_Output/inhweSNPs_results/inhwe_mac2_CLUMPP/pop_K5-combined-aligned.txt")

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
nooutlier_mac2_SNPs_aligned_K2 <- readQ("Data/STRUCTURE_Output/nooutlierswRDA_results/nooutlierswRDA_mac2_CLUMPP/pop_K2-combined-aligned.txt")
nooutlier_mac2_SNPs_aligned_K3 <- readQ("Data/STRUCTURE_Output/nooutlierswRDA_results/nooutlierswRDA_mac2_CLUMPP/pop_K3-combined-aligned.txt")
nooutlier_mac2_SNPs_aligned_K4 <- readQ("Data/STRUCTURE_Output/nooutlierswRDA_results/nooutlierswRDA_mac2_CLUMPP/pop_K4-combined-aligned.txt")
nooutlier_mac2_SNPs_aligned_K5 <- readQ("Data/STRUCTURE_Output/nooutlierswRDA_results/nooutlierswRDA_mac2_CLUMPP/pop_K5-combined-aligned.txt")

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
outlier_SNPs_aligned_K2 <- readQ("Data/STRUCTURE_Output/outliersonlywRDA_results/outliersonlywRDA_mac2_CLUMPP/pop_K2-combined-aligned.txt")
outlier_SNPs_aligned_K3 <- readQ("Data/STRUCTURE_Output/outliersonlywRDA_results/outliersonlywRDA_mac2_CLUMPP/pop_K3-combined-aligned.txt")
outlier_SNPs_aligned_K4 <- readQ("Data/STRUCTURE_Output/outliersonlywRDA_results/outliersonlywRDA_mac2_CLUMPP/pop_K4-combined-aligned.txt")
outlier_SNPs_aligned_K5 <- readQ("Data/STRUCTURE_Output/outliersonlywRDA_results/outliersonlywRDA_mac2_CLUMPP/pop_K5-combined-aligned.txt")

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

##############################################################################################################################################

######## Data for noN4 plots ########
#run separately

remove(list = ls())

#read in data
noN4_mac2_SNPs_sfiles <- list.files(path = "Data/STRUCTURE_Output/noN4_results/noN4_mac2_StructureResults/", full.names = TRUE)
  noN4_mac2_SNPs_slist <- readQ(files = noN4_mac2_SNPs_sfiles, filetype = "structure") 
noN4_mac2_HWE_SNPs_sfiles <- list.files(path = "Data/STRUCTURE_Output/noN4_inhwe_results/noN4_inhwe_mac2_StructureResults/", full.names = TRUE)
  noN4_mac2_HWE_SNPs_slist <- readQ(files = noN4_mac2_HWE_SNPs_sfiles, filetype = "structure")
noN4_nooutlier_mac2_sfiles <- list.files(path = "Data/STRUCTURE_Output/noN4_nooutliersnoN4_results/noN4_nooutliersnoN4_mac2_StructureResults/", full.names = TRUE)
  noN4_nooutlier_mac2_slist <- readQ(files = noN4_nooutlier_mac2_sfiles, filetype = "structure")
noN4_outlier_SNPs_sfiles <- list.files(path = "Data/STRUCTURE_Output/noN4_outliersonlynoN4_results/noN4_outliersonlynoN4_mac2_StructureResults/", full.names = TRUE)
  noN4_outlier_SNPs_slist <- readQ(files = noN4_outlier_SNPs_sfiles, filetype = "structure")

#create group labels 
grplab <- c(rep("Japan", 8), rep("Indonesia", 6), rep("Philippines", 10))#create group labels
meta.data <- data.frame(loc = grplab)
meta.data$loc <- as.character(meta.data$loc)

################################################################################################################################################

######## noN4 SNPs w/mac > 2 Only Set-up ########

#summarize mac2_SNPs results
noN4_mac2_SNPs_table <- tabulateQ(qlist = noN4_mac2_SNPs_slist) #pulls parameters (including mean ln likelihood - mvli) from mac2_SNPs runs and puts into dataframe
noN4_mac2_SNPs_table_sum <- summariseQ(noN4_mac2_SNPs_table) #summarizes mac2_SNPs_table by K

#evanno method
noN4_mac2_SNPs_em <- evannoMethodStructure(data = noN4_mac2_SNPs_table_sum) #evanno method to pick best value of K
noN4_mac2_SNPs_em_plot <- evannoMethodStructure(data = noN4_mac2_SNPs_table_sum, exportplot = TRUE)

#clumpp
clumppExport(qlist = noN4_mac2_SNPs_slist, useexe = TRUE) #run clumpp to order clusters properly

######## Create STRUCTURE Plot ########

#read in CLUMPP data
noN4_mac2_SNPs_aligned_K2 <- readQ("Data/STRUCTURE_Output/noN4_results/noN4_mac2_CLUMPP/pop_K2-combined-aligned.txt")
noN4_mac2_SNPs_aligned_K3 <- readQ("Data/STRUCTURE_Output/noN4_results/noN4_mac2_CLUMPP/pop_K3-combined-aligned.txt")
noN4_mac2_SNPs_aligned_K4 <- readQ("Data/STRUCTURE_Output/noN4_results/noN4_mac2_CLUMPP/pop_K4-combined-aligned.txt")
noN4_mac2_SNPs_aligned_K5 <- readQ("Data/STRUCTURE_Output/noN4_results/noN4_mac2_CLUMPP/pop_K5-combined-aligned.txt")

#create plots
noN4_mac2_SNPs_K2 <- plotQ(noN4_mac2_SNPs_aligned_K2[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, 
                      clustercol = c("#2121D9", "#FF9329"),
                      showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
                      showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
                      grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
                      showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "SNPs w/mac > 2 included")

noN4_mac2_SNPs_K3 <- plotQ(noN4_mac2_SNPs_aligned_K3[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,  
                      clustercol = c("#9999FF", "#2121D9", "#FF9329"),
                      showsp = TRUE, spbgcol = "white", splab = "K = 3", splabsize = 6,
                      showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
                      grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1,
                      showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "SNPs w/mac > 2 included")

noN4_mac2_SNPs_K4 <- plotQ(noN4_mac2_SNPs_aligned_K4[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,  
                      clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23"),
                      showsp = TRUE, spbgcol = "white", splab = "K = 4", splabsize = 6,
                      showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
                      grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1,
                      showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "SNPs w/mac > 2 included")

noN4_mac2_SNPs_K5 <- plotQ(noN4_mac2_SNPs_aligned_K5[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,  
                      clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23", "#610B5E"),
                      showsp = TRUE, spbgcol = "white", splab = "K = 5", splabsize = 6,
                      showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
                      grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1,
                      showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "SNPs w/mac > 2 included")

################################################################################################################################################

######## noN4 HWE SNPs w/mac > 2 Only Set-up ########

#summarize mac2_HWE_SNPs results
noN4_mac2_HWE_SNPs_table <- tabulateQ(qlist = noN4_mac2_HWE_SNPs_slist) #pulls parameters (including mean ln likelihood - mvli) from mac2_HWE_SNPs runs and puts into dataframe
noN4_mac2_HWE_SNPs_table_sum <- summariseQ(noN4_mac2_HWE_SNPs_table) #summarizes mac2_HWE_SNPs_table by K

#evanno method
noN4_mac2_HWE_SNPs_em <- evannoMethodStructure(data = noN4_mac2_HWE_SNPs_table_sum) #evanno method to pick best value of K
noN4_mac2_HWE_SNPs_em_plot <- evannoMethodStructure(data = noN4_mac2_HWE_SNPs_table_sum, exportplot = TRUE)

#clumpp
clumppExport(qlist = noN4_mac2_HWE_SNPs_slist, useexe = TRUE) #run clumpp to order clusters properly

######## Create STRUCTURE Plot ########

#read in CLUMPP data
noN4_mac2_HWE_SNPs_aligned_K2 <- readQ("Data/STRUCTURE_Output/noN4_inhwe_results/noN4_inhwe_mac2_CLUMPP/pop_K2-combined-aligned.txt")
noN4_mac2_HWE_SNPs_aligned_K3 <- readQ("Data/STRUCTURE_Output/noN4_inhwe_results/noN4_inhwe_mac2_CLUMPP/pop_K3-combined-aligned.txt")
noN4_mac2_HWE_SNPs_aligned_K4 <- readQ("Data/STRUCTURE_Output/noN4_inhwe_results/noN4_inhwe_mac2_CLUMPP/pop_K4-combined-aligned.txt")
noN4_mac2_HWE_SNPs_aligned_K5 <- readQ("Data/STRUCTURE_Output/noN4_inhwe_results/noN4_inhwe_mac2_CLUMPP/pop_K5-combined-aligned.txt")

#create plots
noN4_mac2_HWE_SNPs_K2 <- plotQ(noN4_mac2_HWE_SNPs_aligned_K2[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, 
                          clustercol = c("#2121D9", "#9999FF"),
                          showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
                          showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
                          grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
                          showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "HWE SNPs w/mac > 2 included")

noN4_mac2_HWE_SNPs_K3 <- plotQ(noN4_mac2_HWE_SNPs_aligned_K3[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,  
                          clustercol = c("#FF9329", "#9999FF", "#2121D9"),
                          showsp = TRUE, spbgcol = "white", splab = "K = 3", splabsize = 6,
                          showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
                          grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1,
                          showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "HWE SNPs w/mac > 2 included")

noN4_mac2_HWE_SNPs_K4 <- plotQ(noN4_mac2_HWE_SNPs_aligned_K4[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,  
                          clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23"),
                          showsp = TRUE, spbgcol = "white", splab = "K = 4", splabsize = 6,
                          showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
                          grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1,
                          showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "HWE SNPs w/mac > 2 included")

noN4_mac2_HWE_SNPs_K5 <- plotQ(noN4_mac2_HWE_SNPs_aligned_K5[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,  
                          clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23", "#610B5E"),
                          showsp = TRUE, spbgcol = "white", splab = "K = 5", splabsize = 6,
                          showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
                          grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1,
                          showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "HWE SNPs w/mac > 2 included")

################################################################################################################################################

######## noN4 No Outlier mac > 2 SNPs Only Set-up ########

#summarize nooutlier_mac2_SNPs results
noN4_nooutlier_mac2_SNPs_table <- tabulateQ(qlist = noN4_nooutlier_mac2_slist) #pulls parameters (including mean ln likelihood -mvli) from nooutlier_mac2_SNPs runs and puts into dataframe
noN4_nooutlier_mac2_SNPs_table_sum <- summariseQ(noN4_nooutlier_mac2_SNPs_table) #summarizes nooutlier_mac2_SNPs_table_sum by K

#evanno method
noN4_nooutlier_mac2_SNPs_em <- evannoMethodStructure(data = noN4_nooutlier_mac2_SNPs_table_sum) #evanno method to pick best value of K
noN4_nooutlier_mac2_SNPs_em_plot <- evannoMethodStructure(data = noN4_nooutlier_mac2_SNPs_table_sum, exportplot = TRUE)

#clumpp
clumppExport(qlist = noN4_nooutlier_mac2_slist, useexe = TRUE) #run clumpp to order clusters properly

######## Create STRUCTURE Plot ########

#read in CLUMPP data
noN4_nooutlier_mac2_SNPs_aligned_K2 <- readQ("Data/STRUCTURE_Output/noN4_nooutliersnoN4_results/noN4_nooutliersnoN4_mac2_CLUMPP/pop_K2-combined-aligned.txt")
noN4_nooutlier_mac2_SNPs_aligned_K3 <- readQ("Data/STRUCTURE_Output/noN4_nooutliersnoN4_results/noN4_nooutliersnoN4_mac2_CLUMPP/pop_K3-combined-aligned.txt")
noN4_nooutlier_mac2_SNPs_aligned_K4 <- readQ("Data/STRUCTURE_Output/noN4_nooutliersnoN4_results/noN4_nooutliersnoN4_mac2_CLUMPP/pop_K4-combined-aligned.txt")
noN4_nooutlier_mac2_SNPs_aligned_K5 <- readQ("Data/STRUCTURE_Output/noN4_nooutliersnoN4_results/noN4_nooutliersnoN4_mac2_CLUMPP/pop_K5-combined-aligned.txt")

#create plots
noN4_nooutlier_mac2_SNPs_K2 <- plotQ(noN4_nooutlier_mac2_SNPs_aligned_K2[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, 
                                clustercol = c("#2121D9", "#9999FF"),
                                showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
                                showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
                                grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
                                showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "no outlier SNPs included")

noN4_nooutlier_mac2_SNPs_K3 <- plotQ(noN4_nooutlier_mac2_SNPs_aligned_K3[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,  
                                clustercol = c("#FF9329", "#9999FF", "#2121D9"),
                                showsp = TRUE, spbgcol = "white", splab = "K = 3", splabsize = 6,
                                showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
                                grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1,
                                showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "no outlier SNPs included")

noN4_nooutlier_mac2_SNPs_K4 <- plotQ(noN4_nooutlier_mac2_SNPs_aligned_K4[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,  
                                clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23"),
                                showsp = TRUE, spbgcol = "white", splab = "K = 4", splabsize = 6,
                                showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
                                grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1,
                                showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "no outlier SNPs included")

noN4_nooutlier_mac2_SNPs_K5 <- plotQ(noN4_nooutlier_mac2_SNPs_aligned_K5[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,  
                                clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23", "#610B5E"),
                                showsp = TRUE, spbgcol = "white", splab = "K = 5", splabsize = 6,
                                showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
                                grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1,
                                showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "no outlier SNPs included")

################################################################################################################################################

######## noN4 Outlier SNPs Only Set-up ########

#summarize outlieronly_SNPs results
noN4_outlier_SNPs_table <- tabulateQ(qlist = noN4_outlier_SNPs_slist) #pulls parameters (including mean ln likelihood - mvli) from outlier_SNPs runs and puts into dataframe
noN4_outlier_SNPs_table_sum <- summariseQ(noN4_outlier_SNPs_table) #summarizes outlier_SNPs_table by K

#evanno method
noN4_outlier_SNPs_em <- evannoMethodStructure(data = noN4_outlier_SNPs_table_sum) #evanno method to pick best value of K
noN4_outlier_SNPs_em_plot <- evannoMethodStructure(data = noN4_outlier_SNPs_table_sum, exportplot = TRUE)

#clumpp
clumppExport(qlist = noN4_outlier_SNPs_slist, useexe = TRUE) #run clumpp to order clusters properly

######## Create STRUCTURE Plot ########

#read in CLUMPP data
noN4_outlier_SNPs_aligned_K2 <- readQ("Data/STRUCTURE_Output/noN4_outliersonlynoN4_results/noN4_outliersonlynoN4_mac2_CLUMPP/pop_K2-combined-aligned.txt")
noN4_outlier_SNPs_aligned_K3 <- readQ("Data/STRUCTURE_Output/noN4_outliersonlynoN4_results/noN4_outliersonlynoN4_mac2_CLUMPP/pop_K3-combined-aligned.txt")
noN4_outlier_SNPs_aligned_K4 <- readQ("Data/STRUCTURE_Output/noN4_outliersonlynoN4_results/noN4_outliersonlynoN4_mac2_CLUMPP/pop_K4-combined-aligned.txt")
noN4_outlier_SNPs_aligned_K5 <- readQ("Data/STRUCTURE_Output/noN4_outliersonlynoN4_results/noN4_outliersonlynoN4_mac2_CLUMPP/pop_K5-combined-aligned.txt")

#create plots
noN4_outlier_SNPs_K2 <- plotQ(noN4_outlier_SNPs_aligned_K2[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE, 
                         clustercol = c("#2121D9", "#9999FF"),
                         showsp = TRUE, spbgcol = "white", splab = "K = 2", splabsize = 6, 
                         showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5, 
                         grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1, 
                         showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "outlier SNPs included")

noN4_outlier_SNPs_K3 <- plotQ(noN4_outlier_SNPs_aligned_K3[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,  
                         clustercol = c("#2121D9", "#9999FF", "#FF9329"),
                         showsp = TRUE, spbgcol = "white", splab = "K = 3", splabsize = 6,
                         showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
                         grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1,
                         showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "outlier SNPs included")

noN4_outlier_SNPs_K4 <- plotQ(noN4_outlier_SNPs_aligned_K4[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,  
                         clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23"),
                         showsp = TRUE, spbgcol = "white", splab = "K = 4", splabsize = 6,
                         showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
                         grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1,
                         showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "outlier SNPs included")

noN4_outlier_SNPs_K5 <- plotQ(noN4_outlier_SNPs_aligned_K5[1], imgoutput = "sep", returnplot = TRUE, exportplot = TRUE,  
                         clustercol = c("#2121D9", "#9999FF", "#FF9329", "#FFFB23", "#610B5E"),
                         showsp = TRUE, spbgcol = "white", splab = "K = 5", splabsize = 6,
                         showyaxis = TRUE, showticks = FALSE, indlabsize = 4, ticksize = 0.5,
                         grplab = meta.data, linesize = 0.2, pointsize = 2, showgrplab = TRUE, grplabspacer = 0.1,
                         showtitle = TRUE, titlelab = "STRUCTURE plot", showsubtitle = TRUE, subtitlelab = "outlier SNPs included")