################################################ Create Pi Chart of Functional Annotations ##################################################

#Pi chart of functional annotations for candidate genes in transcriptome proj
#annotations pulled from GO terms (BLAST best hit for contig of candidate SNP) & literature

################################################################################################################################################

######## Set-up ########

remove(list = ls())

#set working directory
setwd("C:/Users/Rene/Dropbox/Pinsky_Lab/Transcriptome_Proj/R_scripts/A_clarkii_transcriptomics/")
getwd()

#load libraries
library(tidyverse)
library(scales)
library(RColorBrewer)
library(ggsci)

#build dataframe
gene_function <- c("cell death", "cell structure", "metabolism", "protein folding", "protein targeting", "protein degradation", 
                   "response to stress", "transcription", "translation", "signal transduction")
number_involved <- c(2, 3, 3, 5, 4, 4, 1, 2, 6, 3)
function_prop <- c("6%", "9%", "9%", "15%", "12%", "12%", "4%", "6%", "18%", "9%")
annotation_df <- data.frame(gene_function, number_involved, function_prop)
str(annotation_df)

############################################################################################################################################

######## Create pi chart ########

function_pie <- ggplot(annotation_df, aes(x = "", y = number_involved, fill = gene_function)) + 
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start = 0)
function_pie_annotated <- function_pie + theme_minimal() + 
  scale_fill_brewer(palette = "RdYlBu", name = "Function") + 
  theme(axis.title = element_blank(), panel.border = element_blank(), panel.grid = element_blank(),
        axis.ticks = element_blank(), plot.title = element_blank(), axis.text = element_blank(),
        legend.justification = c(1, 0), legend.text = element_text(size = 35), legend.title = element_text(size = 40, face = "bold"))
function_pie_annotated
