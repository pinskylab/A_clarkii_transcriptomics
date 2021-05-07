############################################## Perform RDA ######################################################

#Run RDA for clownfish transcriptomics
#Modifying Jennifer Hoey's script (https://github.com/jahoey/AMPHoutliers/blob/master/clownfishRDA.R)

################################################################################################################################################

######## Set-up ########

remove(list = ls())

#Set working directory
setwd("C:/Users/rclar/Dropbox/Pinsky_Lab/Transcriptome_Proj/R_scripts/A_clarkii_transcriptomics/")
getwd()

#load libraries
library(adegenet)
library(devtools)
library(hierfstat)
library(pegas)
library(vegan)
library(robust)
library(qvalue)
library(ggplot2)

#read in data
mac2_loci <- read.structure("Data/STRUCTURE_mac2.str", n.ind = 25, n.loc = 4212, onerowperind = FALSE,
                            col.lab = 1, col.pop = 2, row.marknames = 1) #run independently bc needs user input
env_data <- read.csv("Data/Env_Data.csv", header = TRUE) #individuals should be same order as in mac2_loci
  #env_data <- env_data[-12,] #remove N4

################################################################################################################################################

######## Prep data ########

#### center genomic data ####
mac2_center <- scaleGen(mac2_loci, center = TRUE, scale = FALSE) #centering the data (centering allele freqs to a mean of 0)
hist(mac2_center) #check centering worked (histogram should center on 0, it does)
sum(is.na(mac2_center)) #make sure no missing data

##### scale environmental data ####
env_scaled <- as.data.frame(scale(data.matrix(env_data[,3:7])))

#scaling SST Min & SST Max together to maintain original scale diff
temp_data <- c(env_data$SST_Min, env_data$SST_Max) #pull out un-scaled temp values into one vector
temp_data_scaled <- c(scale(temp_data))

SST_Min_scaled <- c(temp_data_scaled[1:25]) #pull SST Min from scaled vector
SST_Max_scaled <- c(temp_data_scaled[26:50]) #pull SST Max from scaled vector

env_scaled$SST_Min <- SST_Min_scaled #put back into to scaled dataframe
env_scaled$SST_Max <- SST_Max_scaled

######## Run RDA ########

#perform RDA
rda_results <- rda(mac2_center ~ SSS_Mean + SST_Mean + SST_Min + SST_Max + Latitude, 
                   data = env_scaled, scale = FALSE)
rda_results
summary(rda_results) #issue with colinearity, only showing results for SSS_Mean and one of the other variables
#rearrange order of variables to get the bi-plot scores for the other env variables as well (SST_Min, SST_Max, Latitude) so that can graph in bi-plot

#check percent variance explained
summary(eigenvals(rda_results, model = "constrained")) #RDA1 explains ~60% of total variance explained by RDAs
#PC1 & PC2 explain more variance than RDA2 so means a lot of unexplained variance leftover

#check R2 value --> how much variance is explained by explanatory variables (0.075, 7.5%) --> makes sense, expect most SNPs to be neutral and uncorrelated
R2adjusted_rda <- RsquareAdj(rda_results)$adj.r.squared #use adjusted bc otherwise biased by # of explanatory variables

#permutation tests
#generate reference distribution of RDAs to see if these ones are statistically significant
anova(rda_results, by="axis", step=1000) #looks like significant, do by axis so see each axis independently (residual is variance in PCs)

################################################################################################################################################

######## ID Outliers ########

######## Method 1: Use q-value ########
#From Capblancq et al. (2018) Mol. Ecol. Resources (script in supplemental)

#function for calculating p-values & q-values for all loci
rdadapt<-function(rda,K)
{
  zscores<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}

#chose best number of axes for analysis (should be 2 bc really only have 2 axes)
bestk <- ggplot() +
  geom_line(aes(x=c(1:length(rda_results$CCA$eig)), y=as.vector(rda_results$CCA$eig)), 
            linetype="dotted",size = 1.5, color="darkgrey") +
  geom_point(aes(x=c(1:length(rda_results$CCA$eig)), y=as.vector(rda_results$CCA$eig)), 
             size = 3, color="darkgrey") +
  scale_x_discrete(name = "Ordination axes", limits=c(1:9)) + ylab("Inertia") +
  theme_bw()
bestk

#calculate p-values & q-values
res_rdadapt<-rdadapt(rda_results, 2) #use # axes ID-ed in previous plot

#print list of outlier loci (q.values < 0.1) --> 10% or less of ID-d outliers could be false positives
RDA_outliers <- c(which(res_rdadapt[,2] < 0.1)) #664 loci (really 332, each SNP reported 2X)
write.table(RDA_outliers, "Data/RDA_outliers_TCapblancqmethod_list.txt")

######## Method 2: Use loading scores & p-values ########
#From Forester et al. (2016) Mol. Ecol.
#use loci with scores +/- 3 SD from mean score for that axis to ID outlier loci (2-tailed p-value = 0.0027)
#these are loci that are more strongly associated with the axis (and thus one (or more) of the env variables)
#those with sig relationship (p < 0.001) when regressed against env. variable are considered outlier

#loading scores for each locus
loci.scr <- scores(rda_results, display = "species", scaling = 3, choices = c(1,2)) #species bc vegan package (typically looking at community data)

#### Calculate means & SD ####
mean.rda <- colMeans(loci.scr) #calculate mean per RDA

sd.rda1 <- sd(loci.scr[,"RDA1"]) #calculate SD per RDA
sd.rda2 <- sd(loci.scr[,"RDA2"])

#### Pull outliers ####
#set < 3 SD cut-off
rda1.hi <- mean.rda[1] + 3*sd.rda1
rda2.hi <- mean.rda[2] + 3*sd.rda2

#set > 3 SD cut-off
rda1.lo <- mean.rda[1] - 3*sd.rda1
rda2.lo <- mean.rda[2] - 3*sd.rda2

#pull loci with score < or > 3 SDs
rda.cans <- c((which(loci.scr[,"RDA1"] > rda1.hi)),
              (which(loci.scr[,"RDA1"] < rda1.lo)),
              (which(loci.scr[,"RDA2"] > rda2.hi)),
              (which(loci.scr[,"RDA2"] < rda2.lo)))
rda.cans <- unique(sort(rda.cans)) #remove duplicates, left with 242 candidate loci (really 121, each SNP reported 2X)

#### Check associations with env variables ####
#linear regression between allele frequencies and environmental variables for the RDA outliers

#### Run LMs ####
#make list to populate with LM info for each env variable
SSSMean_lms <- list()
SSTMean_lms <- list()
SSTMin_lms <- list()
SSTMax_lms <- list()
Lat_lms <- list()

#for loop to run LM for each cand loci/env variable combo
for (i in 1:length(rda.cans)){
  SSSMean_lms[[i]] <- lm(mac2_center[,rda.cans[i]] ~ env_scaled[,"SSS_Mean"])
  SSTMean_lms[[i]] <- lm(mac2_center[,rda.cans[i]] ~ env_scaled[,"SST_Mean"])
  SSTMin_lms[[i]] <- lm(mac2_center[,rda.cans[i]] ~ env_scaled[,"SST_Min"])
  SSTMax_lms[[i]] <- lm(mac2_center[,rda.cans[i]] ~ env_scaled[,"SST_Max"])
  Lat_lms[[i]] <- lm(mac2_center[,rda.cans[i]] ~ env_scaled[,"Latitude"])
}

#### Calculate p-values ####
#make list to populate with p-values for each env variable
SSSMean_pval <- list()
SSTMean_pval <- list()
SSTMin_pval <- list()
SSTMax_pval <- list()
Lat_pval <- list()

#function to calculate p-values
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

#for loop to run p-value for each cand loci/env variable combo
for (i in 1:length(SSSMean_lms)){
  SSSMean_pval[[i]] <- lmp(SSSMean_lms[[i]])
  SSTMean_pval[[i]] <- lmp(SSTMean_lms[[i]])
  SSTMin_pval[[i]] <- lmp(SSTMin_lms[[i]])
  SSTMax_pval[[i]] <- lmp(SSTMax_lms[[i]])
  Lat_pval[[i]] <- lmp(Lat_lms[[i]])
}

#identify which loci have significant p-value for each env variable
SSSMean_cands <- which(SSSMean_pval < 0.0001)
length(SSSMean_cands)#42 (21) cands
SSTMean_cands <- which(SSTMean_pval < 0.0001)
length(SSTMean_cands) #118 (59) cands
SSTMin_cands <- which(SSTMin_pval < 0.0001)
length(SSTMin_cands) #120 (60) cands
SSTMax_cands <- which(SSTMax_pval < 0.0001)
length(SSTMax_cands) #112 (56) cands
Lat_cands <- which(Lat_pval < 0.0001)
length(Lat_cands) #106 (53) cands

#combine p-values and write out
SSSMean_pval <- unlist(SSSMean_pval)
SSTMean_pval <- unlist(SSTMean_pval)
SSTMin_pval <- unlist(SSTMin_pval)
SSTMax_pval <- unlist(SSTMax_pval)
Lat_pval <- unlist(Lat_pval)

pval_tot <- cbind(rda.cans, SSSMean_pval, SSTMean_pval, SSTMin_pval, SSTMax_pval, Lat_pval)

#write out df of p-values
write.csv(pval_tot, "Data/pval_RDA.csv")

#### Pull correlations between cand loci/env variables ####
#write function to pull loci & loadings +/- 3 SDs (same as above, just pulls more info)
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

#pull loci with score < or > 3 SDs again
cand1 <- outliers(loci.scr[,1],3) #158 (79)
cand2 <- outliers(loci.scr[,2],3) #84 (42)

ncand <- length(cand1) + length(cand2)
ncand

#take cand lists and put into a dataframe
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))

#combine cand dataframes into one
colnames(cand1) <- colnames(cand2) <- c("axis","snp","loading")

cand <- rbind(cand1, cand2)
  cand$snp <- as.character(cand$snp)

#pull correlations btwn cands and env variables
#make matrix to populate with correlations for each env variable
corr_mat <- matrix(nrow=(ncand), ncol=5)  # 5 columns for 5 predictors
colnames(corr_mat) <- c("SSSMean", "SSTMean", "SSTMin", "SSTMax", "Lat")

#pull correlations
for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- mac2_center[,nam]
  corr_mat[i,] <- apply(env_scaled,2,function(x) cor(x,snp.gen))
}

#combine candidate loci names & scores with correlation with each env variable
cand <- cbind.data.frame(cand,corr_mat)  
head(cand)

#add RDA_id to each row
cand_RDAid <- cbind(rda.cans, cand)

#write out df of associations
write.csv(cand_RDAid, "Data/RDA_canloci_associations.csv")

################################################################################################################################################

######## Plot RDA ########

#create table with bi-plot scores for SST_Min, SST_Max & Latitude
SST_Min_scores <- c(1, 0.00859) #(RDA1, RDA2) for N4 included 
SST_Max_scores <- c(0.9719, 0.2355)
Lat_scores <- c(-0.9698, 0.244)
#SST_Min_scores <- c(0.9967, 0.08085) #(RDA1, RDA2) for N4 not included 
#SST_Max_scores <- c(0.9903, -0.1388)
#Lat_scores <- c(-0.9455, -0.3256)

#merge together
env_scores <- as.data.frame(cbind(SST_Min_scores, SST_Max_scores, Lat_scores))
env_scores <- as.data.frame(t(env_scores)) #transpose so RDA scores in columns, env variables in rows
colnames(env_scores) <- c("RDA1", "RDA2")

#plot RDA
plot(rda_results, choices = c(1,2), scaling=3) #only 2 axes (2 RDAs bc only 2 explanatory variables), scaling=3 (symmetrical scaling) scales each SNP and individ by square root of eigenvalues
arrows(0, 0, env_scores[,1], env_scores[,2], length = 0, lty = 1, col="blue") #adding other explanatory variables