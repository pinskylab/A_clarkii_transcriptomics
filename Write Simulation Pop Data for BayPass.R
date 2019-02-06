#### Creating Simulation Population Data Sets for Sim Runs of BayPass ####

remove(list = ls())

#Read in environmental data
library(readr)
transcript_individs <- read_csv("Individual_Env_Data.csv", col_names = TRUE) #all 25 fish with environ data divided into 3 pops

#Randomize individuals in population
ran1_individ <- apply(transcript_individs[,1],2, sample)

#Read in csv file
output_hicov2_snps_only <- read_csv("output.hicov2.snps.only.csv", col_names = TRUE)

#Trim to individuals (filters columns with NPJ & digit in the name)
cols <- grepl('[NPJ][[:digit:]]', names(output_hicov2_snps_only)) #creates pattern of logic values (T/F if fit grepl filter)
output_hicov2_snps_only <- output_hicov2_snps_only[,cols] #keeps only columns with T value
dim(output_hicov2_snps_only) #check to make sure correct dimensions (10977 x 25)

#Rename column headers with randomized individual vector (EX: makes J11 column now P1 column to be included with Philippines pop)
names(output_hicov2_snps_only) <- c(ran1_individ)

#Find columns for each population
pcols <- grepl('^P[[:digit:]]', names(output_hicov2_snps_only)) #keeps columns with only P & digit in the name
jcols <- grepl('^J[[:digit:]]', names(output_hicov2_snps_only))
ncols <- grepl('^N[[:digit:]]', names(output_hicov2_snps_only))

#Remove extra information 
clean <- function(x) {
  substitute <- gsub(pattern = ":.*", "", x) #removes everything after : and replaces with nothing
}

jgeno <- apply(output_hicov2_snps_only[,jcols], MARGIN = 1, FUN = clean) #remove extra information for japanese individuals
pgeno <- apply(output_hicov2_snps_only[,pcols], MARGIN = 1, FUN = clean) #remove extra information for philippines individuals
ngeno <- apply(output_hicov2_snps_only[,ncols], MARGIN = 1, FUN = clean) #remove extra information for indonesia individuals

#Transform data frames
trans_jgeno <- t(jgeno) #flip rows and columns so individuals in column and loci in rows
trans_pgeno <- t(pgeno)
trans_ngeno <- t(ngeno)

#Count 0s and 1s by population
count01s <- function(x){
  als <- unlist(strsplit(split='/', fixed=TRUE, as.character(x))) #get alleles (splitting at / and representing objects as characters (instead of numerics))
  out <- c(sum(als==0), sum(als==1)) #sum all 0s and 1s
  return(out)
}

jals <- apply(trans_jgeno, MARGIN = 1, FUN = count01s) #count japanese alleles (applying over rows -- summing 0s and 1s in each row individually)
pals <- apply(trans_pgeno, MARGIN = 1, FUN = count01s) #count philippines alleles (applying over rows)
nals <- apply(trans_ngeno, MARGIN = 1, FUN = count01s) #count indonesia alleles (applying over rows)

#Transform to 3-column baypass format for output
#Uses 1:n indexing to turn allele counts from 2xnloci matrices to nloci*2vectors
be <- data.frame(J=jals[1:length(jals)], P=pals[1:length(pals)], N=nals[1:length(nals)])
dim(be) #check to make sure correct dimensions (1600 x 3)
summary(be)

#Write out
write.table(be, file = "sim10.baypass.format.txt", sep='\t', col.names=TRUE, row.names=FALSE) #table with J P N columns, 2*nloci rows
write.table(ran1_individ, file = "sim10.popassignments.txt", sep='\t', col.names=TRUE, row.names=FALSE) #list of sim pop assignments for individuals