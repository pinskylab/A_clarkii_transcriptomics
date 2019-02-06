#### Creating Simulation Env Data Sets for Sim Runs of BayPass ####

remove(list = ls())

#Read in environmental data
library(readr)
transcript_individs <- read_csv("Individual_Env_Data.csv", col_names = TRUE) #all 25 fish with environ data divided into 3 pops

#Subset to 5 environmental variables of interest
envir <- transcript_individs[,c('pop', 'sss_mean', 'sst_mean', 'sst_min', 'sst_max', 'lat')]

#For loop to randomize, standardize, and write 10 simulation environmental matrices
for (i in 1:10){
  envir2 <- apply(envir[,-1],2, sample) #takes random entry from each column (except col 1) and creates new matrix
  by <- list(transcript_individs$pop) #create list of population of each individual
  pop.avg <- aggregate(envir2, by = by, FUN = mean) #average each new "population" (first 8 individuals in pop 1, etc) 
  stan.pop.avg <- scale(pop.avg[,2:6]) #standardize each new environmental variable so variance is = (can skip this step if standardizing in BayPass)
  rownames(stan.pop.avg) <- c("Pop1", "Pop2", "Pop3")
  t.stan.pop.avg <- t(stan.pop.avg) #transpose matrix
  write.table(t.stan.pop.avg, file = "randomized_env_file.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
}
