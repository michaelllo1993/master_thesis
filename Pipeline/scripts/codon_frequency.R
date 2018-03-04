
# Load libraries ----------------------------------------------------------

library(Biostrings)
library(seqinr)

# Load data ---------------------------------------------------------------

args = commandArgs(trailingOnly = TRUE)

wd = getwd()
# setwd(paste(wd, "/data", sep = ""))
# load("readData.RData")

leucine_codons = names(which(GENETIC_CODE=="L"));
nucleotides = c("A","T","G","C")

# Processing the input mapper file ----------------------------------------

# Filtering only for rows with all available orthologues
original_mapper_name = args[1]
organisms_names = args[-1]
common_mapper_name = paste(strsplit(original_mapper_name,"\\.")[[1]][1],"_common.",strsplit(original_mapper_name,"\\.")[[1]][2],sep = "")
command = paste("grep -Ev 'NULL' ",original_mapper_name," > ",common_mapper_name,sep = "")
system(command)

orthologues = read.csv(common_mapper_name,header = F)[,-1]
colnames(orthologues) = organisms_names


  























