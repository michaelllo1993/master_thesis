

# Load libraries ----------------------------------------------------------

library(Biostrings)

# Get the directory and arguments -------------------------------------------------------
wd = getwd()

# Read the predicted signal peptide cleavage sites ------------------------

setwd(paste(wd, "/results/SP", sep = ""))
myFiles3 <- list.files(pattern = "*signalp_positives.out")
SP = list()
SP_names = c()
for (k in 1:length(myFiles3)) {
  SP_names[k] = strsplit(myFiles3, "[.]")[[k]][1]
  SP[[k]] = read.table(myFiles3[k], stringsAsFactors = F)
  SP[[k]] = SP[[k]][,c(1,2,5)]
}

names(SP) = SP_names

setwd(paste(wd,"/data/readData",sep = ""))
saveRDS(SP, file = "readData_SP.rds")

# Read SAAR data within signal peptides -----------------------------------

setwd(paste(wd, "/results/SP", sep = ""))
myFiles4 <- list.files(pattern = "extracted_sigp_*")
SAAR = list()
SAAR_names = c()
for (k in 1:length(myFiles4)) {
  SAAR_names[k] = strsplit(myFiles4, "[.]")[[k]][1]
  SAAR[[k]] = read.csv(myFiles4[k], stringsAsFactors = F,header = F)
  SAAR[[k]] = SAAR[[k]][,-c(3,4)]
  faultyRows=which(!grepl('ENS',SAAR[[k]][,1]))
  if(length(faultyRows)>0){
    SAAR[[k]] = SAAR[[k]][-c(which(!grepl('ENS',SAAR[[k]][,1]))),]
  }
}
SAAR_names[k] = strsplit(myFiles4, "[.]")[[k]][1]
names(SAAR) = SAAR_names

setwd(paste(wd,"/data/readData",sep = ""))
saveRDS(SAAR, file = "readData_SAAR.rds")
