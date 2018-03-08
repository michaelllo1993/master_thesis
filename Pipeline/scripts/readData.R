

# Load libraries ----------------------------------------------------------

library(Biostrings)

# Get the directory -------------------------------------------------------
wd = getwd()  

# Read the protein sequences ----------------------------------------------

setwd(paste(wd, "/data/proteinSequences/", sep = ""))
myFiles <- list.files(pattern = "*.csv")
proteins = list()
prot_names = c()
for (k in 1:length(myFiles)) {
  prot_names[k] = strsplit(myFiles, "[.]")[[k]][1]
  proteins[[k]] = read.csv(myFiles[k], stringsAsFactors = F)
}
prot_names[k] = strsplit(myFiles, "[.]")[[k]][1]
names(proteins) = prot_names

setwd(paste(wd,"/data/readData",sep = ""))
saveRDS(proteins, file = "readData_proteins.rds")

# Read the mappers --------------------------------------------------------

setwd(paste(wd, "/data/mappers/", sep = ""))
myFiles1 <- list.files(pattern = "*.csv")
mapper = list()
map_names = c()
for (k in 1:length(myFiles1)) {
  map_names[k] = strsplit(myFiles1, "[.]")[[k]][1]
  mapper[[k]] = read.csv(myFiles1[k], stringsAsFactors = F)
}
names(mapper) = map_names

setwd(paste(wd,"/data/readData",sep = ""))
saveRDS(mapper, file = "readData_mapper.rds")

# Read the cDNA sequences -------------------------------------------------

setwd(paste(wd, "/data/cDNAsequences/", sep = ""))
myFiles2 <- list.files(pattern = "*.csv")
cDNA = list()
cDNA_names = c()
for (k in 1:length(myFiles2)) {
  cDNA_names[k] = strsplit(myFiles, "[.]")[[k]][1]
  cDNA[[k]] = read.csv(myFiles2[k], stringsAsFactors = F)
}
names(cDNA) = sapply(strsplit(cDNA_names,"_"), function(x) paste(x[1],x[2],sep = "_"))

########### Split START and STOP values ##########

for (org in seq(1, length(cDNA), by = 1)) {
  cDNA[[org]]$START = lapply(cDNA[[org]]$START, function(x)
    sort(as.numeric(strsplit(x, ";")[[1]])))
  #split by semicolon
  cDNA[[org]]$STOP = lapply(cDNA[[org]]$STOP, function(x)
    sort(as.numeric(strsplit(x, ";")[[1]])))
  #split by semicolon
}

setwd(paste(wd,"/data/readData",sep = ""))
saveRDS(cDNA, file = "readData_cDNA.rds")

# Read the predicted signal peptide cleavage sites ------------------------

setwd(wd)
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

setwd(wd)
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
