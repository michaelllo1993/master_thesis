
options(warn = -1)
# Load libraries ----------------------------------------------------------

library(Biostrings)

# Get the directory -------------------------------------------------------
wd = getwd()  

# Get the argument --------------------------------------------------------
args = commandArgs(trailingOnly = TRUE)
OoI = args[1]
# Read the protein sequences ----------------------------------------------
setwd(paste(wd, "/data/proteinSequences/", sep = ""))
myFile <- paste(OoI,".csv",sep = "")
proteins = read.csv(myFile, stringsAsFactors = F)

setwd(paste(wd,"/data/readData",sep = ""))
protein = list()
protein[[1]] = proteins
names(protein) = OoI
saveRDS(protein, file = paste("readData_proteins_",OoI,".rds",sep = ""))

# Read the cDNA sequences -------------------------------------------------

setwd(paste(wd, "/data/cDNAsequences/", sep = ""))
myFile <- paste(OoI,"_cDNA.csv",sep = "")
command = paste("grep -Ev '\\|NA\\|' ", myFile, " > ", "done.csv", sep = "")
system(command,wait = T)
system(paste("mv done.csv ",myFile,sep = ""),wait = T)
cDNAs = read.csv(myFile, stringsAsFactors = F)
cDNAs = cDNAs[!(cDNAs[,1])=="",]

########### Split START and STOP values ##########

cDNAs$START = lapply(cDNAs$START, function(x)
  sort(as.numeric(strsplit(x, ";")[[1]])))
#split by semicolon
cDNAs$STOP = lapply(cDNAs$STOP, function(x)
  sort(as.numeric(strsplit(x, ";")[[1]])))
#split by semicolon

cDNA = list()
cDNA[[1]] = cDNAs
names(cDNA) = OoI
setwd(paste(wd,"/data/readData",sep = ""))
saveRDS(cDNA, file=paste("readData_cDNA_",names(cDNA),".rds",sep = ""))
