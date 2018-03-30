

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
  cDNA_names[k] = strsplit(myFiles2, "[.]")[[k]][1]
  command = paste("grep -Ev '\\|NA\\|' ", myFiles2[k], " > ", "done.csv", sep = "")
  print(command)
  system(command,wait = T)
  system(paste("mv done.csv ",myFiles2[k],sep = ""),wait = T)
  cDNA[[k]] = read.csv(myFiles2[k], stringsAsFactors = F)
  cDNA[[k]] = cDNA[[k]][!(cDNA[[k]][,1])=="",]
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
for(org in seq_len(length(myFiles2))){
  saveRDS(cDNA[[org]], file =paste("readData_cDNA_",names(cDNA)[org],".rds",sep = ""))
}
