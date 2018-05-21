# args = c("data/proteinSequences/homo_sapiens.txt","results/SP/homo_sapiens_signalp_positives.out","results/SP/extracted_sigp_homo_sapiens.out")
# setwd("~/Uczelnia/MGR/master_thesis/Pipeline/")

wd = getwd()
args = commandArgs(trailingOnly = TRUE)
args = paste(wd,args,sep = "/")

filename = strsplit(args[1],split = "/")[[1]][length(strsplit(args[1],split = "/")[[1]])]
org_name = substr(filename,start = 1,stop = nchar(filename) - 4)
protein_count = (strsplit(system(paste("grep '>'",args[1],"| wc -l"),intern = T),split = "/")[[1]][1])
SP_count = (strsplit(system(paste("wc -l",args[2]),intern = T),split = "/")[[1]][1])
SAAR_count = (strsplit(system(paste("wc -l",args[3]),intern = T),split = "/")[[1]][1])

output_mtx = matrix(data = c("Number of proteins",protein_count,"Number of signal peptides",SP_count,"Number of SAARs",SAAR_count),nrow = 3,ncol = 2,byrow = T)
write.csv(output_mtx,file = paste(wd,"results",paste(org_name,"_preliminary_data",".csv",sep = ""),sep = "/"),row.names = F)
