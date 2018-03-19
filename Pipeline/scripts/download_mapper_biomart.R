library(biomaRt)
library(seqinr)
wd = getwd()

orgs = args
#convert orginal organisms names to ones used in the database
orgs = as.vector(sapply(orgs, function(x) paste(s2c(strsplit(x,split = "_")[[1]][1])[1],strsplit(x,split = "_")[[1]][2],sep = "")))
#get the organism of interest - first one 
OoI = orgs[1]
organisms= orgs[-1]

#set the biomart to be used
mart_name <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "ensembl", paste(OoI,"_gene_ensembl",sep = ""),host = "www.ensembl.org")

mappers = list()
#get the count of full orthologues dowlonloads to be done (which contain 6 ortohlogues)
full_batches = length(organisms)%/%6
#perform org number of full downloads 
for(org in seq_len(full_batches)){
  #get proper organism for each batch
  start_element = org*6-5 
  stop_element = org*6
  #set the attributes
  atrbts = append("ensembl_peptide_id",paste(organisms[start_element:stop_element],"_homolog_ensembl_peptide",sep = ""))
  #acquire the data
  mappers[[org]] = as.matrix(getBM(
    attributes = atrbts,
    mart = mart_name
  ))
  #save the data with non-empyt rows in first column (with IDs of organism of interest)
  mappers[[org]] = mappers[[org]][which(mappers[[org]][,1] != ""),]
}

#perform a non-full donload for the organisms that left
leftover = organisms[(stop_element+1):length(organisms)]
atrbts = append("ensembl_peptide_id",paste(leftover,"_homolog_ensembl_peptide",sep = ""))
mappers[[org+1]] = as.matrix(getBM(
  attributes = atrbts,
  mart = mart_name
))
mappers[[org+1]] = mappers[[org+1]][which(mappers[[org+1]][,1] != ""),]

#set the function for mapper list preprocessing
preprocess_mapper_list <- function(x){
  x=x[which(!duplicated(x[,1])),]
  rownames(x) = x[,1]
  x=x[,-1]
  return(x)
}

#preprocess the mappers
mappers = lapply(mappers, preprocess_mapper_list)
#initialize the mapper list
mapper = mappers[[1]]
#merge all mappers together
for(i in seq_len(length(mappers)-1)){
  mapper = merge(mapper,mappers[[i+1]],by="row.names")  
}

#write mappers to the file
write.table(mapper,file = paste(wd,"/",OoI,".tsv",sep=""),sep = "\t",row.names = F)


