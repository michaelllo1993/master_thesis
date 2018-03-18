library(biomaRt)
library(seqinr)

orgs = args
orgs = as.vector(sapply(orgs, function(x) paste(s2c(strsplit(x,split = "_")[[1]][1])[1],strsplit(x,split = "_")[[1]][2],sep = "")))
OoI = orgs[1]
organisms= orgs[-1]

mart_name <- useMart("ensembl", paste(OoI,"_gene_ensembl",sep = ""))

mappers = list()

for(org in seq_len(full_batches)){
  start_element = org*6-5 
  stop_element = org*6
  atrbts = append("ensembl_peptide_id",paste(organisms[start_element:stop_element],"_homolog_ensembl_peptide",sep = ""))
  print(atrbts)
  mappers[[org]] = as.matrix(getBM(
    attributes = atrbts,
    mart = mart_name
  ))
  mappers[[org]] = mappers[[org]][which(mappers[[org]][,1] != ""),]
}


leftover = organisms[(stop_element+1):length(organisms)]
atrbts = append("ensembl_peptide_id",paste(leftover,"_homolog_ensembl_peptide",sep = ""))
print(atrbts)
mappers[[org+1]] = as.matrix(getBM(
  attributes = atrbts,
  mart = mart_name
))
mappers[[org+1]] = mappers[[org+1]][which(mappers[[org+1]][,1] != ""),]

