library(biomaRt)
library(seqinr)
wd = getwd()

args = c("homo_sapiens",       "pan_troglodytes"  ,  "gorilla_gorilla",    "macaca_mulatta"   ,  "mus_musculus"       ,"rattus_norvegicus", "bos_taurus"        , "gallus_gallus"   ,   "xenopus_tropicalis")

organisms = args
organisms = as.vector(sapply(organisms, function(x) paste(s2c(strsplit(x,split = "_")[[1]][1])[1],strsplit(x,split = "_")[[1]][2],sep = "")))


for(organism in seq_len(length(organisms))){
  mart_name <- useMart("ensembl", paste(organisms[organism],"_gene_ensembl",sep = ""))
  protein_sequences[[organism]] = as.matrix(getBM(
    attributes = "peptide",
    mart = mart_name
  ))
}
