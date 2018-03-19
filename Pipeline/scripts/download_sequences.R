library(biomaRt)
library(seqinr)
wd = getwd()

args = c("homo_sapiens",       "pan_troglodytes"  ,  "gorilla_gorilla",    "macaca_mulatta"   ,  "mus_musculus"       ,"rattus_norvegicus", "bos_taurus"        , "gallus_gallus"   ,   "xenopus_tropicalis")

organisms = args
organisms = as.vector(sapply(organisms, function(x) paste(s2c(strsplit(x,split = "_")[[1]][1])[1],strsplit(x,split = "_")[[1]][2],sep = "")))

protein=list()
for(organism in seq_len(length(organisms))){
  mart_name <- useEnsembl(biomart = "ensembl",dataset = paste(organisms[organism],"_gene_ensembl",sep = ""),version = 91)
  protein[[organism]] = as.matrix(getBM(
    attributes = c("ensembl_peptide_id","peptide"),
    mart = mart_name
  ))
}

cDNA=list()
for(organism in seq_len(length(organisms))){
  mart_name <- useMart("ensembl", paste(organisms[organism],"_gene_ensembl",sep = ""))
  cDNA[[organism]] = (getBM(
    attributes = c("ensembl_transcript_id","cdna_coding_start","cdna_coding_end","cdna"),
    mart = mart_name
  ))
}

mart_name <- useEnsembl(biomart = "ensembl",dataset = paste("ggallus","_gene_ensembl",sep = ""),mirror = "useast",verbose = T)
mart_name <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "ensembl", paste(OoI,"_gene_ensembl",sep = ""),host = "www.ensembl.org")



cDNA[[1]] = (getBM(
  attributes = c("ensembl_transcript_id","cdna_coding_start","cdna_coding_end","cdna"),
  mart = mart_name
))
