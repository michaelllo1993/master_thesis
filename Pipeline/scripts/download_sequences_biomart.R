library(biomaRt)
library(seqinr)
library(dplyr)

wd = getwd()
args = commandArgs(trailingOnly = TRUE)
# args = c("homo_sapiens",       "pan_troglodytes"  ,  "gorilla_gorilla",    "macaca_mulatta"   ,  "mus_musculus"       ,"rattus_norvegicus", "bos_taurus"        , "gallus_gallus"   ,   "xenopus_tropicalis")
orgs=args
organisms = orgs
organisms = as.vector(sapply(organisms, function(x) paste(s2c(strsplit(x,split = "_")[[1]][1])[1],strsplit(x,split = "_")[[1]][2],sep = "")))

# read cDNA sequences ----------------------------------------------------------

cDNA=list()
for(organism in seq_len(length(organisms))){
  mart <- useEnsembl(biomart = "ensembl",dataset = paste(organisms[organism],"_gene_ensembl",sep = ""))
  # vector of all gene IDs (only used to improve processing later)
  genes <- getBM(attributes = "ensembl_gene_id", mart = mart)
  # first table with CDNA start/end values
  transcripts <- getBM(values = genes$ensembl_gene_id, filters = "ensembl_gene_id", attributes = c("ensembl_transcript_id","ensembl_peptide_id","cdna_coding_start","cdna_coding_end"), mart = mart)
  # process the dataframe: collapse entries with the same transcript ID
  transcripts_mod = transcripts %>% group_by(ensembl_transcript_id) %>% summarise(START = paste(cdna_coding_start,collapse = ","), STOP = paste(cdna_coding_end,collapse = ","))
  rm(transcripts)
  # second table with sequences
  transcripts2 <- getBM(values = genes$ensembl_gene_id, filters = "ensembl_gene_id", attributes = c("ensembl_transcript_id","cdna"),mart = mart)
  rm(mart,genes)
  # merge the data.frames
  merged = inner_join(transcripts2,transcripts_mod)
  rm(transcripts2,transcripts_mod)
  # assign to the output var
  cDNA[[organism]] = merged
  rm(merged)
}
names(cDNA) = orgs
saveRDS(cDNA,paste(wd,"data/readData/readData_cDNA.rds",sep = ""))


# read protein sequences -------------------------------------------------------

protein=list()
for(organism in seq_len(length(organisms))){
  mart <- useEnsembl(biomart = "ensembl",dataset = paste(organisms[organism],"_gene_ensembl",sep = ""))
  # vector of all gene IDs (only used to improve processing later)
  genes <- getBM(attributes = "ensembl_gene_id", mart = mart)
  # first table with CDNA start/end values
  peptides <- getBM(values = genes$ensembl_gene_id, filters = "ensembl_gene_id", attributes = c("ensembl_peptide_id","peptide"), mart = mart)
  # assign to the output var 
  protein[[organism]] = peptides
  rm(peptides,genes,mart)
}
names(protein) = orgs

saveRDS(protein,paste(wd,"data/readData/readData_proteins.rds",sep = ""))
