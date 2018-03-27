library(biomaRt)
library(seqinr)
library(dplyr)

wd = getwd()
args = commandArgs(trailingOnly = TRUE)
organism_full=args
organisms = paste(s2c(strsplit(organism_full,split = "_")[[1]][1])[1],strsplit(organism_full,split = "_")[[1]][2],sep = "")

# read cDNA sequences ----------------------------------------------------------

cDNA=list()
mart <- useEnsembl(biomart = "ensembl",dataset = paste(organisms,"_gene_ensembl",sep = ""),host = "useast.ensembl.org")
# vector of all gene IDs (only used to improve processing later)
genes <- getBM(attributes = "ensembl_gene_id", mart = mart)
# first table with CDNA start/end values
transcripts <- getBM(values = genes$ensembl_gene_id, filters = "ensembl_gene_id", attributes = c("ensembl_transcript_id","ensembl_peptide_id","cdna_coding_start","cdna_coding_end"), mart = mart)
# process the dataframe: collapse entries with the same transcript ID
transcripts_mod = transcripts %>% group_by(ensembl_transcript_id) %>% summarise(ensembl_peptide_id = unique(ensembl_peptide_id)[1],START = paste(cdna_coding_start,collapse = ","), STOP = paste(cdna_coding_end,collapse = ","))
rm(transcripts)
# second table with sequences
transcripts2 <- getBM(values = genes$ensembl_gene_id, filters = "ensembl_gene_id", attributes = c("ensembl_transcript_id","cdna"),mart = mart)
rm(mart,genes)
# merge the data.frames
merged = inner_join(transcripts2,transcripts_mod)
headers=c()
for(i in seq_len(nrow(merged))){
  starts=gsub(pattern = "NA,",replacement = "",merged$START[i])
  stops=gsub(pattern = "NA,",replacement = "",merged$STOP[i])
  headers[i] = paste(merged$ensembl_transcript_id[i],merged$ensembl_peptide_id[i],gsub(pattern = ",",replacement = ";",starts),gsub(pattern = ",",replacement = ";",stops),sep = "|")
}
rm(transcripts2,transcripts_mod)
# assign to the output var
cDNA[[1]] = merged
write.fasta(as.list(merged$cdna),names = headers,file.out = paste(wd,"/data/cDNAsequences/",organism_full,"_cDNA.txt",sep = ""))
names(cDNA) = organism_full
saveRDS(cDNA,paste(wd,"/data/readData/readData_cDNA_",organism_full,".rds",sep = ""))


# read protein sequences -------------------------------------------------------

protein=list()
mart <- useEnsembl(biomart = "ensembl",dataset = paste(organisms,"_gene_ensembl",sep = ""))
# vector of all gene IDs (only used to improve processing later)
genes <- getBM(attributes = "ensembl_gene_id", mart = mart)
# first table with CDNA start/end values
peptides <- getBM(values = genes$ensembl_gene_id, filters = "ensembl_gene_id", attributes = c("ensembl_peptide_id","peptide"), mart = mart)
write.fasta(sequences = as.list(peptides$peptide),names = peptides$ensembl_peptide_id,as.string = T,file.out = paste(wd,"/data/proteinSequences/",organism_full,".txt",sep = ""))
# assign to the output var 
protein[[1]] = peptides
rm(peptides,genes,mart)
names(protein) = organism_full

saveRDS(protein,paste(wd,"data/readData/readData_proteins",organism_full,".rds",sep = ""))
