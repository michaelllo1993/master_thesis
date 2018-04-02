library(biomaRt)
library(seqinr)
library(dplyr)

wd = getwd()
args = commandArgs(trailingOnly = TRUE)
organism_full = args
organisms = paste(s2c(strsplit(organism_full, split = "_")[[1]][1])[1], strsplit(organism_full, split = "_")[[1]][2], sep = "")

# read cDNA sequences ----------------------------------------------------------

cDNA = list()
# mart <- useEnsembl(biomart = "ensembl",dataset = paste(organisms,"_gene_ensembl",sep = ""),host = "useast.ensembl.org")
mart <-
  useEnsembl(biomart = "ensembl",
             dataset = paste(organisms, "_gene_ensembl", sep = ""))
# vector of all gene IDs (only used to improve processing later)
genes <- getBM(attributes = "ensembl_gene_id", mart = mart)
# first table with CDNA start/end values
tryCatch(
  transcripts <-
    getBM(
      values = genes$ensembl_gene_id,
      filters = "ensembl_gene_id",
      attributes = c(
        "ensembl_transcript_id",
        "ensembl_peptide_id",
        "cdna_coding_start",
        "cdna_coding_end"
      ),
      mart = mart
    ),
  error = function(error_message) {
    message(
      paste(
        "The downloading of the ",
        organism_full,
        " cDNA sequences has failed due to server problems. Please run the rule once again.",
        sep = ""
      )
    )
    return(NA)
  }
)

# process the dataframe: collapse entries with the same transcript ID
transcripts_mod = transcripts %>% group_by(ensembl_transcript_id) %>% summarise(
  ensembl_peptide_id = unique(ensembl_peptide_id)[1],
  START = paste(cdna_coding_start, collapse = ","),
  STOP = paste(cdna_coding_end, collapse = ",")
)
rm(transcripts)
# second table with sequences
tryCatch(
  transcripts2 <-
    getBM(
      values = genes$ensembl_gene_id,
      filters = "ensembl_gene_id",
      attributes = c("ensembl_transcript_id", "cdna"),
      mart = mart
    ),
  error = function(error_message) {
    message(
      paste(
        "The downloading of the ",
        organism_full,
        " cDNA sequences has failed due to server problems. Please run the rule once again.",
        sep = ""
      )
    )
    return(NA)
  }
)
rm(mart, genes)
# merge the data.frames
merged = inner_join(transcripts2, transcripts_mod)
#delete the sequences that do not have the peptide counterpart
merged = merged[-which(merged$ensembl_peptide_id == ""), ]
#rearrange columns order
merged = merged[, c(1, 3, 4, 5, 2)]
names_merged = names(merged)
headers = c()
for (i in seq_len(nrow(merged))) {
  starts = gsub(pattern = "NA,", replacement = "", merged$START[i])
  merged$START[i] = starts
  stops = gsub(pattern = "NA,", replacement = "", merged$STOP[i])
  merged$STOP[i] = stops
  headers[i] = paste(
    merged$ensembl_transcript_id[i],
    merged$ensembl_peptide_id[i],
    gsub(pattern = ",", replacement = ";", starts),
    gsub(pattern = ",", replacement = ";", stops),
    sep = "|"
  )
}
rm(transcripts2, transcripts_mod)
# assign to the output var
cDNA[[1]] = merged
write.fasta(
  as.list(merged$cdna),
  names = headers,
  file.out = paste(wd, "/data/cDNAsequences/", organism_full, "_cDNA.txt", sep = "")
)
#Convert the .txr file to .csv format with a Perl script
command = (paste(
  "perl scripts/seq2csv_cDNA.pl",
  paste(wd, "/data/cDNAsequences/", organism_full, "_cDNA.txt", sep = ""),
  sep = " "
))
system(command)
#remove rows with 'NA'
# command = paste("grep -Ev '\\|NA\\|'", paste(wd,"/data/cDNAsequences/",organism_full,"_cDNA.csv",sep = ""), " > ", paste(wd,"/data/cDNAsequences/",organism_full,"_cDNA_edited.csv",sep = ""), sep = "")
# system(command)
# system("mv ",paste(paste(wd,"/data/cDNAsequences/",organism_full,"_cDNA_edited.csv",sep = ""),paste(wd,"/data/cDNAsequences/",organism_full,"_cDNA.csv",sep = ""),sep=""))
rm(cDNA)
#read the .csv file
cDNA = read.csv(
  paste(wd, "/data/cDNAsequences/", organism_full, "_cDNA.csv", sep = ""),
  stringsAsFactors = F
)
cDNA = cDNA[!(cDNA[, 1]) == "", ]
cDNA$START = lapply(cDNA$START, function(x)
  sort(as.numeric(strsplit(x, ";")[[1]])))
#split by semicolon
cDNA$STOP = lapply(cDNA$STOP, function(x)
  sort(as.numeric(strsplit(x, ";")[[1]])))
#split by semicolon
names(cDNA) = names_merged
cDNA = list(cDNA)
names(cDNA) = organism_full
#save for later use
saveRDS(cDNA,
        paste(
          wd,
          "/data/readData/readData_cDNA_",
          organism_full,
          ".rds",
          sep = ""
        ))


# read protein sequences -------------------------------------------------------

protein = list()
mart <-
  useEnsembl(biomart = "ensembl",
             dataset = paste(organisms, "_gene_ensembl", sep = ""))
# vector of all gene IDs (only used to improve processing later)
genes <- getBM(attributes = "ensembl_gene_id", mart = mart)
# first table with CDNA start/end values
tryCatch(
  peptides <-
    getBM(
      values = genes$ensembl_gene_id,
      filters = "ensembl_gene_id",
      attributes = c("ensembl_peptide_id", "ensembl_transcript_id", "peptide"),
      mart = mart
    ),
  error = function(error_message) {
    message(
      paste(
        "The downloading of the ",
        organism_full,
        " protein sequences has failed due to server problems. Please run the rule once again.",
        sep = ""
      )
    )
    return(NA)
  }
)
headers = c()
for (i in seq_len(nrow(peptides))) {
  headers[i] = paste(peptides$ensembl_peptide_id[i],
                     peptides$ensembl_transcript_id[i],
                     sep = "|")
}
write.fasta(
  sequences = as.list(peptides$peptide),
  names = headers,
  as.string = T,
  file.out = paste(wd, "/data/proteinSequences/", organism_full, ".txt", sep = "")
)
system(paste("sed -i '/Sequence\ unavailable/d' ",paste(wd, "/data/proteinSequences/", organism_full, ".txt", sep = ""),sep = ""))
# assign to the output var
protein[[1]] = peptides
rm(peptides, genes, mart)
names(protein) = organism_full

saveRDS(
  protein,
  paste(
    wd,
    "/data/readData/readData_proteins_",
    organism_full,
    ".rds",
    sep = ""
  )
)
