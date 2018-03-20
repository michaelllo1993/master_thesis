



# Load libraries ----------------------------------------------------------

library(Biostrings)
library(seqinr)

# Load data ---------------------------------------------------------------

args = commandArgs(trailingOnly = TRUE)

wd = getwd()
setwd(paste(wd, "/data/readData", sep = ""))
cDNA = readRDS("readData_cDNA.RData")
SP = readRDS("readData_SP.RData")
SAAR = readRDS("readData_SAAR.RData")
setwd("..")

# Processing the input mapper file ----------------------------------------

# Filtering only for rows with all available orthologues
original_mapper_name = args[1]
repeat_unit = args[2]
organisms_names = args[c(-1, -2)]
OoI = organisms_names[1]
unit_codons = names(which(GENETIC_CODE == repeat_unit))

nucleotides = c("A", "T", "G", "C")
unique_codon_names = c()

# Build the corrected  mapper file name
common_mapper_name = paste(
  strsplit(original_mapper_name, "\\.")[[1]][1],
  "_common.",
  strsplit(original_mapper_name, "\\.")[[1]][2],
  sep = ""
)
# Grep only the lines without 'NULL' in them - ones where all orthologues are available
command = paste("grep -Ev 'NULL' ",
                original_mapper_name,
                " > ",
                common_mapper_name,
                sep = "")
system(command)
# Read corrected file
orthologues = read.csv(common_mapper_name, header = F,skip = 1)[, -1]
colnames(orthologues) = organisms_names

# Very important - order the cDNA list elements as they are in the orthologues mtx
cDNA = cDNA[organisms_names]
SP = SP[paste(organisms_names,"_signalp_positives",sep = "")]
SAAR = SAAR[paste("extracted_sigp_",organisms_names,sep = "")]

# The analysis overall ------------------------------------------------------------
# 
# for (org in seq_len(ncol(orthologues))) {
#   result = c()
#   
#   # Prepare empty result matrix
#   result = matrix("", nrow(orthologues), length(unit_codons) + 1)
#   
#   for (i in seq_len(nrow(orthologues))) {
#     # Find cDNA of interest
#     which_cDNA = which(cDNA[[org]]$ID == orthologues[i, org])
#     if (length(which_cDNA) > 0) {
#       # Get the sequence
#       cdna = s2c(cDNA[[org]]$SEQUENCE[which_cDNA])
#       # Get the starting and stopping positions of the translated sequence
#       codon_start = c()
#       codon_start = cDNA[[org]]$START[which_cDNA][[1]]
#       codon_stop = c()
#       codon_stop = cDNA[[org]]$STOP[which_cDNA][[1]]
#       vec = c()
#       # Concatenate all translated region
#       for (j in seq_len(length(codon_start))) {
#         vec = append(x = vec,
#                      values = seq(codon_start[j], codon_stop[j], by = 1))
#         
#       }
#       # Get only the translated nucleotides
#       translated = cdna[vec]
#       
#       # In case of faulty ("N") nucleotide - replace it with a random one
#       translated[which(translated == "N")] = nucleotides[sample(1:4, 1)]
#       # Sequence to codons and get the frequencies
#       codons_freq = table(as.character(codons(DNAString(
#         c2s(translated)
#       ))))
#       
#       # Just unit of interest codons codons
#       avail_unit_codons_freq = codons_freq[unit_codons]
#       
#       # Combine results into a matrix
#       result[i, ] = append(cDNA[[org]]$cDNA_ID[which_cDNA], as.vector(avail_unit_codons_freq))
#     }
#   }
#   colnames(result) = append("cDNA_ID",unit_codons)
#   # Write the result to .csv file
#   write.csv(
#     result,
#     file = paste(
#       sep = "",
#       getwd(),
#       "/codon_frequency/",
#       OoI,
#       "_orthologues/",
#       organisms_names[org],
#       "_codon_frequency_overall.csv"
#     ),
#     row.names = F
#   )
#   
# }

# Within signal peptides --------------------------------------------------

for (org in seq_len(ncol(orthologues))) {
  result = c()
  
  # Prepare empty result matrix
  result = matrix("", nrow(orthologues), length(unit_codons) + 1)
  
  for (i in seq_len(nrow(orthologues))) {
    # Find cDNA of interest
    which_cDNA = which(cDNA[[org]]$ID == orthologues[i, org])
    if (length(which_cDNA) > 0) {
      # Get the sequence
      cdna = s2c(cDNA[[org]]$SEQUENCE[which_cDNA])
      cdna_id = cDNA[[org]]$cDNA_ID[which_cDNA]
      # Get the predicted length of signal peptide and check if it exists
      sp_length = SP[[org]][which(SP[[org]][, 2] == cdna_id), 3]
      if (length(sp_length) > 0) {
        # Get the starting and stopping positions of the translated sequence
        codon_start = c()
        codon_start = cDNA[[org]]$START[which_cDNA][[1]]
        codon_stop = c()
        codon_stop = cDNA[[org]]$STOP[which_cDNA][[1]]
        vec = c()
        # Concatenate all translated region
        for (j in seq_len(length(codon_start))) {
          vec = append(x = vec,
                       values = seq(codon_start[j], codon_stop[j], by = 1))
          
        }
        # Get only the translated nucleotides
        translated = cdna[vec]
        # Trim the sequence to the preticted length of signal peptide * 3
        translated = translated[1:(sp_length * 3)]
        saar_start = SAAR[[org]][which(SAAR[[org]][, 2] == cdna_id), 4]
        saar_length = SAAR[[org]][which(SAAR[[org]][, 2] == cdna_id), 3]
        if (length(saar_start) > 0 && saar_start > 0) {
          translated = translated[c(seq(1, (saar_start * 3), 1), seq((saar_start * 3 + 1) + (saar_length * 3), sp_length *
                                                                       3))]
        }
        # In case of faulty ("N") nucleotide - replace it with a random one
        translated[which(translated == "N")] = nucleotides[sample(1:4, 1)]
        # Sequence to codons and get the frequencies
        codons_freq = table(as.character(codons(DNAString(
          c2s(translated)
        ))))
        
        # Just unit of interest codons codons
        avail_unit_codons_freq = codons_freq[unit_codons]
        
        # Combine results into a matrix
        result[i, ] = append(cDNA[[org]]$cDNA_ID[which_cDNA], as.vector(avail_unit_codons_freq))
      }
    }
  }
  colnames(result) = append("cDNA_ID",unit_codons)
  # Write the result to .csv file
  write.csv(
    result,
    file = paste(
      sep = "",
      getwd(),
      "/results/codon_frequency/",
      OoI,
      "_orthologues/",
      organisms_names[org],
      "_codon_frequency_SP.csv"
    ),
    row.names = F
  )
  
}


# Overall excluding signal peptides ---------------------------------------

for (org in seq_len(ncol(orthologues))) {
  result = c()
  
  # Prepare empty result matrix
  result = matrix("", nrow(orthologues), length(unit_codons) + 1)
  
  for (i in seq_len(nrow(orthologues))) {
    # Find cDNA of interest
    which_cDNA = which(cDNA[[org]]$ID == orthologues[i, org])
    if (length(which_cDNA) > 0) {
      # Get the sequence
      cdna = s2c(cDNA[[org]]$SEQUENCE[which_cDNA])
      cdna_id = cDNA[[org]]$cDNA_ID[which_cDNA]
      # Get the predicted length of signal peptide and check if it exists
      # Get the starting and stopping positions of the translated sequence
      codon_start = c()
      codon_start = cDNA[[org]]$START[which_cDNA][[1]]
      codon_stop = c()
      codon_stop = cDNA[[org]]$STOP[which_cDNA][[1]]
      vec = c()
      # Concatenate all translated region
      for (j in seq_len(length(codon_start))) {
        vec = append(x = vec,
                     values = seq(codon_start[j], codon_stop[j], by = 1))
        
      }
      # Get only the translated nucleotides
      translated = cdna[vec]
      sp_length = SP[[org]][which(SP[[org]][, 2] == cdna_id), 3]
      if (length(sp_length) > 0) {
        # Cut the signal peptide from the sequence
        translated = translated[(1 + 3 * sp_length):length(translated)]
      }
      # In case of faulty ("N") nucleotide - replace it with a random one
      translated[which(translated == "N")] = nucleotides[sample(1:4, 1)]
      # Sequence to codons and get the frequencies
      codons_freq = table(as.character(codons(DNAString(
        c2s(translated)
      ))))
      
      # Just unit of interest codons codons
      avail_unit_codons_freq = codons_freq[unit_codons]
      
      # Combine results into a matrix
      result[i, ] = append(cDNA[[org]]$cDNA_ID[which_cDNA], as.vector(avail_unit_codons_freq))
    }
  }
  colnames(result) = append("cDNA_ID",unit_codons)
  # Write the result to .csv file
  write.csv(
    result,
    file = paste(
      sep = "",
      wd,
      "/results/codon_frequency/",
      OoI,
      "_orthologues/",
      organisms_names[org],
      "_codon_frequency_excluding_SP.csv"
    ),
    row.names = F
  )
  
}


# Within SAARs ------------------------------------------------------------

for (org in seq_len(ncol(orthologues))) {
  result = c()
  
  # Prepare empty result matrix
  result = matrix("", nrow(orthologues) * 2, length(unit_codons) + 1)
  
  for (i in seq_len(nrow(orthologues))) {
    # Find cDNA of interest
    which_cDNA = which(cDNA[[org]]$ID == orthologues[i, org])
    if (length(which_cDNA) > 0) {
      # Get the sequence
      cdna = s2c(cDNA[[org]]$SEQUENCE[which_cDNA])
      cdna_id = cDNA[[org]]$cDNA_ID[which_cDNA]
      # Get the start point and length of SAAR in the predicted signal peptide and check if it exists
      saar_start = SAAR[[org]][which(SAAR[[org]][, 2] == cdna_id), 4]
      saar_length = SAAR[[org]][which(SAAR[[org]][, 2] == cdna_id), 3]
      if (length(saar_start) > 0) {
        # Get the starting and stopping positions of the translated sequence
        codon_start = c()
        codon_start = cDNA[[org]]$START[which_cDNA][[1]]
        codon_stop = c()
        codon_stop = cDNA[[org]]$STOP[which_cDNA][[1]]
        vec = c()
        # Concatenate all translated region
        for (j in seq_len(length(codon_start))) {
          vec = append(x = vec,
                       values = seq(codon_start[j], codon_stop[j], by = 1))
          
        }
        # Get only the translated nucleotides
        translated = cdna[vec]
        # Cut the signal peptide from the sequence
        translated = translated[(3 * saar_start + 1):(3 * saar_start + 1 + (saar_length * 3 - 1))]
        # In case of faulty ("N") nucleotide - replace it with a random one
        translated[which(translated == "N")] = nucleotides[sample(1:4, 1)]
        
        # Sequence to codons and get the frequencies
        codons_freq = table(as.character(codons(DNAString(
          c2s(translated)
        ))))
        
        # Just unit of interest codons codons
        avail_unit_codons_freq = codons_freq[unit_codons]
        
        # Combine results into a matrix
        result[i, ] = append(cDNA[[org]]$cDNA_ID[which_cDNA], as.vector(avail_unit_codons_freq))
      }
    }
  }
  colnames(result) = append("cDNA_ID",unit_codons)
  # Write the result to .csv file
  write.csv(
    result,
    file = paste(
      sep = "",
      getwd(),
      "/results/codon_frequency/",
      OoI,
      "_orthologues/",
      organisms_names[org],
      "_codon_frequency_SAAR.csv"
    ),
    row.names = F
  )
  
}

gc()
