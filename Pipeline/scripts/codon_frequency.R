
# Load libraries ----------------------------------------------------------

library(Biostrings)
library(seqinr)

# Load data ---------------------------------------------------------------

args = commandArgs(trailingOnly = TRUE)

wd = getwd()
setwd(paste(wd, "/data", sep = ""))
load("readData.RData")
setwd("..")


# Processing the input mapper file ----------------------------------------

# Filtering only for rows with all available orthologues
original_mapper_name = args[1]
repeat_unit = args[2]
organisms_names = args[c(-1,-2)]
unit_codons = names(which(GENETIC_CODE==repeat_unit));
nucleotides = c("A","T","G","C")

# Build the corrected  mapper file name
common_mapper_name = paste(strsplit(original_mapper_name,"\\.")[[1]][1],"_common.",strsplit(original_mapper_name,"\\.")[[1]][2],sep = "")
# Grep only the lines without 'NULL' in them - ones where all orthologues are available
command = paste("grep -Ev 'NULL' ",original_mapper_name," > ",common_mapper_name,sep = "")
system(command)
# Read corrected file
orthologues = read.csv(common_mapper_name,header = F)[,-1]
colnames(orthologues) = organisms_names 


# The analysis ------------------------------------------------------------

for (org in seq_len(ncol(orthologues))) {
  result = c();
  # Prepare empty result matrix
  result = matrix("",nrow(orthologues) * 2,length(unit_codons)+1);
  for (i in seq_len(nrow(orthologues))) {
    # Find cDNA of interest
    which_cDNA = which(cDNA[[org]]$ID == orthologues[i,org])
    if(length(which_cDNA)>0){
      # Get the sequence
      cdna = s2c(cDNA[[org]]$SEQUENCE[which_cDNA])
      # Get the starting and stopping positions of the translated sequence 
      codon_start = c(); codon_start = cDNA[[org]]$START[which_cDNA][[1]]
      codon_stop = c(); codon_stop = cDNA[[org]]$STOP[which_cDNA][[1]]
      vec = c()
      # Concatenate all translated region
      for (j in seq_len(length(codon_start))) {
        vec = append(x = vec,values = seq(codon_start[j],codon_stop[j],by = 1));
      }
      # Get only the translated nucleotides
      translated = cdna[vec];
      # In case of faulty ("N") nucleotide - replace it with a random one
      translated[which(translated == "N")] = nucleotides[sample(1:4, 1)]
      # Sequence to codons and get the frequencies
      codons_freq = table(as.character(codons(DNAString(c2s(translated)))));
      # Just unit of interest codons codons
      avail_unit_codons_freq = codons_freq[unit_codons];
      # Combine results into a matrix
      result[(i * 2) - 1,] = append(names(avail_unit_codons_freq),cDNA[[org]]$cDNA_ID[which_cDNA])
      result[i * 2,] = append(as.vector(avail_unit_codons_freq),0)
    }
  }
  # Write the result to .csv file
  write.csv(result,file = paste(sep = "",getwd(),"/",organisms_names[org],"_corrected_ortho_only_result.csv"));
}
  























