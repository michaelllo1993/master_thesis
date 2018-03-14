library(seqinr)
options(warn=-1)
# args = commandArgs(trailingOnly = TRUE)
# setwd("~/Uczelnia/MGR/master_thesis/Pipeline/")
wd = getwd()

test <- function(x) {
  row = strsplit(x, " ")[[1]]
  row = as.numeric(row)
  return(row)
}

proteins = readRDS(paste(wd, "/data/readData/readData_proteins.rds", sep = ""))
mapper = as.matrix(read.csv("bos_taurus_ids_mapper_common.csv", header = T)[,-1])
no_of_proteins = dim(mapper)[1]
list_distances = list()
system("rm toAlignMEGA.fasta Alignment* Distance*")
for (rep in seq_len(100)) {
  row_number = sample(1:no_of_proteins, 1)
  ids = mapper[row_number,]
  ids_names = names(ids)
  
  sequences = list()
  i = 1
  while (i <= length(ids)) {
    number = (which(proteins[[ids_names[i]]]$ID == ids[i]))
    if (length(number) < 1) {
      row_number = sample(1:no_of_proteins, 1)
      ids = mapper[row_number,]
      sequences = list()
      i = 1
    } else{
      sequences[[i]] = s2c(proteins[[ids_names[i]]]$SEQUENCE[number])
      i = i + 1
    }
  }
  
  names(sequences) = ids_names
  write.fasta(sequences = sequences,
              names = ids_names,
              file.out = "toAlignMEGA.fasta")
  
  system("megacc -a clustal_align_protein.mao -d toAlignMEGA.fasta -o Alignment -n -s",wait = T)
  system(
    "megacc -a distance_estimation_pairwise_protein.mao -d Alignment.meg -o Distance",wait = T
  )
  system("grep '^\\[[0-9]\\]\\s#' Distance.meg  > Distance_names.meg",wait = T)
  system("grep '^\\[[0-9]\\]\\  [0-9]' Distance.meg  > Distance_values.meg",wait = T)
  
  names_file = read.delim("Distance_names.meg",
                          header = F,
                          stringsAsFactors = F)[, 1]
  names = as.character(sapply(names_file, function(x)
    gsub(
      pattern = '\\[\\d\\]\\s#',
      replacement = "",
      x = x,
      perl = T
    )))
  distance_file = read.delim("Distance_values.meg",
                             header = F,
                             stringsAsFactors = F)[, 1]
  pre_distances = as.character(sapply(distance_file, function(x)
    gsub(
      pattern = '\\[\\d\\]\\ \\ ',
      replacement = "",
      x = x,
      perl = T
    )))
  
  distances_list = lapply(pre_distances, test)
  no_org = length(distances_list)
  dists = lapply(distances_list, function(x) x[1:no_org])
  distances = matrix(unlist(dists), ncol = no_org, byrow = TRUE)
  distances = cbind(distances,(rep(NA,no_org-1)))
  distances = rbind(t(rep(NA,no_org+1)),distances)
  rownames(distances) = names;colnames(distances)=names
  list_distances[[rep]] = distances
  system("rm toAlignMEGA.fasta Alignment* Distance*")
}
summed = list_distances[[1]]
summed[length(summed) > 0] = 0
for(e in seq_len(length(list_distances))){
  summed = list_distances[[e]] + summed
}

averaged=summed/length(list_distances)

write.csv(x = averaged,file = paste(wd,"/mean_pairwise_evolutionary_distances.csv",sep = ""),row.names = T,col.names = T)


