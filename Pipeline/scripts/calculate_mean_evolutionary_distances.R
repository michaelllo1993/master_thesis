
# Load libraries ----------------------------------------------------------

library(seqinr)
library(ape)
require(svglite)

#suppress printing warning messages (which are caused by cbind in this case)
options(warn=-1)
#get working directory
wd = getwd()
#define function to be used using apply
test <- function(x) {
  row = strsplit(x, " ")[[1]]
  row = as.numeric(row)
  return(row)
}
#read protein sequences
proteins = readRDS(paste(wd, "/data/readData/readData_proteins.rds", sep = ""))
#read mapper with common orthologoues protein IDs
mapper = as.matrix(read.csv("bos_taurus_ids_mapper_common.csv", header = T)[,-1])
#get number of rows of the mapper 
no_of_proteins = dim(mapper)[1]
#prepare an empty list for final results
list_distances = list()
#remove any remainders after previous runs
system("rm toAlignMEGA.fasta Alignment* Distance*")
#main loop - 100 iteration = 100 random protein samples
for (rep in seq_len(100)) {
  #draw row number
  row_number = sample(1:no_of_proteins, 1)
  #get the orthologous IDs
  ids = mapper[row_number,]
  #get the names of organisms associated with them
  ids_names = names(ids)
  
  sequences = list()
  i = 1
  #crucial loop - extracts sequences of the ortohlogous proteins from the row which was drawn. Notably, assures all the orthologues are available
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
  #write the  extracted sequences to the FASTA file
  write.fasta(sequences = sequences,
              names = ids_names,
              file.out = "toAlignMEGA.fasta")
  #align the sequences written to the FASTA file
  system("megacc -a clustal_align_protein.mao -d toAlignMEGA.fasta -o Alignment -n -s",wait = T)
  #calculate pairwise distances among the sequences
  system(
    "megacc -a distance_estimation_pairwise_protein.mao -d Alignment.meg -o Distance",wait = T
  )
  #filter just the names in the distance matrix and save to another file
  system("grep '^\\[[0-9]\\]\\s#' Distance.meg  > Distance_names.meg",wait = T)
  #filter just the values in the distance matrix and save to another file
  system("grep '^\\[[0-9]\\]\\s*[0-9]' Distance.meg  > Distance_values.meg",wait = T)
  #parse the output
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
  #parse the output
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
  #save the distance matrix as a list element
  list_distances[[rep]] = distances
  #delete the files from the current iteration
  system("rm toAlignMEGA.fasta Alignment* Distance*")
}

#average the results obtained from all iterations
summed = list_distances[[1]]
summed[length(summed) > 0] = 0
for(e in seq_len(length(list_distances))){
  summed = list_distances[[e]] + summed
}

averaged=summed/length(list_distances)

#write the averaged distance matrix to CSV file
write.csv(x = averaged,file = paste(wd,"/mean_pairwise_evolutionary_distances.csv",sep = ""),row.names = T,col.names = T)
#generate a phylogenetic tree based on averged pairwise distances 
svglite(file="phylogenetic_tree.svg",
        width=40, 
        height=30, 
        pointsize=50)
plot(as.phylo(hclust(as.dist(averaged),method = "average")), cex = 1, label.offset = 0.01,main = "Phylogenetic tree; hierarchical clusetring using UPGMA metod",sub = 
       paste("Based on alignment of", rep ,"randomly selected protein sequences"),edge.width = 3)
dev.off()
