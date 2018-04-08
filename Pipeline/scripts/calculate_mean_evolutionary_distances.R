
# Load libraries ----------------------------------------------------------

library(seqinr)
library(ape)
require(svglite)
require(msa)

args = commandArgs(trailingOnly = TRUE)

#suppress printing warning messages (which are caused by cbind in this case)
options(warn=-1)
#get command line arguments
algorithm = args[1];args = args[-1]
no_samples = as.numeric(args[1]);args=args[-1]
OoI = args[1]
organisms_names=args
#get working directory
wd = getwd()
#read protein sequences
proteins=list()
setwd(paste(wd, "/data/readData", sep = ""))
for(i in seq_len(length(organisms_names))){
  proteins[[i]] = readRDS(paste("readData_proteins_",organisms_names[i],".rds",sep = ""))[[1]]
}
names(proteins) = organisms_names
# proteins = readRDS(paste(wd, "/data/readData/readData_proteins.rds", sep = ""))
#read mapper with common orthologoues protein IDs
mapper = as.matrix(read.csv(paste(wd,"/data/orthoMappers/",OoI,"_ids_mapper_common.csv",sep = ""), header = T)[,-1])
#get number of rows of the mapper 
no_of_proteins = dim(mapper)[1]
#prepare an empty list for final results
list_distances = list()
row_number=c()
#main loop - 100 iteration = 100 random protein samples
for (rep in seq_len(no_samples)) {
  #draw row number
  row_number[rep] = sample(1:no_of_proteins, 1)
  #get the orthologous IDs
  ids = mapper[row_number[rep],]
  #get the names of organisms associated with them
  ids_names = names(ids)
  
  sequences = c()
  i = 1
  #crucial loop - extracts sequences of the ortohlogous proteins from the row which was drawn. Notably, assures all the orthologues are available
  while (i <= length(ids)) {
    number = (which(proteins[[ids_names[i]]]$ensembl_peptide_id == ids[i]))
    if (length(number) < 1) {
      row_number[rep] = sample(1:no_of_proteins, 1)
      ids = mapper[row_number[rep],]
      sequences = c()
      i = 1
    } else{
      sequences[i] = proteins[[ids_names[i]]]$peptide[number]
      i = i + 1
    }
  }
  
  names(sequences) = ids_names
  #convert output vector with sequences to AAStringSet object with names
  sequences_converted = AAStringSet(sequences,use.names = T)
  #align the sequences
  alignment = msa(sequences_converted,method = algorithm)
  #calculate pairwise distances from the alignment 
  list_distances[[rep]] = dist.alignment(msaConvert(alignment,type = "seqinr::alignment"))
}

#average the results obtained from all iterations
summed = list_distances[[1]]
for(e in seq_len(length(list_distances))){
  if(any(is.nan(list_distances[[e]]))){
    summed = summed
  } else{
    summed = list_distances[[e]] + summed  
  }
}

averaged=summed/length(list_distances)

#write the averaged distance matrix to CSV file
write.csv(x = as.matrix(averaged),file = paste(wd,"/results/",OoI,"_mean_pairwise_evolutionary_distances.csv",sep = ""),row.names = T,col.names = T)
#generate a phylogenetic tree based on averged pairwise distances 
svglite(file=paste(wd,"/results/",OoI,"_phylogenetic_tree.svg",sep = ""),
        width=40, 
        height=30, 
        pointsize=50)
plot(as.phylo(hclust((averaged),method = "ward.D2")), cex = 1, label.offset = 0.01,main = "Phylogenetic tree; hierarchical clusetring using UPGMA metod",sub = 
       paste("Based on alignment of", rep ,"randomly selected protein sequences using",algorithm,"algorithm"),edge.width = 3)
dev.off()
