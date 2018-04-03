#USAGE: RScript prep2revtrans.R <organisms_Latin_names> 

library(seqinr)
args = commandArgs(trailingOnly = TRUE)
wd = getwd()
organisms_names = args
setwd(paste(wd, "/data/readData", sep = ""))
#read the saved cDNA data
cDNA=list()
for(i in seq_len(length(organisms_names))){
  cDNA[[i]] = readRDS(paste("readData_cDNA_",organisms_names[i],".rds",sep = ""))[[1]]
}
names(cDNA) = organisms_names
setwd(wd)
#extract just the codong sequences based on start and stop codons 
for(org in seq(1,length(cDNA))){
  print(org)
  for (i in seq(1,length(cDNA[[org]]$STOP))){
    codon_start = c()
    codon_start = cDNA[[org]]$START[i][[1]]
    codon_stop = c()
    codon_stop = cDNA[[org]]$STOP[i][[1]]
    vec = c()
    for (j in seq_len(length(codon_start))) {
      vec = append(x = vec,
                   values = seq(codon_start[j], codon_stop[j], by = 1))
      
    }
    # Get only the translated nucleotides
    cDNA[[org]]$cdna[i] = c2s(s2c(cDNA[[org]]$cdna[i])[vec])
  }
}
#write to file which will be used as a RevTrans input
setwd(wd)
prep2revtrans = list();
for(org in 1:length(cDNA)){
  prep2revtrans[[org]] = do.call(cbind,cDNA[[org]][,c(2,5)])
  for (i in seq(1,length(prep2revtrans[[org]][,2]),1)){
    prep2revtrans[[org]][i,2] = prep2revtrans[[org]][i,2]
  }
  for (i in seq(1, length(prep2revtrans),1)){
    write.csv(row.names = F,x = prep2revtrans[[i]],file = paste(wd,"/results/revtrans_",organisms_names[1],"/prep2revtrans_full","_",organisms_names[org],".csv",sep = ""))
  }
}
