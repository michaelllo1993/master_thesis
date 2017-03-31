setwd("~/Uczelnia/MGR/praca_magisterska/codon_analysis/cDNA/");
library(seqinr)
myFiles <- list.files(pattern = "*.csv");
cDNA = list();
cDNA_names = c();
for (k in 1:length(myFiles)) {
  cDNA_names[k] = strsplit(myFiles,"[.]")[[k]][1];
  cDNA[[k]] = read.csv(myFiles[k],stringsAsFactors = F);
  print(paste("cDNA data read for:",cDNA_names[k]))
}
cDNA_names[k] = strsplit(myFiles,"[.]")[[k]][1];
names(cDNA) = cDNA_names;


for (org in seq(1,length(cDNA))){
    cDNA[[org]][[3]] = lapply(cDNA[[org]][[3]], function(x) strsplit(x = x, split = ";")[[1]]);#split by semicolon
    cDNA[[org]][[4]] = lapply(cDNA[[org]][[4]], function(x) strsplit(x = x, split = ";")[[1]])
    print(paste("cDNA data processed for:",cDNA_names[org]))
}

trimm2coding = function(cDNA){
  return(c2s(s2c(cDNA[[5]])[seq(min(as.numeric(unlist(cDNA[[3]]))),max(as.numeric(unlist(cDNA[[4]]))),by = 1)]))
}
a = c();
cDNA_test = list(a,a,a,a,a,a,a,a,a);
for(org in seq(1,length(cDNA))){
  for (i in seq(1,length(cDNA[[org]]$STOP))){
    vec=c();
    starting=min(as.numeric(unlist(cDNA[[org]]$START[i])))#coding sequence start
    stopping=max(as.numeric(unlist(cDNA[[org]]$STOP[i])))#coding sequence end
    vec = seq(starting,stopping,by = 1)#numbers of coding nucleotides
    cseq=c2s(s2c(cDNA[[org]]$SEQUENCE[i])[vec]); #proper coding sequence
    cDNA[[org]]$SEQUENCE[i] = cseq
  }
  print(paste("cDNA data processed for",cDNA_names[org]))
}

a = c();
prep2revtrans = list(a,a,a,a,a,a,a,a,a);
for(org in 1:length(cDNA)){
  prep2revtrans[[org]] = as.matrix(cDNA[[org]])[,c(2,5)]
  for (i in seq(1,length(prep2revtrans[[org]][,2]),1)){
    if (length(s2c(prep2revtrans[[org]][i,2]$SEQUENCE)) >= 3*69){
      prep2revtrans[[org]][i,2] = c2s(s2c(prep2revtrans[[org]][i,2]$SEQUENCE)[1:(3*69)])
  }
    else{
      prep2revtrans[[org]][i,2] = prep2revtrans[[org]][i,2]$SEQUENCE
    }
  }
  names(prep2revtrans) = cDNA_names
  for (i in seq(1, length(prep2revtrans),1)){
    write.csv(row.names = F,x = prep2revtrans[[i]],file = paste(getwd(),"/","prep2revtrans","_",cDNA_names[i],".csv",sep = ""))
  }
}
