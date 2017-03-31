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

#import packages
library(foreach)
library(doParallel)

# prepare cDNA sequences (0-69*3) to RevTrans -----------------------------
#setup parallel backend to use 3processors
cl <- makeCluster(3)
registerDoParallel(cl)

#start timecDN
strt <- Sys.time()

########### Split START and STOP values - time consuming ##########
for (org in seq(1,length(cDNA))){
  for(s in seq(1,length(cDNA[[org]][[3]]))) {
    if (length(grep(";",cDNA[[org]][[3]][s])) == 1) {
      cDNA[[org]][[3]][s] = strsplit(cDNA[[org]][[3]][s][[1]],";");#split by semicolon
      cDNA[[org]][[4]][s] = strsplit(cDNA[[org]][[4]][s][[1]],";");#split by semicolon
    }
    else{
      cDNA[[org]][[3]][s] = cDNA[[org]][[3]][s][[1]];
      cDNA[[org]][[4]][s] = cDNA[[org]][[4]][s][[1]];
    }
  }
  print(paste("cDNA data processed for:",cDNA_names[org]))
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
  names(prep2revtrans) = prot_names
  for (i in seq(1, length(prep2revtrans),1)){
    write.csv(row.names = F,x = prep2revtrans[[i]],file = paste(getwd(),"/","prep2revtrans","_",cDNA_names[i],".csv",sep = ""))
  }
}
