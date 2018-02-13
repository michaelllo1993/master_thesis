#USAGE: RScript prep2revtrans.R <organism_of_interest_Latin_name> 

library(seqinr)
wd = getwd()
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1){
  print("Wrong number of arguments passed!")
  stop()
}
search_string = paste("cDNA_",args[1],".csv",sep = "")
myFile <- list.files(pattern = search_string);
if(length(myFile) != 1){
  print(paste("There must be a single csv file named",search_string,"in the current directory"))
  stop()
}
cDNA = list();
cDNA_names = c();
for (k in 1:length(myFile)) {
  cDNA_names[k] = strsplit(myFile,"[.]")[[k]][1];
  cDNA[[k]] = read.csv(myFile[k],stringsAsFactors = F);
  print(paste("cDNA data read for:",args[1]))
}
cDNA_names[k] = strsplit(myFile,"[.]")[[k]][1];
names(cDNA) = cDNA_names;

for (org in seq(1,length(cDNA))){
    cDNA[[org]][[3]] = lapply(cDNA[[org]][[3]], function(x) strsplit(x = x, split = ";")[[1]]);#split by semicolon
    cDNA[[org]][[4]] = lapply(cDNA[[org]][[4]], function(x) strsplit(x = x, split = ";")[[1]])
}

for(org in seq(1,length(cDNA))){
  for (i in seq(1,length(cDNA[[org]]$STOP))){
    indices = seq(min(as.numeric(unlist(cDNA[[org]]$START[i]))),max(as.numeric(unlist(cDNA[[org]]$STOP[i]))),by = 1)
    cDNA[[org]]$SEQUENCE[i] = c2s(s2c(cDNA[[org]]$SEQUENCE[i])[indices])
  }
  print(paste("cDNA data processed for",args[1]))
}

a = c();
prep2revtrans_trimmed = list(a);
for(org in 1:length(cDNA)){
  prep2revtrans_trimmed[[org]] = as.matrix(cDNA[[org]])[,c(2,5)]
for (i in seq(1,length(prep2revtrans_trimmed[[org]][,2]),1)){
    if (length(s2c(prep2revtrans_trimmed[[org]][i,2]$SEQUENCE)) >= 3*69){
      prep2revtrans_trimmed[[org]][i,2] = c2s(s2c(prep2revtrans_trimmed[[org]][i,2]$SEQUENCE)[1:(3*69)])
  }
    else{
      prep2revtrans_trimmed[[org]][i,2] = prep2revtrans_trimmed[[org]][i,2]$SEQUENCE
    }
  }
  names(prep2revtrans_trimmed) = cDNA_names
  for (i in seq(1, length(prep2revtrans_trimmed),1)){
    write.csv(row.names = F,x = prep2revtrans_trimmed[[i]],file = paste(getwd(),"/","prep2revtrans_trimmed","_",args[1],".csv",sep = ""))
  }
  print(paste("trimmed cDNA data written for",args[1]))
}

a = c();
prep2revtrans = list(a);
for(org in 1:length(cDNA)){
  prep2revtrans[[org]] = as.matrix(cDNA[[org]])[,c(2,5)]
  for (i in seq(1,length(prep2revtrans[[org]][,2]),1)){
    prep2revtrans[[org]][i,2] = prep2revtrans[[org]][i,2]$SEQUENCE
  }
  names(prep2revtrans) = cDNA_names
  for (i in seq(1, length(prep2revtrans),1)){
    write.csv(row.names = F,x = prep2revtrans[[i]],file = paste(getwd(),"/","prep2revtrans_full","_",args[1],".csv",sep = ""))
  }
  print(paste("full cDNA data written for",args[1]))
}
