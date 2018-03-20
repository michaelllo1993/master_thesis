#USAGE: RScript prep2revtrans.R <organisms_Latin_names> 

library(seqinr)
wd = getwd()
cDNA = readRDS(paste(wd,"/data/readData/readData_cDNA.rds",sep = ""))
organisms=names(cDNA)

for(org in seq(1,length(cDNA))){
  for (i in seq(1,length(cDNA[[org]]$STOP))){
    indices = seq(min(as.numeric(unlist(cDNA[[org]]$START[i]))),max(as.numeric(unlist(cDNA[[org]]$STOP[i]))),by = 1)
    cDNA[[org]]$SEQUENCE[i] = c2s(s2c(cDNA[[org]]$SEQUENCE[i])[indices])
  }
}

setwd(wd)
prep2revtrans = list();
for(org in 1:length(cDNA)){
  prep2revtrans[[org]] = do.call(cbind,cDNA[[org]][,c(2,5)])
  for (i in seq(1,length(prep2revtrans[[org]][,2]),1)){
    prep2revtrans[[org]][i,2] = prep2revtrans[[org]][i,2]
  }
  for (i in seq(1, length(prep2revtrans),1)){
    write.csv(row.names = F,x = prep2revtrans[[i]],file = paste(getwd(),"/revtrans_ready/prep2revtrans_full","_",organisms[org],".csv",sep = ""))
  }
}
