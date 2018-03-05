wd = getwd()
require(seqinr);
require(tidyverse)
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2){
  print("Wrong number of arguments passed!")
  stop()
}

filename = args[1]
saar = args[2];

#USAGE: RScprit name input_file_name SAARtyp
# reading data
data_orig=read.csv(filename, sep=",",header=FALSE);
organismName=paste(strsplit(strsplit(filename,"/")[[1]][length(strsplit(filename,"/")[[1]])],"_")[[1]][c(1,2)][1],strsplit(strsplit(filename,"/")[[1]][length(strsplit(filename,"/")[[1]])],"_")[[1]][c(1,2)][2],sep = "_")
dirOrganismName=paste(organismName,"_orthoAAstats",sep = "")
unique_IDs=unique(sapply(data_orig[c(1:1000),3], function(x) gsub('[[:digit:]]+', '', x)))
for(Ortho in unique_IDs){
  # extracting just the data of interest
  pattern = paste("^",Ortho,sep = "")
  Ortho_rows = which(grepl(as.vector(data_orig[,3]),pattern = pattern,perl = T) == TRUE)
  data = data_orig[Ortho_rows,]
  OoI = strsplit(as.vector(data[1,1]),split = "0")[[1]][1]
  changes=c("");
  
  for (i in seq(2,length(data[,1]))) {
    dane=(as.character(data[i,2]));
    dane.1=s2c(as.character(data[i,4]));
    index=grepRaw(c2s(rep(saar,5)),(dane));#lrun position
    if (length(index)!=0) {
    a=rle(s2c(dane));
    length=a$lengths[which(a$values==saar)][which.max(a$lengths[which(a$values==saar)])]#Lrun length
     change=dane.1[index:(index+length-1)];
    changes=append(changes,change);
  }}
  
  changes_table_percentage=round((sort(table(changes))/length(changes)*100),2);
  changes_table=sort(table(changes));
  
  df_results_percentage = as.data.frame(changes_table_percentage)
  df_results = as.data.frame(changes_table)
  df_results_binded = cbind(df_results,df_results_percentage)
  write.csv(df_results_binded,file = paste(wd,"/",dirOrganismName,"/",OoI,"_",Ortho,"_saar_sigp_aa_on_",saar,"_position.csv",sep = ""),row.names = F)
  changes_all=c("");
  for (i in seq(2,length(data[,1]))) {
    dane=s2c(as.character(data[i,2]));
    dane.1=s2c(as.character(data[i,4]));
    
      change_all=dane.1[which(dane==saar)];
      changes_all=append(changes_all,change_all);
  }
  
  all_changes_table_percentage=round((sort(table(changes_all))/length(changes_all)*100),2);
  all_changes_table=(sort(table(changes_all)))
  
  df_all_results_percentage = as.data.frame(all_changes_table_percentage)
  df_all_results = as.data.frame(all_changes_table)
  df_all_results_binded = (cbind(df_all_results,df_all_results_percentage))
  write.csv(df_all_results_binded,file = paste(wd,"/",dirOrganismName,"/",OoI,"_",Ortho,"_aa_on_",saar,"_position.csv",sep = ""),row.names = F)
}
