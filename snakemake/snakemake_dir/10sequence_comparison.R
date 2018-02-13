wd = getwd()
require(seqinr);
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3){
  print("Wrong number of arguments passed!")
  stop()
}

filename = args[1]
saar = args[2];
Ortho = args[3]

#USAGE: RScprit name input_file_name SAARtype
# reading data
data=read.csv(filename, sep=",",header=FALSE);

# extracting just the data of interest
pattern = paste("^",Ortho,sep = "")
Ortho_rows = which(grepl(as.vector(data[,3]),pattern = pattern,perl = T) == TRUE)
data = data[Ortho_rows,]
OoI = strsplit(as.vector(data[1,1]),split = "0")[[1]][1]
changes=c("");

for (i in seq(1,length(data[,1]))) {
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
#write.csv(df_results_binded,file = paste(wd,"/",OoI,"_",Ortho,"saar_sigp_aa_on_l_position.csv",sep = ""),row.names = F)
print(paste(saar,"-SAAR in Signal Peptide, amino acids on ", saar, " position"))
print(df_results_binded, quote = TRUE, row.names = FALSE)
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
df_all_results_binded = cbind(df_all_results,df_all_results_percentage)
#write.csv(df_all_results_binded,file = paste(wd,"/",OoI,"_",Ortho,"_aa_on_l_position.csv",sep = ""),row.names = F)
print(paste("Amino acids on", saar, " position"))
print(df_all_results_binded, quote = TRUE, row.names = FALSE)
