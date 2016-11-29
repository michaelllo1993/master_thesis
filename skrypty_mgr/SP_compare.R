require(seqinr);
wd = getwd()
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3){
  print("Wrong number of arguments passed!")
  stop()
}

#USAGE: Rscript name filename (prep2comp_<OoI>_all.csv) OoI(ensembl exclusive code) 1Ortho(ensembl exclusive code)

filename = args[1]
OoI = args[2]
Ortho = args[3]

# reading data
data=read.csv(filename, sep=",",header=FALSE);

# extracting just the data of interest
pattern = paste("^",Ortho,sep = "")
Ortho_rows = which(grepl(as.vector(data[,3]),pattern = pattern,perl = T) == TRUE)
data = data[Ortho_rows,]

changes=c("");
for (i in seq(1,length(data[,1]))) {
  dane=(as.character(data[i,2]));
  dane.1=s2c(as.character(data[i,4]));
  
  index=grepRaw(c2s(rep("L",5)),(dane));#lrun position
  if (length(index)!=0) {
  a=rle(s2c(dane));
  length=a$lengths[which(a$values=="L")][which.max(a$lengths[which(a$values=="L")])]#Lrun length
   change=dane.1[index:(index+length-1)];
  changes=append(changes,change);
}}

changes_table_percentage=round((sort(table(changes))/length(changes)*100),2);
changes_table_percentage

changes_table=sort(table(changes));
changes_table

changes_all=c("");
for (i in seq(2,length(data[,1]))) {
  dane=s2c(as.character(data[i,2]));
  dane.1=s2c(as.character(data[i,4]));
  
    change_all=dane.1[which(dane=="L")];
    changes_all=append(changes_all,change_all);
}

all_changes_table_percentage=round((sort(table(changes_all))/length(changes_all)*100),2);
all_changes_table_percentage

all_changes_table=(sort(table(changes_all)))
all_changes_table
