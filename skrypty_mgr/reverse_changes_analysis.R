wd = getwd()
require(seqinr)
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2){
  print("Wrong number of arguments passed!")
  stop()
}
#USAGE: RScprit name infile_name organism_exclusive_code
pre_OoI = strsplit(strsplit(args[1],"/")[[1]][length(strsplit(args[1],"/")[[1]])],"_")[[1]][c(2,3)]
OoI = paste(pre_OoI[1],pre_OoI[2],sep = "_")
command = paste("/home/mstolarczyk/Uczelnia/MGR/praca_magisterska/skrypty_mgr/do_analizy_odwrotnej.sh",sep = " ",args[1], args[2])
system(command = command)

data=read.csv(paste(wd,"/tmp_",args[2],".csv",sep=""),sep=",",header=FALSE);
changes=c("");
# Columns in tmp_<organism_exclusive_code> file: lower_name,lower_seq,higher_name,higher_seq
# Consequently dane contains sequence from lower organism and dane.1 contains sequence from higher organism (reverse changes are anayzed; changes frequency from leucine in lower organism to other AA in higher one)
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
print(paste("Result saved to: ",wd,"/",args[2],"_",OoI,"_changes_percentage.csv",sep=""))
write.table(changes_table_percentage,file = paste(wd,"/",args[2],"_",OoI,"_changes_percentage.csv",sep=""),sep = ",",row.names = T,col.names = F)


