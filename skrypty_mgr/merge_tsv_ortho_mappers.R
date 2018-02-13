# USAGE Rscript merge_tsv_ortho_mappers.R file_1 file_2 output_name
wd = getwd()
require(methods)
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3){
  print("Wrong number of arguments passed!")
  stop()
}
print("Reading files.")
file = read.csv(paste(paste(wd,"/",sep = ""),args[1],sep=""),sep = "\t",header = T);
file_1 = read.csv(paste(paste(wd,"/",sep = ""),args[2],sep = ""),sep = "\t",header = T);
print("Files read.")
file = as.matrix(file); 
file_1 = as.matrix(file_1);
A = file[which(!duplicated(file[,1])),1]; # find non duplicated IDs in 1st column of the larger file
file=file[which(!duplicated(file[,1])),]; # leave only these rows
indices = c();
print("Identifying unique rows.")
for (i in seq(1, length(A),1)){
  indices = append(indices,which(file_1[,1] == A[i])[1]) # which rows of the other file are these?
}
print("Unique rows identified.")
file_1 = file_1[indices,] # leave only these rows in smaller file
if(any(!(file_1[,1] == file[,1]))){
  print("Mismatch encountered!");
  stop()
}; # check for any mismatches
to_merge = file_1[,-1]; # leave first column (already in the larger file)
merged = cbind2(file, to_merge) # merge 'em
write.table(merged,file = (paste(paste(wd,"/",sep = ""),args[3],sep="")),sep = "\t",row.names = F)


