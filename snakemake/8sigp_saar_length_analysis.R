wd = getwd()
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1){
  print("Wrong number of arguments passed!")
  stop()
}

#USAGE: RScprit name input_file_name

SPL=read.csv(args[1],sep=",",header = F);#Read csv file
number_of_orgs = ((dim(SPL)[2])/3)
OoI_Ortho_LSAAR_length_diff =list()
OoI_Ortho_LSAAR_data =list()
OoI_Ortho_SP_length_diff =list()
OoI_Ortho_SP_data =list()
lengthened_indices=list()
correlations=list()
myNames = c()
for(i in seq(1,number_of_orgs,by = 1)){
  test_index = which(as.vector(SPL[,(i*3-2)])!="NULL")[1]
  if(is.na(test_index)){
    test_index = "UNKNOWN"
  }
  myNames[i] = strsplit(as.character(SPL[test_index,(i*3-2)]),split = "0")[[1]][1]
}
for(org in seq(2,number_of_orgs,by = 1)){
  OoI_Ortho_LSAAR_length_diff[[org]]= SPL[,3] - SPL[,org*3]; #length(homo sapiens LSAAR) - length(othologous LSAAR) -> positive number = LSAAR got longer -> negative number = LSAAR got shorter
  OoI_Ortho_SP_length_diff[[org]]= SPL[,2] - SPL[,(org*3)-1]; #length(homo sapiens SP) - length(othologous SP) -> positive number = SP got longer -> negative number = SP got shorter
  lengthened_indices[[org]] = which(OoI_Ortho_LSAAR_length_diff[[org]]>0)
  OoI_Ortho_LSAAR_data[[org]]=OoI_Ortho_LSAAR_length_diff[[org]][lengthened_indices[[org]]]
  OoI_Ortho_SP_data[[org]]=OoI_Ortho_SP_length_diff[[org]][lengthened_indices[[org]]]
  correlations[[org]] = cor(OoI_Ortho_LSAAR_data[[org]],OoI_Ortho_SP_data[[org]],method = "pearson")
}
correlations[[1]] = "Not applicable"
names(correlations)=myNames
print("Pearson correlation values for lengthened SAARs and signal peptides for following organisms:")
print(correlations)
