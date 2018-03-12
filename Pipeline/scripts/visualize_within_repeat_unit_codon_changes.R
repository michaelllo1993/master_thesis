
# Load libraries ----------------------------------------------------------

require(gplots)
require(seqinr)
require(Biostrings)
require(fmsb)
require(svglite)

# Get the arguments -------------------------------------------------------

#get working directory
wd = getwd()
#get command line args
args = commandArgs(trailingOnly=TRUE)
#get the last argument and delete it
repeat_unit = args[length(args)]; args = args[-length(args)]
#get the last argument and delete it
OoI =  args[length(args)]; args = args[-length(args)]
#all other arguments are files with data to be visualized
codon_changes_files = args
#infer the output directory from one input files
dir=strsplit(codon_changes_files[1],"/")[[1]][1]
#get all codons 
all_codon_names = append("---",sort(names(GENETIC_CODE)))
#get repeat unit codons
unit_codons = names(which(GENETIC_CODE == repeat_unit))
#create an empty list for the results
within_repeat_unit_list = list()
#create regions vector
regions = c("SAAR","repeat_unit")
i=1

# Visualization -----------------------------------------------------------

for (file_name in codon_changes_files){
  #get the other organism name which matches the input file name - other organism of interest
  other_OoI = gsub(".csv","",file_name)
  splitted=strsplit(other_OoI,"_")[[1]]
  other_OoI=paste(splitted[seq(length(splitted)-1,length(splitted))],collapse = "_")
  print(other_OoI)
  #get the name of the region of interest which matches to input file name 
  RoI = names(which(sapply(regions, function(x) grepl(x,file_name))))
  #construct the name of the output image
  image_name = gsub(".csv",".svg",basename(file_name))
  print(image_name)
  #read the data in
  occurrences=data.matrix(read.csv(file_name,header = T,sep = ",",row.names = all_codon_names))[,-1]
  #change 0 to NAs
  occurrences[which(occurrences==0)] = NA
  #keep just the codons encoding the repeat unit
  within_repeat_unit = occurrences[unit_codons,]
  within_repeat_unit_list[[i]] = within_repeat_unit
  setwd(paste(dir,"/visualization_",OoI,sep = ""))
  #create image and save it 
  svglite(file=image_name,
      width=40, 
      height=30, 
      pointsize=50)
  heatmap.2(dendrogram = "none",within_repeat_unit,cexRow=1.1,cexCol=1.1, Rowv = NA,Colv = NA,scale = "row", key = T, key.xtickfun = NULL,density.info = "none",trace = "none",symkey = F,col=redblue,na.color = "red", colsep = c(1,2,3,4,5,6),sepwidth = c(0.021, 0.021),srtCol=90, offsetRow=0, offsetCol=1,rowsep = seq(1,dim(occurrences)[1],1),xlab = paste("Codons in",OoI),ylab = paste("Codons in",other_OoI),cellnote = within_repeat_unit,notecol="black",notecex=1.5, main = paste("Codon changes within ",repeat_unit," in all ",RoI,"s",sep = ""),keysize = 1)
  dev.off()
  setwd(wd)
  # Cohen's kappa estimation ------------------------------------------------
  if(any(is.na(within_repeat_unit))){
    within_repeat_unit[which(is.na(within_repeat_unit))] = 0
    to_kappa_est = within_repeat_unit
  } else{
    to_kappa_est = within_repeat_unit
  }
  kappa_result = Kappa.test(to_kappa_est)
  string = paste(sep = ",",OoI, other_OoI,RoI,round(kappa_result$Result$estimate,digits = 2),paste(sep = "",round(kappa_result$Result$conf.int[1],digits = 2)," - ",round(kappa_result$Result$conf.int[2],digits = 2)),kappa_result$Result$p.value)
  if(file.exists(paste(dir,"/visualization_",OoI,"/kappa_results_",other_OoI,".csv",sep = ""))){
    write(string,file = paste(dir,"/visualization_",OoI,"/kappa_results_",other_OoI,".csv",sep = ""),append = T);
  } else{
    header = paste(sep=",","higher_organism","lower_organism","region","Kappa estimate","CI","p-value")
    write(header,file = paste(dir,"/visualization_",OoI,"/kappa_results_",other_OoI,".csv",sep = ""),append = T);
    write(string,file = paste(dir,"/visualization_",OoI,"/kappa_results_",other_OoI,".csv",sep = ""),append = T);
  }
  i=i+1
}


 