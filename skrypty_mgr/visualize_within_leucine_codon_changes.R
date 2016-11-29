require(gplots)
require(seqinr)
require(Biostrings)
# USAGE Rscript visualize_within_leucine_codon_changes.R <higher organism of interest(full latin name)> <lower organism of interest(full latin name)>
# parameters --------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2){
  print("Wrong number of arguments passed!")
  stop()
}

OoI = args[1] # higher organism of interest
lower_OoI = args[2] # lower organism of interest
all_codon_names = append("---",sort(names(GENETIC_CODE)))
leucine_codons = c("TTA","TTG","CTT","CTC","CTA","CTG");
within_leucine_list = list()
regions = c("L","LSAAR")
i=1
# comparison --------------------------------------------------------------
for (RoI in regions){
  occurrences=data.matrix(read.csv(paste("/home/mstolarczyk/Uczelnia/MGR/praca_magisterska/codon_analysis/within_leucine/",OoI,"/codon_changes_within_",RoI,"_",lower_OoI,".csv",sep = ""),header = T,sep = ",",row.names = all_codon_names))[,-1]
  occurrences[which(occurrences==0)] = NA
  within_leucine = occurrences[leucine_codons,]
  within_leucine_list[[i]] = within_leucine
  svg(filename=paste("/home/mstolarczyk/Uczelnia/MGR/praca_magisterska/codon_analysis/within_leucine/",OoI,"/codon_changes_within_",RoI,"_",lower_OoI,".svg",sep = ""), 
      width=40, 
      height=30, 
      pointsize=50)
  heatmap.2(dendrogram = "none",within_leucine,cexRow=1.1,cexCol=1.1, Rowv = NA,Colv = NA,scale = "row", key = T, key.xtickfun = NULL,density.info = "none",trace = "none",symkey = F,col=redblue,na.color = "red", colsep = c(1,2,3,4,5,6),sepwidth = c(0.021, 0.021),srtCol=90, offsetRow=0, offsetCol=1,rowsep = seq(1,dim(occurrences)[1],1),xlab = paste("Codons in",OoI),ylab = paste("Codons in",lower_OoI),cellnote = within_leucine,notecol="black",notecex=1.5, main = paste("Codon changes within leucine in all ",RoI,"s",sep = ""),keysize = 1)
  dev.off()
  
  # slippage vs. point mutations --------------------------------------------
  
  output = matrix(data=NA,nrow = length(leucine_codons),ncol = 2)
  counter = 1
  for(codon in leucine_codons){
    to_consider = occurrences[-which(rownames(occurrences)==codon),]
    slippage = to_consider[which(rownames(to_consider)=="---"),which(colnames(to_consider) == codon)]
    point_mutations = sum(to_consider[-which(rownames(to_consider)=="---"),which(colnames(to_consider) == codon)],na.rm = T)
    output[counter,1] = point_mutations
    output[counter,2] = slippage
    counter = counter + 1
  }
  rownames(output) = leucine_codons
  colnames(output) = c("point_mutations","slippage")
  write.csv(x = output,file = paste(sep = "","/home/mstolarczyk/Uczelnia/MGR/praca_magisterska/codon_analysis/to_leucine_transformation/",OoI,"/",lower_OoI,"_",OoI,"_codon_changes_in_",RoI,".csv"))
  print(paste(sep = "","Output file written in: ","/home/mstolarczyk/Uczelnia/MGR/praca_magisterska/codon_analysis/to_leucine_transformation/",OoI,"/",lower_OoI,"_",OoI,"_codon_changes_in_",RoI,".csv"))
  i=i+1
}

within_leucine_no_LSAAR = within_leucine_list[[1]] - within_leucine_list[[2]]
within_leucine_no_LSAAR[which(within_leucine_no_LSAAR==0)] = NA
svg(filename=paste("/home/mstolarczyk/Uczelnia/MGR/praca_magisterska/codon_analysis/within_leucine/",OoI,"/codon_changes_within_","L\\LSAAR","_",lower_OoI,".svg",sep = ""), 
    width=40, 
    height=30, 
    pointsize=50)
heatmap.2(dendrogram = "none",within_leucine_no_LSAAR,cexRow=1.1,cexCol=1.1, Rowv = NA,Colv = NA,scale = "row", key = T, key.xtickfun = NULL,density.info = "none",trace = "none",symkey = F,col=redblue,na.color = "red", colsep = c(1,2,3,4,5,6),sepwidth = c(0.021, 0.021),srtCol=90, offsetRow=0, offsetCol=1,rowsep = seq(1,dim(occurrences)[1],1),xlab = paste("Codons in",OoI),ylab = paste("Codons in",lower_OoI),cellnote = within_leucine_no_LSAAR,notecol="black",notecex=1.5, main = paste("Codon changes within leucine in ","L\\LSAAR","s",sep = ""),keysize = 1)
dev.off()
 