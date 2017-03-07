library(gplots)

# parameters --------------------------------------------------------------

OoI = "homo_sapiens" # higher organism of interest
lower_OoI = "xenopus_tropicalis" # lower organism of interest
RoI = "L" # region of interest

# comparison --------------------------------------------------------------

occurrences=data.matrix(read.csv(paste("/home/mstolarczyk/Uczelnia/MGR/praca_magisterska/analiza_kodonow/within_leucine/",OoI,"/codon_changes_within_",RoI,"_",lower_OoI,".csv",sep = ""),header = T,sep = ",",row.names = all_codon_names))[,-1]
occurrences[which(occurrences==0)] = NA
within_leucine = occurrences[leucine_codons,]
heatmap.2(within_leucine, Rowv = NA,Colv = NA,scale = "row", key = T, key.xtickfun = NULL,density.info = "none",trace = "none",symkey = F,col=redblue,na.color = "red", colsep = c(1,2,3,4,5,6),sepwidth = c(0.001, 0.001),rowsep = seq(1,dim(occurrences)[1],1),xlab = paste("Codons in",OoI),ylab = paste("Codons in",lower_OoI),cellnote = within_leucine,notecol="black",notecex=1.5, main = paste("Codon changes within leucine in all ",RoI,"s",sep = ""),keysize = 1)

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
write.csv(x = output,file = paste(sep = "","/home/mstolarczyk/Uczelnia/MGR/praca_magisterska/analiza_kodonow/to_leucine_transformation/",OoI,"/",lower_OoI,"_",OoI,"_codon_changes_in_",RoI,".csv"))
