#USAGE: RScprit name organism_of_inerest_name organism_of_interest_exclusive_code

wd = getwd()
require(seqinr)
require(Biostrings)
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3){
  print("Wrong number of arguments passed!")
  stop()
}

organism_of_interest_name = args[1]
organisms = c("bos_taurus","gallus_gallus","gorilla_gorilla","homo_sapiens","macaca_mulatta","mus_musculus","pan_troglodytes","rattus_norvegicus","xenopus_tropicalis")
tmp_organisms = organisms[-which(prot_names == organism_of_interest_name)]
ensembl_exclusive_codes = c("ENSBTAP","ENSGALP","ENSGGOP","ENSP","ENSMMUP","ENSMUSP","ENSPTRP","ENSRNOP","ENSXETP");
tmp_codes = ensembl_exclusive_codes[-which(ensembl_exclusive_codes == args[2])]

# Reading the data (output from run_revtrans.pl) --------------------------
organism_of_interest_data = as.matrix(read.csv(file = paste(organism_of_interest_name,"_revtrans.csv",sep = ""),header = F))
organism_of_interest = lapply(1:9,function(x) matrix(NaN,nrow = dim(organism_of_interest_data)[1],ncol = dim(organism_of_interest_data)[2]));
for (org in seq(1,length(tmp_codes),by = 1)){
  reg_ex = paste(tmp_codes[org],"[0-9]+",sep = "");
  organism_of_interest[[org]] = organism_of_interest_data[grep(pattern = reg_ex ,organism_of_interest_data[,3],perl = T),]
}
names(organism_of_interest) = tmp_organisms

# Loop constructing thee lists with codons on leucine codons positions in organism of interest cDNA sequences --------
OoI = lapply(1:length(v), function(x) lapply(1:length(leucine_codons), function(x) NaN))
other_codons = c();
for (org in seq(1,length(v),by = 1)){
  print(paste("Processing orgsnism: ",sep = "", current_names[org]))
  for (l in seq(1, length(leucine_codons),by = 1)){
    print(paste("Codon: ",sep = "",leucine_codons[l]))
    for (i in seq(1, dim(v[[org]])[1],by = 1)){
      organism_of_interest_codons = strsplit(organism_of_interest[[org]][i,2], "(?<=.{3})", perl = TRUE)[[1]]
      organism_of_interest_L_indices = which(organism_of_interest_codons == leucine_codons[l])
      other_codons_tmp = strsplit(organism_of_interest[[org]][i,4], "(?<=.{3})", perl = TRUE)[[1]][organism_of_interestL_indices]
      other_codons = append(other_codons,other_codons_tmp)
    }
    OoI[[org]][[l]] = table(other_codons)
    other_codons = c();
  }
}
names(OoI) = current_names
for (org in seq(1,length(OoI),by = 1)){
  names(OoI[[org]]) = leucine_codons
}

# Saving codons of other organisms into matrix - NEEDS CORRECTION ---------------------------


output = c()
for (org in seq(1,length(OoI),by = 1)){
  for (l in seq(1,length(leucine_codons),by = 1)){
    row.names(output) = output[,1]
    output = merge(output,as.matrix(OoI[[org]][[l]]),by="row.names",all=T)
  }
  row_names = output[,1]
  output = output[,-c(1:length(leucine_codons))]
  rownames(output) = row_names
  colnames(output) = leucine_codons
}
