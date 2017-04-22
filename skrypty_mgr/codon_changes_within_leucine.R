#USAGE: RScprit name organism_of_inerest_name organism_of_interest_exclusive_code
options(warn=-1)
wd = getwd()
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2){
  print("Wrong number of arguments passed!")
  stop()
}
require(seqinr)
require(Biostrings)

organism_of_interest_name = args[1]
organisms = c("bos_taurus","gallus_gallus","gorilla_gorilla","homo_sapiens","macaca_mulatta","mus_musculus","pan_troglodytes","rattus_norvegicus","xenopus_tropicalis")
tmp_organisms = organisms[-which(organisms == organism_of_interest_name)]
ensembl_exclusive_codes = c("ENSBTAP","ENSGALP","ENSGGOP","ENSP","ENSMMUP","ENSMUSP","ENSPTRP","ENSRNOP","ENSXETP");
tmp_codes = ensembl_exclusive_codes[-which(ensembl_exclusive_codes == args[2])]
leucine_codons = c("TTA","TTG","CTT","CTC","CTA","CTG");

# Reading the data (output from run_revtrans.pl) --------------------------
organism_of_interest_data = as.matrix(read.csv(file = paste(organism_of_interest_name,"_revtrans.csv",sep = ""),header = F))
organism_of_interest = lapply(1:length(tmp_codes),function(x) matrix(NaN,nrow = dim(organism_of_interest_data)[1],ncol = dim(organism_of_interest_data)[2]));
for (org in seq(1,length(tmp_codes),by = 1)){
  reg_ex = paste(tmp_codes[org],"[0-9]+",sep = "");
  organism_of_interest[[org]] = organism_of_interest_data[grep(pattern = reg_ex ,organism_of_interest_data[,3],perl = T),]
}
names(organism_of_interest) = tmp_organisms

# Loop constructing thee lists with codons on leucine codons positions in organism of interest cDNA sequences --------
OoI = lapply(1:length(organism_of_interest), function(x) lapply(1:length(leucine_codons), function(x) NaN))
OoI_LSAAR = lapply(1:length(organism_of_interest), function(x) lapply(1:length(leucine_codons), function(x) NaN))
other_codons = c();
other_codons_LSAAR = c();
for (org in seq(1,length(tmp_codes),by = 1)){
  print(paste("Processing orgsnism: ",sep = "", tmp_organisms[org]))
  for (l in seq(1, length(leucine_codons),by = 1)){
    for (i in seq(1, dim(organism_of_interest[[org]])[1],by = 1)){
      check_for_ns = s2c(organism_of_interest[[org]][i,2]);
      check_for_ns[which(check_for_ns == "N")] = "C"
      organism_of_interest[[org]][i,2] = c2s(check_for_ns)
      organism_of_interest_codons = strsplit(organism_of_interest[[org]][i,2], "(?<=.{3})", perl = TRUE)[[1]]
      organism_of_interest_L_indices = which(organism_of_interest_codons == leucine_codons[l])
      organism_of_interest_codons[which(organism_of_interest_codons == "---")] = "TGG"
      aas = as.vector(translate(DNAStringSet(organism_of_interest_codons)))
      index=grepRaw(c2s(rep("L",5)),(c2s(aas)))
      if (length(index) >= 1){
        a=rle(aas)
        length=a$lengths[which(a$values=="L")][which.max(a$lengths[which(a$values=="L")])]#Lrun length
        organism_of_interest_LSAAR_indices = seq(from = index,to = (index+length-1),by = 1)
        other_codons_LSAAR_tmp1 = strsplit(organism_of_interest[[org]][i,2], "(?<=.{3})", perl = TRUE)[[1]][organism_of_interest_LSAAR_indices]
        which_lcodon_analyzed = which(other_codons_LSAAR_tmp1 == leucine_codons[l])
        proper_indices = (index + which_lcodon_analyzed)  - 1
        other_codons_LSAAR_tmp = strsplit(organism_of_interest[[org]][i,4], "(?<=.{3})", perl = TRUE)[[1]][proper_indices]
        other_codons_LSAAR = append(other_codons_LSAAR,other_codons_LSAAR_tmp)
      }
      other_codons_tmp = strsplit(organism_of_interest[[org]][i,4], "(?<=.{3})", perl = TRUE)[[1]][organism_of_interest_L_indices]
      other_codons = append(other_codons,other_codons_tmp)
    }
    OoI[[org]][[l]] = sort(table(other_codons),decreasing = T)
    names(OoI[[org]][[l]]) = toupper(names(OoI[[org]][[l]]))
    other_codons = c();
    OoI_LSAAR[[org]][[l]] = sort(table(other_codons_LSAAR),decreasing = T)
    names(OoI_LSAAR[[org]][[l]]) = toupper(names(OoI_LSAAR[[org]][[l]]))
    other_codons_LSAAR = c();
  }
}
names(OoI) = tmp_organisms
for (org in seq(1,length(OoI),by = 1)){
  names(OoI[[org]]) = leucine_codons
}

names(OoI_LSAAR) = tmp_organisms
for (org in seq(1,length(OoI_LSAAR),by = 1)){
  names(OoI_LSAAR[[org]]) = leucine_codons
}

# results saving -------------------------------------------------------------
all_codon_names = append("---",sort(names(GENETIC_CODE)))
nms = names(OoI)
for (org in seq(1,length(OoI),by = 1)){
  output=matrix(data=NA,nrow = length(all_codon_names),ncol = 1)
  rownames(output) = all_codon_names
  for (l in seq(1,length(leucine_codons),by = 1)){
    #tmp = as.matrix(OoI[[org]][[l]])
    tmp = OoI[[org]][[l]]
    if(any(grep("N",rownames(tmp),perl = T))){
      rows_to_delete = grep("N",rownames(tmp),perl = T)
      tmp = tmp[-rows_to_delete]
    }
    if(length(setdiff(names(tmp), all_codon_names))!=0){
      rows_to_delete = which(names(tmp)==setdiff(names(tmp), all_codon_names))
      tmp = tmp[-rows_to_delete]
    }
    output=merge(output,tmp,by="row.names",all=T)
    output=as.matrix(output[,-1])
    rownames(output) = all_codon_names
  }
  output=output[,-1]
  output[which(is.na(output))] = 0;
  colnames(output) = leucine_codons
  write.csv(x = output,file = paste(wd,"/codon_changes_within_L_",nms[org],".csv",sep = ""),row.names = T)
}
print(paste("Results saved to:", paste(wd,"/codon_changes_within_L_<organism name>")))


# results saving L-SAARs --------------------------------------------------

nms = names(OoI_LSAAR)
for (org in seq(1,length(OoI_LSAAR),by = 1)){
  output=matrix(data=NA,nrow = length(all_codon_names),ncol = 1)
  rownames(output) = all_codon_names
  for (l in seq(1,length(leucine_codons),by = 1)){
    tmp = as.matrix(OoI_LSAAR[[org]][[l]])
    if(any(grep("N",rownames(tmp),perl = T))){
      rows_to_delete = grep("N",rownames(tmp),perl = T)
      tmp = tmp[-rows_to_delete,]
    }
    output=merge(output,tmp,by="row.names",all=T)
    output=as.matrix(output[,-1])
    rownames(output) = all_codon_names
  }
  output=output[,-1]
  output[which(is.na(output))] = 0;
  colnames(output) = leucine_codons
  write.csv(x = output,file = paste(wd,"/codon_changes_within_LSAAR_",nms[org],".csv",sep = ""),row.names = T)
}
print(paste("Results saved to:", paste(wd,"/codon_changes_within_LSAAR_<organism name>",sep="")))


