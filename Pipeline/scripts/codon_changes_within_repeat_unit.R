#USAGE: RScprit name organism_of_inerest_name organism_of_interest_exclusive_code
require(seqinr)
require(Biostrings)

options(warn=-1)
wd = getwd()
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3){
  print("Wrong number of arguments passed!")
  stop()
}

unit_codons = names(which(GENETIC_CODE == repeat_unit))

revtrans_file = args[1]
codes_dict=read.csv(args[2],stringsAsFactors = F)
organism_of_interest_name = args[3]
repeat_unit = args[4]
organis_of_interest_code=as.character(codes_dict[organis_of_interest_name])
organisms = colnames(codes_dict)
tmp_organisms = organisms[-which(organisms == organism_of_interest_name)]
ensembl_exclusive_codes = as.character(codes_dict)
tmp_codes = ensembl_exclusive_codes[-which(ensembl_exclusive_codes == organis_of_interest_code)]

# Reading the data (output from run_revtrans.pl) --------------------------
organism_of_interest_data = as.matrix(read.csv(file = revtrans_file,header = F))
organism_of_interest = lapply(1:length(tmp_codes),function(x) matrix(NaN,nrow = dim(organism_of_interest_data)[1],ncol = dim(organism_of_interest_data)[2]));
for (org in seq(1,length(tmp_codes),by = 1)){
  reg_ex = paste(tmp_codes[org],"[0-9]+",sep = "");
  organism_of_interest[[org]] = organism_of_interest_data[grep(pattern = reg_ex ,organism_of_interest_data[,3],perl = T),]
}
names(organism_of_interest) = tmp_organisms

# Loop constructing thee lists with codons on leucine codons positions in organism of interest cDNA sequences --------
OoI = lapply(1:length(organism_of_interest), function(x) lapply(1:length(unit_codons), function(x) NaN))
OoI_SAAR = lapply(1:length(organism_of_interest), function(x) lapply(1:length(unit_codons), function(x) NaN))
other_codons = c();
other_codons_SAAR = c();
for (org in seq(1,length(tmp_codes),by = 1)){
  print(paste("Processing orgsnism: ",sep = "", tmp_organisms[org]))
  for (l in seq(1, length(unit_codons),by = 1)){
    for (i in seq(1, dim(organism_of_interest[[org]])[1],by = 1)){
      check_for_ns = s2c(organism_of_interest[[org]][i,2]);
      check_for_ns[which(check_for_ns == "N")] = "C"
      organism_of_interest[[org]][i,2] = c2s(check_for_ns)
      organism_of_interest_codons = strsplit(organism_of_interest[[org]][i,2], "(?<=.{3})", perl = TRUE)[[1]]
      organism_of_interest_repeat_unit_indices = which(organism_of_interest_codons == unit_codons[l])
      organism_of_interest_codons[which(organism_of_interest_codons == "---")] = "TGG"
      aas = as.vector(translate(DNAStringSet(organism_of_interest_codons)))
      index=grepRaw(c2s(rep(repeat_unit,5)),(c2s(aas)))
      if (length(index) >= 1){
        a=rle(aas)
        length=a$lengths[which(a$values==repeat_unit)][which.max(a$lengths[which(a$values==repeat_unit)])]#Lrun length
        organism_of_interest_SAAR_indices = seq(from = index,to = (index+length-1),by = 1)
        other_codons_SAAR_tmp1 = strsplit(organism_of_interest[[org]][i,2], "(?<=.{3})", perl = TRUE)[[1]][organism_of_interest_SAAR_indices]
        which_lcodon_analyzed = which(other_codons_SAAR_tmp1 == unit_codons[l])
        proper_indices = (index + which_lcodon_analyzed)  - 1
        other_codons_SAAR_tmp = strsplit(organism_of_interest[[org]][i,4], "(?<=.{3})", perl = TRUE)[[1]][proper_indices]
        other_codons_SAAR = append(other_codons_SAAR,other_codons_SAAR_tmp)
      }
      other_codons_tmp = strsplit(organism_of_interest[[org]][i,4], "(?<=.{3})", perl = TRUE)[[1]][organism_of_interest_repeat_unit_indices]
      other_codons = append(other_codons,other_codons_tmp)
    }
    OoI[[org]][[l]] = sort(table(other_codons),decreasing = T)
    other_codons = c();
    OoI_SAAR[[org]][[l]] = sort(table(other_codons_SAAR),decreasing = T)
    other_codons_SAAR = c();
  }
}
names(OoI) = tmp_organisms
for (org in seq(1,length(OoI),by = 1)){
  names(OoI[[org]]) = unit_codons
}

names(OoI_SAAR) = tmp_organisms
for (org in seq(1,length(OoI_SAAR),by = 1)){
  names(OoI_SAAR[[org]]) = unit_codons
}

# results saving -------------------------------------------------------------
all_codon_names = append("---",sort(names(GENETIC_CODE)))
nms = names(OoI)
for (org in seq(1,length(OoI),by = 1)){
  output=matrix(data=NA,nrow = length(all_codon_names),ncol = 1)
  rownames(output) = all_codon_names
  for (l in seq(1,length(unit_codons),by = 1)){
    tmp = as.matrix(OoI[[org]][[l]])
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
  colnames(output) = unit_codons
  write.csv(x = output,file = paste(wd,"/",organism_of_interest_name,"_changes_within_repeatUnit/codon_changes_within_repeat_unit_",nms[org],".csv",sep = ""),row.names = T)
}
print(paste("Results saved to:", paste(wd,"/",organism_of_interest_name,"_changes_within_repeatUnit/codon_changes_within_repeat_unit_<organism name>")))


# results saving L-SAARs --------------------------------------------------

nms = names(OoI_SAAR)
for (org in seq(1,length(OoI_SAAR),by = 1)){
  output=matrix(data=NA,nrow = length(all_codon_names),ncol = 1)
  rownames(output) = all_codon_names
  for (l in seq(1,length(unit_codons),by = 1)){
    tmp = as.matrix(OoI_SAAR[[org]][[l]])
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
  colnames(output) = unit_codons
  write.csv(x = output,file = paste(wd,"/",organism_of_interest_name,"_changes_within_repeatUnit/codon_changes_within_SAAR_",nms[org],".csv",sep = ""),row.names = T)
}

