library(Biostrings);library(seqinr)
setwd(dir = "/home/mstolarczyk/Uczelnia/MGR/praca_magisterska/analiza_kodonow/within_leucine/")
homo_sapiens_data = as.matrix(read.csv(file = "homo_sapiens_revtrans.csv",header = F))
a = matrix(NaN,nrow = dim(homo_sapiens_data)[1],ncol = dim(homo_sapiens_data)[2])
homo_sapiens = list(a,a,a,a,a,a,a,a);
homo_sapiens[[6]] = homo_sapiens_data[grep(pattern = "ENSPTRP[0-9]+",homo_sapiens_data[,3],perl = T),]
homo_sapiens[[4]] = homo_sapiens_data[grep(pattern = "ENSMMUP[0-9]+",homo_sapiens_data[,3],perl = T),]
homo_sapiens[[3]] = homo_sapiens_data[grep(pattern = "ENSGGOP[0-9]+",homo_sapiens_data[,3],perl = T),]
homo_sapiens[[5]] = homo_sapiens_data[grep(pattern = "ENSMUSP[0-9]+",homo_sapiens_data[,3],perl = T),]
homo_sapiens[[7]] = homo_sapiens_data[grep(pattern = "ENSRNOP[0-9]+",homo_sapiens_data[,3],perl = T),]
homo_sapiens[[1]] = homo_sapiens_data[grep(pattern = "ENSBTAP[0-9]+",homo_sapiens_data[,3],perl = T),]
homo_sapiens[[2]] = homo_sapiens_data[grep(pattern = "ENSGALP[0-9]+",homo_sapiens_data[,3],perl = T),]
homo_sapiens[[8]] = homo_sapiens_data[grep(pattern = "ENSXETP[0-9]+",homo_sapiens_data[,3],perl = T),]

names(homo_sapiens) = prot_names[-which(prot_names == "homo_sapiens")]


for (org in seq(1,length(homo_sapiens),by = 1)){
  for (i in seq(1, dim(homo_sapiens)[1],by = 1)){
    for (l in seq(1, length(leucine_codons),by = 1)){
      human_codons = strsplit(homo_sapiens[i,2], "(?<=.{3})", perl = TRUE)[[1]]
      human_L_indices = which(human_codons == leucine_codons[l])
      other_codons = strsplit(homo_sapiens[i,4], "(?<=.{3})", perl = TRUE)[[1]][human_L_indices]
    }
  }
}
