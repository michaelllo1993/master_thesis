library(Biostrings);library(seqinr)
load("compare_region/data.RData")
load("compare_region/with_protein.RData")
load("compare_region/with_sigp.RData")
# all ---------------------------------------------------------------------
a = matrix(data = NaN,nrow = length(cDNA$cDNA_homo_sapiens$SEQUENCE),ncol = length(leucine_codons));
codons_count_protein = list(a,a,a,a,a,a,a,a,a);
leucine_codons = c("TTA","TTG","CTT","CTC","CTA","CTG");
for(org in seq(1,9,1)){
  for (i in seq(1,length(cDNA[[org]]$SEQUENCE))) {
    if (!length(sigp[[org]]$cDNA_ID[which(sigp[[org]]$cDNA_ID == cDNA[[org]]$cDNA_ID[i])]) == 0){
      sp_end = sigp[[org]]$SP_END[which(sigp[[org]]$cDNA_ID == cDNA[[org]]$cDNA_ID[i])]
    }
    else {
      sp_end = 1;
    }
    for (codon in seq(1,length(leucine_codons),by = 1)){
      seq=c2s(s2c(cDNA[[org]]$SEQUENCE[i])[seq(sp_end*3,length(s2c(cDNA[[org]]$SEQUENCE[i])),by=1)])
      codons_count_protein[[org]][i,codon]=length(which(strsplit(seq, "(?<=.{3})", perl = TRUE)[[1]] == leucine_codons[codon]))
    }
  }
  print(paste("ORGANISM",org,"done"))
}
# sigp --------------------------------------------------------------------
a = matrix(data = NaN,nrow = length(sigp$sigp_ALL_homo_sapiens$SEQ),ncol = length(leucine_codons));
codons_count_sigp = list(a,a,a,a,a,a,a,a,a);
leucine_codons = c("TTA","TTG","CTT","CTC","CTA","CTG");
for(org in seq(1,9,1)) {
  for (i in seq(1,length(sigp[[org]]$SEQ))) {
    if (!length(sigpL[[org]]$cDNA_ID[which(sigpL[[org]]$cDNA_ID == sigp[[org]]$cDNA_ID[i])]) == 0){
      sp_end = sigpL[[org]]$SP_END[which(sigpL[[org]]$cDNA_ID == sigp[[org]]$cDNA_ID[i])]
      lsaar_start = sigpL[[org]]$LSAAR_START[which(sigpL[[org]]$cDNA_ID == sigp[[org]]$cDNA_ID[i])]
      lsaar_length = sigpL[[org]]$LSAAR_LENGTH[which(sigpL[[org]]$cDNA_ID == sigp[[org]]$cDNA_ID[i])]
      for (codon in seq(1,length(leucine_codons),by = 1)) {
        seq=c2s(s2c(cDNA[[org]]$SEQ[which(cDNA[[org]]$cDNA_ID == sigp[[org]]$cDNA_ID[i])])[c(seq(1,(lsaar_start*3+1),1),seq((lsaar_start*3+1)+(lsaar_length*3),length(s2c(sigp[[org]]$SEQ[i]))*3))])
        codons_count_sigp[[org]][i,codon]=length(which(strsplit(seq, "(?<=.{3})", perl = TRUE)[[1]] == leucine_codons[codon]))
      }
    }
    else {
      for (codon in seq(1,length(leucine_codons),by = 1)) {
        seq=c2s(s2c(cDNA[[org]]$SEQ[which(cDNA[[org]]$cDNA_ID == sigp[[org]]$cDNA_ID[i])])[seq(1,sigp[[org]]$SP_END[i]*3)])
        codons_count_sigp[[org]][i,codon]=length(which(strsplit(seq, "(?<=.{3})", perl = TRUE)[[1]] == leucine_codons[codon]))
      }
  }
  }
  print(paste("ORGANISM",org,"done"))
}
# L-SAAR ------------------------------------------------------------------
a = matrix(data = NaN,nrow = length(sigpL$sigpL_ALL_homo_sapiens$SEQ),ncol = length(leucine_codons));
codons_count_lsaar = list(a,a,a,a,a,a,a,a,a);
leucine_codons = c("TTA","TTG","CTT","CTC","CTA","CTG");
for(org in seq(1,9,1)) {
  for (i in seq(1,length(sigpL[[org]]$SEQ))) {
      for (codon in seq(1,length(leucine_codons),by = 1)) {
        if(!is.na(sigpL[[org]]$LSAAR_LENGTH[i])){
        seq=c2s(s2c(cDNA[[org]]$SEQ[which(cDNA[[org]]$cDNA_ID == sigpL[[org]]$cDNA_ID[i])])[seq((sigpL[[org]]$LSAAR_START[i]*3+1),(sigpL[[org]]$LSAAR_START[i]*3+1)+(sigpL[[org]]$LSAAR_LENGTH[i]*3))])
        codons_count_lsaar[[org]][i,codon]=length(which(strsplit(seq, "(?<=.{3})", perl = TRUE)[[1]] == leucine_codons[codon]))
        }
      }
  }
  print(paste("ORGANISM",org,"done"))
}
# data manipualtion -------------------------------------------------------
for (i in seq(1,length(codons_count_lsaar),1)){
  codons_count_lsaar[[i]] = codons_count_lsaar[[i]][seq(1,length(which(!is.na(codons_count_lsaar[[i]][,1]))),1),]
  codons_count_protein[[i]] = codons_count_protein[[i]][seq(1,length(which(!is.na(codons_count_protein[[i]][,1]))),1),]
  codons_count_sigp[[i]] = codons_count_sigp[[i]][seq(1,length(which(!is.na(codons_count_sigp[[i]][,1]))),1),]
  colnames(codons_count_lsaar[[i]]) = leucine_codons
  colnames(codons_count_sigp[[i]]) = leucine_codons
  colnames(codons_count_protein[[i]]) = leucine_codons
}
names(codons_count_lsaar)=prot_names
names(codons_count_protein)=prot_names
names(codons_count_sigp)=prot_names
# data into matrices ------------------------------------------------------
codons_count_lsaar_fraction = matrix(NaN,nrow = length(codons_count_lsaar) ,ncol = length(leucine_codons))
for (org in seq(1,length(codons_count_lsaar),1)){
  for (cod in seq(1,length(leucine_codons),1)){
    codons_count_lsaar_fraction[org,cod] = sum(codons_count_lsaar[[org]][,cod],na.rm = T)/sum(codons_count_lsaar[[org]],na.rm = T)
  }
}
codons_count_sigp_fraction = matrix(NaN,nrow = length(codons_count_sigp) ,ncol = length(leucine_codons))
for (org in seq(1,length(codons_count_sigp),1)){
  for (cod in seq(1,length(leucine_codons),1)){
    codons_count_sigp_fraction[org,cod] = sum(codons_count_sigp[[org]][,cod],na.rm = T)/sum(codons_count_sigp[[org]],na.rm = T)
  }
}
codons_count_protein_fraction = matrix(NaN,nrow = length(codons_count_protein) ,ncol = length(leucine_codons))
for (org in seq(1,length(codons_count_protein),1)) {
  for (cod in seq(1,length(leucine_codons),1)) {
    codons_count_protein_fraction[org,cod] = sum(codons_count_protein[[org]][,cod],na.rm = T)/sum(codons_count_protein[[org]],na.rm = T)
  }
}
rownames(codons_count_protein_fraction)=prot_names
rownames(codons_count_sigp_fraction)=prot_names
rownames(codons_count_lsaar_fraction)=prot_names
colnames(codons_count_protein_fraction)=leucine_codons
colnames(codons_count_sigp_fraction)=leucine_codons
colnames(codons_count_lsaar_fraction)=leucine_codons


# Statistical analysis of CTG codon usage in different regions ---------------------------------------
# Normality test
# all
normality_test_CTG_lsaars_all = shapiro.test(x = codons_count_lsaar_fraction[,6])
normality_test_CTG_sigp_all = shapiro.test(x = codons_count_sigp_fraction[,6])
normality_test_CTG_protein_all = shapiro.test(x = codons_count_protein_fraction[,6])

qqnorm(codons_count_lsaar_fraction[,6])
qqline(codons_count_lsaar_fraction[,6],col="red")
qqnorm(codons_count_sigp_fraction[,6])
qqline(codons_count_sigp_fraction[,6],col="red")
qqnorm(codons_count_protein_fraction[,6])
qqline(codons_count_protein_fraction[,6],col="red")

# w/o frog
normality_test_CTG_lsaars_nofrog = shapiro.test(x = codons_count_lsaar_fraction[1:8,6])
normality_test_CTG_sigp_nofrog = shapiro.test(x = codons_count_sigp_fraction[1:8,6])
normality_test_CTG_protein_nofrog = shapiro.test(x = codons_count_protein_fraction[1:8,6])

qqnorm(codons_count_lsaar_fraction[1:8,6])
qqline(codons_count_lsaar_fraction[1:8,6],col="red")
qqnorm(codons_count_sigp_fraction[1:8,6])
qqline(codons_count_sigp_fraction[1:8,6],col="red")
qqnorm(codons_count_protein_fraction[1:8,6])
qqline(codons_count_protein_fraction[1:8,6],col="red")

# variance equality test
var_equality_lsaar_sigp = var.test(codons_count_lsaar_fraction[1:8,6],codons_count_sigp_fraction[1:8,6])
var_equality_sigp_protein = var.test(codons_count_sigp_fraction[1:8,6],codons_count_protein_fraction[1:8,6])

# mean equality test
t_test_lsaar_sigp = t.test(x = codons_count_lsaar_fraction[1:9,6], y = codons_count_sigp_fraction[1:8,6], alternative = "greater",var.equal = T)
t_test_sigp_protein = t.test(x = codons_count_sigp_fraction[1:9,6], y = codons_count_protein_fraction[1:8,6], alternative = "greater",var.equal = F)

# saving environment ------------------------------------------------------
save.image(file = paste(sep = "",date(),".RData"))