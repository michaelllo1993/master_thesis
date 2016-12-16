library(Biostrings);library(seqinr);library(parallel);  library(foreach);library(doParallel)
setwd("~/Uczelnia/MGR/praca_magisterska/analiza_kodonow/");
load(file = "newest.RData");

########### Loading Data ###################

setwd("~/Uczelnia/MGR/praca_magisterska/analiza_kodonow/proteins/");
myFiles <- list.files(pattern = "*.csv");
proteins = list();
prot_names = c();
for (k in 1:length(myFiles)) {
  prot_names[k] = strsplit(myFiles,"[.]")[[k]][1];
  proteins[[k]] = read.csv(myFiles[k],stringsAsFactors = F);
}
prot_names[k] = strsplit(myFiles,"[.]")[[k]][1];
names(proteins) = prot_names;

setwd("~/Uczelnia/MGR/praca_magisterska/analiza_kodonow/mappers/");
myFiles1 <- list.files(pattern = "*.csv");
mapper = list();
map_names = c();
for (k in 1:length(myFiles1)) {
  map_names[k] = strsplit(myFiles1,"[.]")[[k]][1];
  mapper[[k]] = read.csv(myFiles1[k],stringsAsFactors = F);
}
names(mapper) = map_names;

setwd("~/Uczelnia/MGR/praca_magisterska/analiza_kodonow/cDNA/");
myFiles <- list.files(pattern = "*.csv");
cDNA = list();
cDNA_names = c();
for (k in 1:length(myFiles)) {
  cDNA_names[k] = strsplit(myFiles,"[.]")[[k]][1];
  cDNA[[k]] = read.csv(myFiles[k],stringsAsFactors = F);
}
cDNA_names[k] = strsplit(myFiles,"[.]")[[k]][1];
names(cDNA) = cDNA_names;

########### Split START and STOP values - time consuming ##########
for (org in seq(1,length(cDNA),by = 1)) {
  print(org)
  for (s in seq(1,length(cDNA[[org]][[3]]),by = 1)) {
    if (length(grep(";",cDNA[[org]][[3]][s])) == 1) {
      cDNA[[org]][[3]][s] = strsplit(cDNA[[org]][[3]][s][[1]],";");#split by semicolon
      cDNA[[org]][[4]][s] = strsplit(cDNA[[org]][[4]][s][[1]],";");#split by semicolon
    }
    else{
      cDNA[[org]][[3]][s] = cDNA[[org]][[3]][s][[1]];
      cDNA[[org]][[4]][s] = cDNA[[org]][[4]][s][[1]];
    }
  }
  
  ########### Codons to .csv - time consuming #############
  
  most_freq_L_codon = matrix("",length(mapper),length(mapper$protein_cdna_mapper_homo_sapiens$Protein_ID));
  most_freq_L_codon_names = matrix("",length(mapper),length(mapper$protein_cdna_mapper_homo_sapiens$Protein_ID));
  
  leucine_codons = c("TTA","TTG","CTT","CTC","CTA","CTG");
  
  for (org in seq(1,length(mapper),by = 1)) {
    result = c();
    result = matrix("",length(mapper[[org]]$Protein_ID) * 2,7);
    for (i in seq(1,length(mapper[[org]][[1]]),by = 1)) {
      if (i %% 1000 == 1) {
        print(paste(round(i / (
          length(mapper[[org]][[1]])
        ) * 100),"% DONE, organism: ",prot_names[org],sep = ""));
      }
      cdna = cDNA[[org]]$SEQUENCE[which(cDNA[[org]]$cDNA_ID == mapper[[org]]$cDNA_ID[i])]
      protein = proteins[[org]]$SEQUENCE[which(proteins[[org]]$ID == mapper[[org]]$Protein_ID[i])]
      startt = c();
      stopp = c();
      startt = sort(as.numeric(cDNA[[org]]$START[which(cDNA[[org]]$cDNA_ID ==
                                                         mapper[[org]]$cDNA_ID[i])][[1]]))
      stopp = sort(as.numeric(cDNA[[org]]$STOP[which(cDNA[[org]]$cDNA_ID ==
                                                       mapper[[org]]$cDNA_ID[i])][[1]]))
      vec = c();
      cdna = s2c(cdna);
      for (j in seq(1,length(stopp),by = 1)) {
        vec = append(x = vec,values = seq(startt[j],stopp[j],by = 1));
      }
      odp = cdna[vec];
      if (length(which(odp == "N")) != 0) {
        odp[which(odp == "N")] = "A"
      }
      if (!is.na(odp[length(odp)])) {
        codons_freq = table(as.character(codons(DNAString(c2s(
          odp
        )))));
        avail_L_codons_freq = codons_freq[leucine_codons];
        result[(i * 2) - 1,] = append(
          names(avail_L_codons_freq),mapper$protein_cdna_mapper_homo_sapiens$cDNA_ID[i]
        )
        result[i * 2,] = append(as.vector(avail_L_codons_freq),0)
      }
      else {
        result[(i * 2) - 1,] = append(rep("-",6),mapper$protein_cdna_mapper_homo_sapiens$cDNA_ID[i])
        result[i * 2,] = append(rep("-",6),0)
      }
    }
    write.csv(result,file = paste(sep = "",getwd(),"/",prot_names[org],"_result.csv"));
    print(paste("ORGANISM: ",org,"-> DONE.",sep = ""))
  }
  
  ########### Analysis all available ###########
  ########### ALL ###########
  setwd("~/Uczelnia/MGR/praca_magisterska/analiza_kodonow/cDNA/codons_freq_results/");
  myFiles <- list.files(pattern = "*.csv");
  leucine_codons = c("TTA","TTG","CTT","CTC","CTA","CTG");
  codon = list();
  codon_names = c();
  # Loop reading .csv files containing L codons frequencies
  for (k in 1:length(myFiles)) {
    codon_names[k] = strsplit(myFiles,"[.]")[[k]][1];
    codon[[k]] = read.csv(myFiles[k],stringsAsFactors = F);
    codon[[k]] = codon[[k]][,-1];
    names(codon[[k]]) = append(leucine_codons,"SEQ_ID");
  }
  names(codon) = codon_names;
  
  # Loop extracting the frequencies and seq_ids only
  for (org in seq(1,length(codon_names),by = 1)) {
    for (cod in seq(1,6,by = 1)) {
      codon[[org]][[cod]] = codon[[org]][[cod]][seq(2,length(codon[[org]][[cod]]),by = 2)];
    }
    codon[[org]][[7]] = codon[[org]][[7]][seq(1,length(codon[[org]][[7]]),by = 2)]
  }
  
  # Loop changing NAs to zeros, since that's what NA means in this case. Also, it converts chracters to numeric
  for (org in seq(1,length(codon_names),by = 1)) {
    for (cod in seq(1,6,by = 1)) {
      NAs = which(is.na(as.numeric(codon[[org]][[cod]])));
      codon[[org]][[cod]][NAs] = 0;
      codon[[org]][[cod]] = as.numeric(codon[[org]][[cod]]);
    }
  }
  
  mean_results = matrix("",nrow = length(cDNA_names),ncol = 6);
  sd_results = matrix("",nrow = length(cDNA_names),ncol = 6);
  sum_results = matrix(0,nrow = length(cDNA_names),ncol = 7);
  per_results = matrix(0,nrow = length(cDNA_names),ncol = 6);
  
  
  for (org in seq(1,length(codon_names),by = 1)) {
    for (cod in seq(1,6,by = 1)) {
      mean_results[org,cod] = mean(codon[[org]][[cod]]);
      sd_results[org,cod] = sd(codon[[org]][[cod]]);
      sum_results[org,cod] = sum(as.numeric(codon[[org]][[cod]]));
    }
    sum_results[org,7] = sum(as.numeric(sum_results[org,]));
  }
  
  for (org in seq(1,length(codon_names),by = 1)) {
    for (cod in seq(1,6,by = 1)) {
      per_results[org,cod] = sum_results[org,cod] / sum_results[org,7];
    }
  }
  
  rownames(mean_results) = codon_names;colnames(mean_results) = leucine_codons;
  rownames(sd_results) = codon_names;colnames(sd_results) = leucine_codons;
  rownames(sum_results) = codon_names;colnames(sum_results) = append(leucine_codons,"TOTAL");
  rownames(per_results) = codon_names;colnames(per_results) = leucine_codons;
  
  all_results = list(
    mean = mean_results,percentage = per_results,standard_deviation = sd_results,sum =
      sum_results
  )
  
  ########### SIGP ###########
  setwd("~/Uczelnia/MGR/praca_magisterska/analiza_kodonow/proteins/sigp/");
  myFiles_L <- list.files(pattern = "^sigpL_*");
  myFiles <-
    list.files(pattern = "^sigp_*");myFiles <-
    myFiles[!myFiles %in% myFiles_L];
  
  sigpL = list();
  sigpL_names = c();
  for (k in 1:length(myFiles_L)) {
    sigpL[[k]] = read.csv(myFiles_L[k],stringsAsFactors = F);
    sigpL_names[k] = strsplit(myFiles_L,"[.]")[[k]][1];
  }
  names(sigpL) = sigpL_names;
  
  sigp = list();
  sigp_names = c();
  for (k in 1:length(myFiles)) {
    sigp[[k]] = read.csv(myFiles[k],stringsAsFactors = F);
    sigp_names[k] = strsplit(myFiles,"[.]")[[k]][1];
  }
  names(sigp) = sigp_names;
  
  leucine_codons = c("TTA","TTG","CTT","CTC","CTA","CTG");
  rs = 1;
  err = 0;
  temp_avail_L_codons_freq = matrix("",nrow = length(cDNA[[org]]$ID), length(leucine_codons) +
                                      1)
  avail_L_codons_freq = list();
  for (org in seq(1,length(sigp_names),by = 1)) {
    temp_avail_L_codons_freq = matrix("",nrow = length(cDNA[[org]]$ID), length(leucine_codons) +
                                        1)
    for (s in seq(1,length(cDNA[[org]]$ID),by = 1)) {
      if (length(sigp[[org]]$SP_END[which(cDNA[[org]]$ID[s] == sigp[[org]]$PROTEIN_ID)]) !=
          0) {
        #if exists: length of signal peptide
        cDNA_corr = (s2c(cDNA[[org]]$SEQUENCE[s])[min(as.numeric(cDNA[[org]]$START[s][[1]])):max(as.numeric(cDNA[[org]]$STOP[s][[1]]))]); #cDNA - sequence of the whole peptide that contains the SP
        cDNA_corr_proper = cDNA_corr[seq(1,(length(s2c(sigp[[org]]$SEQ[which(cDNA[[org]]$ID[s] ==
                                                                               sigp[[org]]$PROTEIN_ID)]))) * 3)]; #cDNA - sequence corresponding exactly to the SP
        if (length(which(cDNA_corr_proper == "N") != 0)) {
          cDNA_corr_proper[which(cDNA_corr_proper == "N")] = "A"
          err = err + 1;
        }
        if (length(which(is.na(cDNA_corr_proper)) != 0)) {
          cDNA_corr_proper[which(is.na(cDNA_corr_proper))] = "A"
          err = err + 1;
        }
        if (rs %% 100 == 0) {
          print(paste((round(
            rs / length(sigp[[org]]$PROTEIN_ID) * 100
          )),"% DONE, organism: ",sigp_names[org],sep = ""))
        }
        codons_freq = table(as.character(codons(DNAString(
          c2s(cDNA_corr_proper)
        ))));
        temp_avail_L_codons_freq[rs,] = append(codons_freq[leucine_codons],sigp[[org]]$PROTEIN_ID[which(cDNA[[org]]$ID[s] ==
                                                                                                          sigp[[org]]$PROTEIN_ID)]);
        temp_avail_L_codons_freq[rs,][which(is.na(temp_avail_L_codons_freq[rs,]))] =
          0;
        rs = rs + 1;
      }
    }
    rs = 0;
    avail_L_codons_freq[[org]] = temp_avail_L_codons_freq;
  }
  
  names(avail_L_codons_freq) = sigp_names;
  for (i in seq(1,length(sigp_names),by = 1)) {
    colnames(avail_L_codons_freq[[i]]) = append(leucine_codons,"PROTEIN_ID");
  }
  
  setwd(
    "~/Uczelnia/MGR/praca_magisterska/analiza_kodonow/cDNA/codons_freq_results/sigp/"
  );
  for (org in seq(1,length(prot_names),by = 1)) {
    write.csv(
      avail_L_codons_freq[[org]],file = paste(sep = "",getwd(),"/",prot_names[org],"_sigp_result.csv")
    );
  }
  
  mean_sigp_results = matrix("",nrow = length(cDNA_names),ncol = 6);
  sd_sigp_results = matrix("",nrow = length(cDNA_names),ncol = 6);
  sum_sigp_results = matrix(0,nrow = length(cDNA_names),ncol = 7);
  per_sigp_results = matrix(0,nrow = length(cDNA_names),ncol = 6);
  
  
  for (org in seq(1,length(codon_names),by = 1)) {
    for (cod in seq(1,6,by = 1)) {
      mean_sigp_results[org,cod] = mean(as.numeric(avail_L_codons_freq[[org]][,cod])[which(!is.na(as.numeric(avail_L_codons_freq[[org]][,cod])))]);
      sd_sigp_results[org,cod] = sd(as.numeric(avail_L_codons_freq[[org]][,cod])[which(!is.na(as.numeric(avail_L_codons_freq[[org]][,cod])))]);
      sum_sigp_results[org,cod] = sum(as.numeric(avail_L_codons_freq[[org]][,cod])[which(!is.na(as.numeric(avail_L_codons_freq[[org]][,cod])))]);
    }
    sum_sigp_results[org,7] = sum(as.numeric(sum_sigp_results[org,]));
  }
  
  
  rownames(mean_sigp_results) = codon_names;colnames(mean_sigp_results) =
    leucine_codons;
  rownames(sd_sigp_results) = codon_names;colnames(sd_sigp_results) = leucine_codons;
  rownames(sum_sigp_results) = codon_names;colnames(sum_sigp_results) = append(leucine_codons,"suma");
  rownames(per_sigp_results) = codon_names;colnames(per_sigp_results) = leucine_codons;
  
  for (org in seq(1,length(codon_names),by = 1)) {
    for (cod in seq(1,6,by = 1)) {
      per_sigp_results[org,cod] = sum_sigp_results[org,cod] / sum_sigp_results[org,7];
    }
  }
  
  sigp_results = list(
    mean = mean_sigp_results,percentage = per_sigp_results,standard_deviation =
      sd_sigp_results,sum = sum_sigp_results
  )
  
  setwd(
    "~/Uczelnia/MGR/praca_magisterska/analiza_kodonow/cDNA/codons_freq_results/statresults"
  );
  for (i in seq(1,length(sigp_results))) {
    write.csv(sigp_results[[i]],file = paste(
      sep = "",getwd(),"/","sigp_results_",names(sigp_results)[i],".csv"
    ));
  }
  
  for (i in seq(1,length(all_results))) {
    write.csv(all_results[[i]],file = paste(
      sep = "",getwd(),"/","all_results_",names(all_results)[i],".csv"
    ));
  }
  
  ########### LSAAR ###########
  leucine_codons = c("TTA","TTG","CTT","CTC","CTA","CTG");
  rs = 1;
  err = 0;
  avail_LSAAR_codons_freq = list();
  for (org in seq(1,length(sigpL_names),by = 1)) {
    temp_avail_LSAAR_codons_freq = matrix("",nrow = length(sigpL[[org]]$SEQ), length(leucine_codons) +
                                            1)
    for (s in seq(1,length(cDNA[[org]]$ID),by = 1)) {
      if (length(sigpL[[org]]$LSAAR_LENGTH[which(cDNA[[org]]$ID[s] == sigpL[[org]]$PROTEIN_ID)]) !=
          0) {
        #if exists: length of L-SAAR in signal peptide
        cDNA_corr = (s2c(cDNA[[org]]$SEQUENCE[s])[min(as.numeric(cDNA[[org]]$START[s][[1]])):max(as.numeric(cDNA[[org]]$STOP[s][[1]]))]); #cDNA - sequence of the whole peptide that contains the SP
        cDNA_corr_proper = cDNA_corr[seq(((sigpL[[org]]$LSAAR_START[which(cDNA[[org]]$ID[s] ==
                                                                            sigpL[[org]]$PROTEIN_ID)]) * 3 + 1),((sigpL[[org]]$LSAAR_LENGTH[which(cDNA[[org]]$ID[s] ==
                                                                                                                                                    sigpL[[org]]$PROTEIN_ID)]) * 3) + ((sigpL[[org]]$LSAAR_START[which(cDNA[[org]]$ID[s] ==
                                                                                                                                                                                                                         sigpL[[org]]$PROTEIN_ID)]) * 3))]; #cDNA - sequence corresponding exactly to the SP
        if (length(which(cDNA_corr_proper == "N") != 0)) {
          cDNA_corr_proper[which(cDNA_corr_proper == "N")] = "A"
          err = err + 1;
        }
        if (length(which(is.na(cDNA_corr_proper)) != 0)) {
          cDNA_corr_proper[which(is.na(cDNA_corr_proper))] = "A"
          err = err + 1;
        }
        if (rs %% 100 == 0) {
          print(paste((round(
            rs / length(sigpL[[org]]$PROTEIN_ID) * 100
          )),"% DONE, organism: ",sigpL_names[org],sep = ""))
        }
        codons_freq = table(as.character(codons(DNAString(
          c2s(cDNA_corr_proper)
        ))));
        temp_avail_LSAAR_codons_freq[rs,] = append(codons_freq[leucine_codons],sigpL[[org]]$PROTEIN_ID[which(cDNA[[org]]$ID[s] ==
                                                                                                               sigpL[[org]]$PROTEIN_ID)]);
        temp_avail_LSAAR_codons_freq[rs,][which(is.na(temp_avail_LSAAR_codons_freq[rs,]))] =
          0;
        rs = rs + 1;
      }
    }
    rs = 0;
    avail_LSAAR_codons_freq[[org]] = temp_avail_LSAAR_codons_freq;
  }
  
  names(avail_LSAAR_codons_freq) = sigpL_names;
  for (i in seq(1,length(sigpL_names),by = 1)) {
    colnames(avail_LSAAR_codons_freq[[i]]) = append(leucine_codons,"PROTEIN_ID");
  }
  
  
  setwd(
    "~/Uczelnia/MGR/praca_magisterska/analiza_kodonow/cDNA/codons_freq_results/LSAARsigp/"
  );
  for (org in seq(1,length(prot_names),by = 1)) {
    write.csv(
      avail_LSAAR_codons_freq[[org]],file = paste(sep = "",getwd(),"/",prot_names[org],"_LSAARsigp_result.csv")
    );
  }
  
  mean_LSAAR_results = matrix("",nrow = length(cDNA_names),ncol = 6);
  sd_LSAAR_results = matrix("",nrow = length(cDNA_names),ncol = 6);
  sum_LSAAR_results = matrix(0,nrow = length(cDNA_names),ncol = 7);
  per_LSAAR_results = matrix(0,nrow = length(cDNA_names),ncol = 6);
  
  
  for (org in seq(1,length(codon_names),by = 1)) {
    for (cod in seq(1,6,by = 1)) {
      mean_LSAAR_results[org,cod] = mean(as.numeric(avail_LSAAR_codons_freq[[org]][,cod])[which(!is.na(as.numeric(avail_LSAAR_codons_freq[[org]][,cod])))]);
      sd_LSAAR_results[org,cod] = sd(as.numeric(avail_LSAAR_codons_freq[[org]][,cod])[which(!is.na(as.numeric(avail_LSAAR_codons_freq[[org]][,cod])))]);
      sum_LSAAR_results[org,cod] = sum(as.numeric(avail_LSAAR_codons_freq[[org]][,cod])[which(!is.na(as.numeric(avail_LSAAR_codons_freq[[org]][,cod])))]);
    }
    sum_LSAAR_results[org,7] = sum(as.numeric(sum_LSAAR_results[org,]));
  }
  
  
  rownames(mean_LSAAR_results) = codon_names;colnames(mean_LSAAR_results) =
    leucine_codons;
  rownames(sd_LSAAR_results) = codon_names;colnames(sd_LSAAR_results) = leucine_codons;
  rownames(sum_LSAAR_results) = codon_names;colnames(sum_LSAAR_results) =
    append(leucine_codons,"suma");
  rownames(per_LSAAR_results) = codon_names;colnames(per_LSAAR_results) =
    leucine_codons;
  
  for (org in seq(1,length(codon_names),by = 1)) {
    for (cod in seq(1,6,by = 1)) {
      per_LSAAR_results[org,cod] = sum_LSAAR_results[org,cod] / sum_LSAAR_results[org,7];
    }
  }
  
  LSAAR_results = list(
    mean = mean_LSAAR_results,percentage = per_LSAAR_results,standard_deviation =
      sd_LSAAR_results,sum = sum_LSAAR_results
  )
  setwd(
    "~/Uczelnia/MGR/praca_magisterska/analiza_kodonow/cDNA/codons_freq_results/statresults"
  );
  for (i in seq(1,length(LSAAR_results))) {
    write.csv(LSAAR_results[[i]],file = paste(
      sep = "",getwd(),"/","LSAAR_results_",names(LSAAR_results)[i],".csv"
    ));
  }
  
  ########### Analysis only for avaialable ortho ######### 
  
  setwd("~/Uczelnia/MGR/praca_magisterska/analiza_kodonow/mappers/");
  myFiles_ortho <-
    list.files(pattern = "^homo_sapiens*"); #divided into separate files because ensemble enables only 6 orthologous organisms at once.
  ortho_mapper = list();
  ortho_mapper_names = c();
  for (k in 1:length(myFiles_ortho)) {
    ortho_mapper[[k]] = read.csv(myFiles_ortho[k],stringsAsFactors = F);
    ortho_mapper_names[k] = strsplit(myFiles_ortho,"[.]")[[k]][1];
  }
  names(ortho_mapper) = ortho_mapper_names;
  
  #ortho consists of rows which cosnsist each organism id (of two: human and ortgologous)
  a = vector("list", 2)
  ortho = list(a,a,a,a,a,a,a,a);
  for (org in seq(1,length(ortho_mapper),by = 1)) {
    print(org)
    num = 1;
    for (i in seq(1,length(ortho_mapper[[org]][[1]]),by = 1)) {
      if (ortho_mapper[[org]][[1]][i] != "" &&
          ortho_mapper[[org]][[2]][i] != "") {
        ortho[[org]][[2]][num] = ortho_mapper[[org]][[2]][i];
        ortho[[org]][[1]][num] = ortho_mapper[[org]][[1]][i];
        num = num + 1;
      }
    }
  }
  
  names(ortho) = ortho_mapper_names;
  
  #orthologoues conists of rows containing only unique human ids
  # time consuming as hell
  a = vector("list", 2)
  orthologues = list(a,a,a,a,a,a,a,a);
  for (org in seq(1,length(ortho))) {
    print(org)
    for (i in seq(2,length(unique(ortho_mapper[[1]][[1]])))) {
      ktory = which(ortho[[org]][[1]] == unique(ortho_mapper[[1]][[1]])[i])[1];
      orthologues[[org]][[1]][i] = ortho[[org]][[1]][ktory];
      orthologues[[org]][[2]][i] = ortho[[org]][[2]][ktory];
      if (i %% 1000 == 0) {
        print(paste((round(
          i / length(unique(ortho_mapper[[1]][[1]])) * 100
        )),"% DONE, organism: ",sigpL_names[org],sep = ""))
      }
    }
  }
  
  #usunięcie wierszy z NA
  for (org in seq(1,length(orthologues))) {
    orthologues[[org]][[1]] = orthologues[[org]][[1]][-which((is.na(orthologues[[org]][[1]])))];
    orthologues[[org]][[2]] = orthologues[[org]][[2]][-which((is.na(orthologues[[org]][[2]])))];
  }
  
  names(orthologues) = ortho_mapper_names;
  
  unique_common_ids = unique(Reduce(
    intersect, list(
      orthologues[[1]][[1]],orthologues[[2]][[1]],orthologues[[3]][[1]],orthologues[[4]][[1]],orthologues[[5]][[1]],orthologues[[6]][[1]],orthologues[[7]][[1]],orthologues[[8]][[1]]
    )
  ))
  2
  #import packages
  library(foreach)
  library(doParallel)
  
  #setup parallel backend to use 8 processors
  cl <- makeCluster(3)
  registerDoParallel(cl)
  
  #start time
  strt <- Sys.time()
  
  a = c();
  common_orthologues = list(a,a,a,a,a,a,a,a,a);
  foreach(org = 1:8) %dopar% {
    for (i in seq(1,length(unique_common_ids),by = 1)) {
      common_orthologues[[org]][i] = orthologues[[org]][[2]][which(orthologues[[org]][[1]] == unique_common_ids[i])];
    }
  }
  
  print(Sys.time() - strt)
  stopCluster(cl)
  
  #Adding artificicial element to the list containing homo_sapiens ensids
  common_orthologues[[9]] = list(a);
  common_orthologues[4:9] = c(list(a),common_orthologues[4:8])
  common_orthologues[[4]] = list(common_orthologues[[1]][[1]],common_orthologues[[1]][[1]])
  ortho_mapper_names[4:9] = c("homo_sapiens_homo_sapiens",ortho_mapper_names[4:8])
  names(common_orthologues) = ortho_mapper_names
  
  most_freq_L_codon_ortho = matrix("",length(common_orthologues),length(common_orthologues[[org]][[2]]));
  most_freq_L_codon_names_ortho = matrix("",length(common_orthologues),length(common_orthologues[[org]][[2]]));
  leucine_codons = c("TTA","TTG","CTT","CTC","CTA","CTG");
  for (org in seq(1,length(common_orthologues),by = 1)) {
    result = c();
    result = matrix("",length(common_orthologues[[org]][[2]]) * 2,7);
    for (i in seq(1,length(common_orthologues[[org]][[2]]),by = 1)) {
      if (i %% 1000 == 1) {
        print(paste(round(i / (
          length(common_orthologues[[org]][[2]])
        ) * 100),"% DONE, organism: ",prot_names[org],sep = ""));
      }
      cdna = cDNA[[org]]$SEQUENCE[cDNA[[org]]$cDNA_ID == mapper[[org]]$cDNA_ID[which(mapper[[org]]$Protein_ID ==
                                                                                       common_orthologues[[org]][[2]][i])]]
      #protein = proteins[[org]]$SEQUENCE[which(proteins[[org]]$ID == mapper[[org]]$Protein_ID[i])]
      startt = c();
      stopp = c();
      startt = sort(as.numeric(cDNA[[org]]$START[cDNA[[org]]$cDNA_ID == mapper[[org]]$cDNA_ID[which(mapper[[org]]$Protein_ID ==
                                                                                                      common_orthologues[[org]][[2]][i])]][[1]]))
      stopp = sort(as.numeric(cDNA[[org]]$STOP[cDNA[[org]]$cDNA_ID == mapper[[org]]$cDNA_ID[which(mapper[[org]]$Protein_ID ==
                                                                                                    common_orthologues[[org]][[2]][i])]][[1]]))
      cdna = s2c(cdna);
      vec = c();
      for (j in seq(1,length(stopp),by = 1)) {
        vec = append(x = vec,values = seq(startt[j],stopp[j],by = 1));
      }
      odp = cdna[vec];
      if (length(which(odp == "N")) != 0) {
        odp[which(odp == "N")] = "A"
      }
      if (!is.na(odp[length(odp)])) {
        codons_freq = table(as.character(codons(DNAString(c2s(
          odp
        )))));
        avail_L_codons_freq = codons_freq[leucine_codons];
        result[(i * 2) - 1,] = append(
          names(avail_L_codons_freq),mapper$protein_cdna_mapper_homo_sapiens$cDNA_ID[which(mapper[[org]]$Protein_ID ==
                                                                                             common_orthologues[[org]][[2]][i])]
        )
        result[i * 2,] = append(as.vector(avail_L_codons_freq),0)
      }
      else {
        result[(i * 2) - 1,] = append(rep("-",6),mapper$protein_cdna_mapper_homo_sapiens$cDNA_ID[which(mapper[[org]]$Protein_ID ==
                                                                                                         common_orthologues[[org]][[2]][i])])
        result[i * 2,] = append(rep("-",6),0)
      }
    }
    write.csv(result,file = paste(sep = "",getwd(),"/",prot_names[org],"_ortho_only_result.csv"));
    print(paste("ORGANISM: ",org,"-> DONE.",sep = ""))
  }
  
  setwd(
    "~/Uczelnia/MGR/praca_magisterska/analiza_kodonow/cDNA/codons_freq_results/ortho_only/"
  );
  myFiles <- list.files(pattern = "*.csv");
  leucine_codons = c("TTA","TTG","CTT","CTC","CTA","CTG");
  codon_ortho = list();
  codon_names_ortho = c();
  # Loop reading .csv files containing L codons frequencies
  for (k in 1:length(myFiles)) {
    codon_names_ortho[k] = strsplit(myFiles,"[.]")[[k]][1];
    codon_ortho[[k]] = read.csv(myFiles[k],stringsAsFactors = F);
    codon_ortho[[k]] = codon_ortho[[k]][,-1];
    names(codon_ortho[[k]]) = append(leucine_codons,"SEQ_ID");
  }
  names(codon_ortho) = codon_names_ortho;
  
  # Loop extracting the frequencies and seq_ids only
  for (org in seq(1,length(codon_names_ortho),by = 1)) {
    for (cod in seq(1,6,by = 1)) {
      codon_ortho[[org]][[cod]] = codon_ortho[[org]][[cod]][seq(2,length(codon_ortho[[org]][[cod]]),by = 2)];
    }
    codon_ortho[[org]][[7]] = codon_ortho[[org]][[7]][seq(1,length(codon_ortho[[org]][[7]]),by = 2)]
  }
  
  # Loop changing NAs to zeros, since that's what NA means in this case. Also, it converts chracters to numeric
  for (org in seq(1,length(codon_names_ortho),by = 1)) {
    for (cod in seq(1,6,by = 1)) {
      NAs = which(is.na(as.numeric(codon_ortho[[org]][[cod]])));
      codon_ortho[[org]][[cod]][NAs] = 0;
      codon_ortho[[org]][[cod]] = as.numeric(codon_ortho[[org]][[cod]]);
    }
  }
  
  mean_results_ortho = matrix("",nrow = length(codon_names_ortho),ncol = 6);
  sd_results_ortho = matrix("",nrow = length(codon_names_ortho),ncol = 6);
  sum_results_ortho = matrix(0,nrow = length(codon_names_ortho),ncol = 7);
  per_results_ortho = matrix(0,nrow = length(codon_names_ortho),ncol = 6);
  
  
  for (org in seq(1,length(codon_names_ortho),by = 1)) {
    for (cod in seq(1,6,by = 1)) {
      mean_results_ortho[org,cod] = mean(codon_ortho[[org]][[cod]]);
      sd_results_ortho[org,cod] = sd(codon_ortho[[org]][[cod]]);
      sum_results_ortho[org,cod] = sum(as.numeric(codon_ortho[[org]][[cod]]));
    }
    sum_results_ortho[org,7] = sum(as.numeric(sum_results_ortho[org,]));
  }
  
  for (org in seq(1,length(codon_names_ortho),by = 1)) {
    for (cod in seq(1,6,by = 1)) {
      per_results_ortho[org,cod] = sum_results_ortho[org,cod] / sum_results_ortho[org,7];
    }
  }
  
  rownames(mean_results_ortho) = codon_names_ortho;colnames(mean_results_ortho) = leucine_codons;
  rownames(sd_results_ortho) = codon_names_ortho;colnames(sd_results_ortho) = leucine_codons;
  rownames(sum_results_ortho) = codon_names_ortho;colnames(sum_results_ortho) = append(leucine_codons,"TOTAL");
  rownames(per_results_ortho) = codon_names_ortho;colnames(per_results_ortho) = leucine_codons;
  
  all_results_ortho = list(
    mean = mean_results_ortho,percentage = per_results_ortho,standard_deviation = sd_results_ortho,sum =
      sum_results_ortho
  )
  
  setwd(
    "~/Uczelnia/MGR/praca_magisterska/analiza_kodonow/cDNA/codons_freq_results/statresults"
  );
  for (i in seq(1,length(all_results_ortho))) {
    write.csv(
      all_results_ortho[[i]],file = paste(
        sep = "",getwd(),"/","ortho_only_all_results_",names(LSAAR_results)[i],".csv"
      )
    );
  }

  
  ####### CORRECTED COMMON ORTHOLOGUES ###########
  
  
  indices = list(0,0,0,0,0,0,0,0,0)
  #Improvement - no mendacious results?
  for (org in seq(1,length(codon_names),by = 1)) {
    print(org)
    for (i in seq(1,length(unique(common_orthologues[[org]][[2]])),by = 1)) {
      indices[[org]] = append(indices[[org]],which(common_orthologues[[org]][[2]] ==
                                                     unique(common_orthologues[[org]][[2]])[i])[1])
    }
  }
  
  unique_indices = unique(Reduce(
    intersect, list(
      indices[[1]],indices[[2]],indices[[3]],indices[[4]],indices[[5]],indices[[6]],indices[[7]],indices[[8]],indices[[9]]
    )
  ))
  
  a = c();
  corrected_common_orthologues = list(a,a,a,a,a,a,a,a,a);
  for (org in seq(1,length(codon_names),by = 1)) {
    for (i in seq(1,2,by = 1)) {
      corrected_common_orthologues[[org]][[i]] = common_orthologues[[org]][[i]][unique_indices];
    }
  }
  names(corrected_common_orthologues) = ortho_mapper_names
  
  
  most_freq_L_codon_ortho_corrected = matrix("",length(corrected_common_orthologues),length(corrected_common_orthologues[[org]][[2]]));
  most_freq_L_codon_names_ortho_corrected = matrix("",length(corrected_common_orthologues),length(corrected_common_orthologues[[org]][[2]]));
  leucine_codons = c("TTA","TTG","CTT","CTC","CTA","CTG");
  for (org in seq(1,length(corrected_common_orthologues),by = 1)) {
    result = c();
    result = matrix("",length(corrected_common_orthologues[[org]][[2]]) * 2,7);
    for (i in seq(1,length(corrected_common_orthologues[[org]][[2]]),by = 1)) {
      if (i %% 1000 == 1) {
        print(paste(round(i / (
          length(corrected_common_orthologues[[org]][[2]])
        ) * 100),"% DONE, organism: ",prot_names[org],sep = ""));
      }
      cdna = cDNA[[org]]$SEQUENCE[cDNA[[org]]$cDNA_ID == mapper[[org]]$cDNA_ID[which(mapper[[org]]$Protein_ID ==
                                                                                       corrected_common_orthologues[[org]][[2]][i])]]
      #protein = proteins[[org]]$SEQUENCE[which(proteins[[org]]$ID == mapper[[org]]$Protein_ID[i])]
      startt = c();
      stopp = c();
      startt = sort(as.numeric(cDNA[[org]]$START[cDNA[[org]]$cDNA_ID == mapper[[org]]$cDNA_ID[which(mapper[[org]]$Protein_ID ==
                                                                                                      corrected_common_orthologues[[org]][[2]][i])]][[1]]))
      stopp = sort(as.numeric(cDNA[[org]]$STOP[cDNA[[org]]$cDNA_ID == mapper[[org]]$cDNA_ID[which(mapper[[org]]$Protein_ID ==
                                                                                                    corrected_common_orthologues[[org]][[2]][i])]][[1]]))
      cdna = s2c(cdna);
      vec = c();
      for (j in seq(1,length(stopp),by = 1)) {
        vec = append(x = vec,values = seq(startt[j],stopp[j],by = 1));
      }
      odp = cdna[vec];
      if (length(which(odp == "N")) != 0) {
        odp[which(odp == "N")] = "A"
      }
      if (!is.na(odp[length(odp)])) {
        codons_freq = table(as.character(codons(DNAString(c2s(
          odp
        )))));
        avail_L_codons_freq = codons_freq[leucine_codons];
        result[(i * 2) - 1,] = append(
          names(avail_L_codons_freq),mapper$protein_cdna_mapper_homo_sapiens$cDNA_ID[which(mapper[[org]]$Protein_ID ==
                                                                                             corrected_common_orthologues[[org]][[2]][i])]
        )
        result[i * 2,] = append(as.vector(avail_L_codons_freq),0)
      }
      else {
        result[(i * 2) - 1,] = append(rep("-",6),mapper$protein_cdna_mapper_homo_sapiens$cDNA_ID[which(mapper[[org]]$Protein_ID ==
                                                                                                         corrected_common_orthologues[[org]][[2]][i])])
        result[i * 2,] = append(rep("-",6),0)
      }
    }
    write.csv(result,file = paste(sep = "",getwd(),"/",prot_names[org],"corrected_ortho_only_result.csv"));
    print(paste("ORGANISM: ",org,"-> DONE.",sep = ""))
  }
  
  
  
  
  setwd(
    "~/Uczelnia/MGR/praca_magisterska/analiza_kodonow/cDNA/codons_freq_results/corrected_ortho_only/"
  );
  myFiles <- list.files(pattern = "*.csv");
  leucine_codons = c("TTA","TTG","CTT","CTC","CTA","CTG");
  codon_ortho_corrected = list();
  codon_names_ortho_corrected = c();
  # Loop reading .csv files containing L codons frequencies
  for (k in 1:length(myFiles)) {
    codon_names_ortho_corrected[k] = strsplit(myFiles,"[.]")[[k]][1];
    codon_ortho_corrected[[k]] = read.csv(myFiles[k],stringsAsFactors = F);
    codon_ortho_corrected[[k]] = codon_ortho_corrected[[k]][,-1];
    names(codon_ortho_corrected[[k]]) = append(leucine_codons,"SEQ_ID");
  }
  names(codon_ortho_corrected) = codon_names_ortho_corrected;
  
  # Loop extracting the frequencies and seq_ids only
  for (org in seq(1,length(codon_names_ortho_corrected),by = 1)) {
    for (cod in seq(1,6,by = 1)) {
      codon_ortho_corrected[[org]][[cod]] = codon_ortho_corrected[[org]][[cod]][seq(2,length(codon_ortho_corrected[[org]][[cod]]),by = 2)];
    }
    codon_ortho_corrected[[org]][[7]] = codon_ortho_corrected[[org]][[7]][seq(1,length(codon_ortho_corrected[[org]][[7]]),by = 2)]
  }
  
  # Loop changing NAs to zeros, since that's what NA means in this case. Also, it converts chracters to numeric
  for (org in seq(1,length(codon_names_ortho_corrected),by = 1)) {
    for (cod in seq(1,6,by = 1)) {
      NAs = which(is.na(as.numeric(codon_ortho_corrected[[org]][[cod]])));
      codon_ortho_corrected[[org]][[cod]][NAs] = 0;
      codon_ortho_corrected[[org]][[cod]] = as.numeric(codon_ortho_corrected[[org]][[cod]]);
    }
  }
  
  mean_results_ortho_corrected = matrix("",nrow = length(codon_names_ortho_corrected),ncol = 6);
  sd_results_ortho_corrected = matrix("",nrow = length(codon_names_ortho_corrected),ncol = 6);
  sum_results_ortho_corrected = matrix(0,nrow = length(codon_names_ortho_corrected),ncol = 7);
  per_results_ortho_corrected = matrix(0,nrow = length(codon_names_ortho_corrected),ncol = 6);
  
  
  for (org in seq(1,length(codon_names_ortho_corrected),by = 1)) {
    for (cod in seq(1,6,by = 1)) {
      mean_results_ortho_corrected[org,cod] = mean(codon_ortho_corrected[[org]][[cod]]);
      sd_results_ortho_corrected[org,cod] = sd(codon_ortho_corrected[[org]][[cod]]);
      sum_results_ortho_corrected[org,cod] = sum(as.numeric(codon_ortho_corrected[[org]][[cod]]));
    }
    sum_results_ortho_corrected[org,7] = sum(as.numeric(sum_results_ortho_corrected[org,]));
  }
  
  for (org in seq(1,length(codon_names_ortho_corrected),by = 1)) {
    for (cod in seq(1,6,by = 1)) {
      per_results_ortho_corrected[org,cod] = sum_results_ortho_corrected[org,cod] / sum_results_ortho_corrected[org,7];
    }
  }
  
  rownames(mean_results_ortho_corrected) = codon_names_ortho_corrected;colnames(mean_results_ortho_corrected) = leucine_codons;
  rownames(sd_results_ortho_corrected) = codon_names_ortho_corrected;colnames(sd_results_ortho_corrected) = leucine_codons;
  rownames(sum_results_ortho_corrected) = codon_names_ortho_corrected;colnames(sum_results_ortho_corrected) = append(leucine_codons,"TOTAL");
  rownames(per_results_ortho_corrected) = codon_names_ortho_corrected;colnames(per_results_ortho_corrected) = leucine_codons;
  
  all_results_ortho = list(
    mean = mean_results_ortho_corrected,percentage = per_results_ortho_corrected,standard_deviation = sd_results_ortho_corrected,sum =
      sum_results_ortho_corrected
  )
  
  setwd(
    "~/Uczelnia/MGR/praca_magisterska/analiza_kodonow/cDNA/codons_freq_results/statresults"
  );
  for (i in seq(1,length(all_results_ortho))) {
    write.csv(
      all_results_ortho[[i]],file = paste(
        sep = "",getwd(),"/","ortho_corrected_only_all_results_",names(LSAAR_results)[i],".csv"
      )
    );
  }
  
  
  bos_taurus_codon_usage = matrix(NaN,nrow = 3,ncol = 7)
  bos_taurus_codon_usage[1,] = sum_LSAAR_results[1,] / sum_LSAAR_results[1,7];
  bos_taurus_codon_usage[2,] = sum_results[1,] / sum_results[1,7];
  bos_taurus_codon_usage[3,] = sum_sigp_results[1,] / sum_sigp_results[1,7];
  colnames(bos_taurus_codon_usage) = colnames(sum_LSAAR_results)
  rownames(bos_taurus_codon_usage) = c(
    "Leucine codons in L-SAARs", "Leucine codons in general", "Leucine codons in signal peptides"
  )
  
  gallus_gallus_codon_usage = matrix(NaN,nrow = 3,ncol = 7)
  gallus_gallus_codon_usage[1,] = sum_LSAAR_results[2,] / sum_LSAAR_results[2,7];
  gallus_gallus_codon_usage[2,] = sum_results[2,] / sum_results[2,7];
  gallus_gallus_codon_usage[3,] = sum_sigp_results[2,] / sum_sigp_results[2,7];
  colnames(gallus_gallus_codon_usage) = colnames(sum_LSAAR_results)
  rownames(gallus_gallus_codon_usage) = c(
    "Leucine codons in L-SAARs", "Leucine codons in general", "Leucine codons in signal peptides"
  )
  
  gorilla_gorilla_codon_usage = matrix(NaN,nrow = 3,ncol = 7)
  gorilla_gorilla_codon_usage[1,] = sum_LSAAR_results[3,] / sum_LSAAR_results[3,7];
  gorilla_gorilla_codon_usage[2,] = sum_results[3,] / sum_results[3,7];
  gorilla_gorilla_codon_usage[3,] = sum_sigp_results[3,] / sum_sigp_results[3,7];
  colnames(gorilla_gorilla_codon_usage) = colnames(sum_LSAAR_results)
  rownames(gorilla_gorilla_codon_usage) = c(
    "Leucine codons in L-SAARs", "Leucine codons in general", "Leucine codons in signal peptides"
  )
  
  homo_sapiens_codon_usage = matrix(NaN,nrow = 3,ncol = 7)
  homo_sapiens_codon_usage[1,] = sum_LSAAR_results[4,] / sum_LSAAR_results[4,7];
  homo_sapiens_codon_usage[2,] = sum_results[4,] / sum_results[4,7];
  homo_sapiens_codon_usage[3,] = sum_sigp_results[4,] / sum_sigp_results[4,7];
  colnames(homo_sapiens_codon_usage) = colnames(sum_LSAAR_results)
  rownames(homo_sapiens_codon_usage) = c(
    "Leucine codons in L-SAARs", "Leucine codons in general", "Leucine codons in signal peptides"
  )
  
  macaca_mulatta_codon_usage = matrix(NaN,nrow = 3,ncol = 7)
  macaca_mulatta_codon_usage[1,] = sum_LSAAR_results[5,] / sum_LSAAR_results[5,7];
  macaca_mulatta_codon_usage[2,] = sum_results[5,] / sum_results[5,7];
  macaca_mulatta_codon_usage[3,] = sum_sigp_results[5,] / sum_sigp_results[5,7];
  colnames(macaca_mulatta_codon_usage) = colnames(sum_LSAAR_results)
  rownames(macaca_mulatta_codon_usage) = c(
    "Leucine codons in L-SAARs", "Leucine codons in general", "Leucine codons in signal peptides"
  )
  
  mus_musculus_codon_usage = matrix(NaN,nrow = 3,ncol = 7)
  mus_musculus_codon_usage[1,] = sum_LSAAR_results[6,] / sum_LSAAR_results[6,7];
  mus_musculus_codon_usage[2,] = sum_results[6,] / sum_results[6,7];
  mus_musculus_codon_usage[3,] = sum_sigp_results[6,] / sum_sigp_results[6,7];
  colnames(mus_musculus_codon_usage) = colnames(sum_LSAAR_results)
  rownames(mus_musculus_codon_usage) = c(
    "Leucine codons in L-SAARs", "Leucine codons in general", "Leucine codons in signal peptides"
  )
  
  pan_troglodytes_codon_usage = matrix(NaN,nrow = 3,ncol = 7)
  pan_troglodytes_codon_usage[1,] = sum_LSAAR_results[7,] / sum_LSAAR_results[7,7];
  pan_troglodytes_codon_usage[2,] = sum_results[7,] / sum_results[7,7];
  pan_troglodytes_codon_usage[3,] = sum_sigp_results[7,] / sum_sigp_results[7,7];
  colnames(pan_troglodytes_codon_usage) = colnames(sum_LSAAR_results)
  rownames(pan_troglodytes_codon_usage) = c(
    "Leucine codons in L-SAARs", "Leucine codons in general", "Leucine codons in signal peptides"
  )
  
  rattus_norvegicus_codon_usage = matrix(NaN,nrow = 3,ncol = 7)
  rattus_norvegicus_codon_usage[1,] = sum_LSAAR_results[8,] / sum_LSAAR_results[8,7];
  rattus_norvegicus_codon_usage[2,] = sum_results[8,] / sum_results[8,7];
  rattus_norvegicus_codon_usage[3,] = sum_sigp_results[8,] / sum_sigp_results[8,7];
  colnames(rattus_norvegicus_codon_usage) = colnames(sum_LSAAR_results)
  rownames(rattus_norvegicus_codon_usage) = c(
    "Leucine codons in L-SAARs", "Leucine codons in general", "Leucine codons in signal peptides"
  )
  
  xenopus_tropicalis_codon_usage = matrix(NaN,nrow = 3,ncol = 7)
  xenopus_tropicalis_codon_usage[1,] = sum_LSAAR_results[9,] / sum_LSAAR_results[9,7];
  xenopus_tropicalis_codon_usage[2,] = sum_results[9,] / sum_results[9,7];
  xenopus_tropicalis_codon_usage[3,] = sum_sigp_results[9,] / sum_sigp_results[9,7];
  colnames(xenopus_tropicalis_codon_usage) = colnames(sum_LSAAR_results)
  rownames(xenopus_tropicalis_codon_usage) = c(
    "Leucine codons in L-SAARs", "Leucine codons in general", "Leucine codons in signal peptides"
  )
  
  codon_usage=list(bos_taurus_codon_usage,gallus_gallus_codon_usage,gorilla_gorilla_codon_usage,homo_sapiens_codon_usage,macaca_mulatta_codon_usage,mus_musculus_codon_usage,pan_troglodytes_codon_usage,rattus_norvegicus_codon_usage,xenopus_tropicalis_codon_usage);
  names(codon_usage)=prot_names
  
  for (i in seq(1,length(codon_usage),by = 1)){
    fname = paste("/home/mstolarczyk/Uczelnia/MGR/praca_magisterska/analiza_kodonow/cDNA/codons_freq_results/codon_usage/",prot_names[i],"_codon_usage.csv",sep = "");
    write.csv(x = codon_usage[[i]],file = fname,row.names = T,col.names = T)
  }
  
  
  
  
  
  
# Codon changes within L --------------------------------------------------
  cDNA_test = list(a,a,a,a,a,a,a,a,a);
for (org in seq(1,length(cDNA),1)){
  print(paste("Organism: ",sep = "",names(cDNA)[org]))
  for (i in seq(1,length(cDNA[[org]]$STOP))){
    vec=c();
    starting=min(as.numeric(unlist(cDNA[[org]]$START[i])))#coding sequence start
    stopping=max(as.numeric(unlist(cDNA[[org]]$STOP[i])))#coding sequence end
    vec = seq(starting,stopping,by = 1)#numbers of coding nucleotides
    cseq=c2s(s2c(cDNA[[org]]$SEQUENCE[i])[vec]); #proper coding sequence
    cDNA[[org]]$SEQUENCE[i] = cseq
  }
}
  # Parallel computing ------------------------------------------------------
  #setup parallel backend to use 3processors
  cl <- makeCluster(3)
  registerDoParallel(cl)
  
  #start timecDN
  strt <- Sys.time()
  
  a = c();
  cDNA_test = list(a,a,a,a,a,a,a,a,a);
  foreach(org = 1:9,.packages='seqinr') %dopar%
    {
      for (i in seq(1,length(cDNA[[org]]$STOP))){
        vec=c();
        starting=min(as.numeric(unlist(cDNA[[org]]$START[i])))#coding sequence start
        stopping=max(as.numeric(unlist(cDNA[[org]]$STOP[i])))#coding sequence end
        vec = seq(starting,stopping,by = 1)#numbers of coding nucleotides
        cseq=c2s(s2c(cDNA[[org]]$SEQUENCE[i])[vec]); #proper coding sequence
        cDNA[[org]]$SEQUENCE[i] = cseq
      }
    }
stopCluster(cl)


# prepare cDNA sequences (0-69*3) to RevTrans -----------------------------
prep2revtrans = list(a,a,a,a,a,a,a,a,a);
for (org in seq(1,length(cDNA),1)){
  prep2revtrans[[org]] = as.matrix(cDNA[[org]])[,c(2,5)]
  for (i in seq(1,length(prep2revtrans[[org]][,2]),1)){
    if (length(s2c(prep2revtrans[[org]][i,2]$SEQUENCE)) >= 3*69){
      prep2revtrans[[org]][i,2] = c2s(s2c(prep2revtrans[[org]][i,2]$SEQUENCE)[1:(3*69)])
    }
    else{
      prep2revtrans[[org]][i,2] = prep2revtrans[[org]][i,2]$SEQUENCE
    }
}
names(prep2revtrans) = prot_names
for (i in seq(1, length(prep2revtrans),1)){
  write.csv(row.names = F,x = prep2revtrans[[i]],file = paste(getwd(),"/","prep2revtrans","_",prot_names[i],".csv",sep = ""))
}
# Save environment --------------------------------------------------------
save(... = ...,file=paste(getwd(),"/",date(),sep = ""))
  



