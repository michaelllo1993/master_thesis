

# Load libraries ----------------------------------------------------------
library(Biostrings)
# Set wd, get arguments ---------------------------------------------------

args = commandArgs(trailingOnly = TRUE)
wd = getwd()
OoI = args[1]
repeat_unit = args[2]
setwd(paste(wd,"/codon_frequency/", OoI, "_orthologues", sep = ""))

# Get the files names -----------------------------------------------------

#myFiles_overall <- list.files(pattern = "*overall.csv")
myFiles_overall_noSP <- list.files(pattern = "*excluding_SP.csv")
myFiles_withinSP <- list.files(pattern = "*frequency_SP.csv")
myFiles_withinSAAR <- list.files(pattern = "*SAAR.csv")

# Read the data into matrices and preprocess them -------------------------
# Overall codon usage
# overall = list()
# overall_names = c()
# for (k in 1:length(myFiles_overall)) {
#   overall_names[k] = paste(strsplit(myFiles_overall[k], "_")[[1]][1],
#                            strsplit(myFiles_overall[k], "_")[[1]][2],
#                            sep = "_")
#   overall[[k]] = read.csv(myFiles_overall[k],
#                           stringsAsFactors = F,
#                           header = T)
#   overall[[k]] = unique(overall[[k]])
#   overall[[k]] = overall[[k]][which(overall[[k]]$cDNA_ID != ""), ]
#   overall[[k]][is.na(overall[[k]])] = 0
# }
# names(overall) = overall_names
# print("DONE")
#Overall\SP
overall_noSP = list()
overall_noSP_names = c()
for (k in 1:length(myFiles_overall_noSP)) {
  overall_noSP_names[k] = paste(
    strsplit(myFiles_overall_noSP[k], "_")[[1]][1],
    strsplit(myFiles_overall_noSP[k], "_")[[1]][2],
    sep = "_"
  )
  overall_noSP[[k]] = read.csv(myFiles_overall_noSP[k],
                               stringsAsFactors = F,
                               header = T)
  overall_noSP[[k]] = unique(overall_noSP[[k]])
  overall_noSP[[k]] = overall_noSP[[k]][which(overall_noSP[[k]]$cDNA_ID !=
                                                ""), ]
  overall_noSP[[k]][is.na(overall_noSP[[k]])] = 0
}
names(overall_noSP) = overall_noSP_names
print("DONE")
# Within SP
withinSP = list()
withinSP_names = c()
for (k in 1:length(myFiles_withinSP)) {
  withinSP_names[k] = paste(strsplit(myFiles_withinSP[k], "_")[[1]][1],
                            strsplit(myFiles_withinSP[k], "_")[[1]][2],
                            sep = "_")
  withinSP[[k]] = read.csv(myFiles_withinSP[k],
                           stringsAsFactors = F,
                           header = T)
  withinSP[[k]] = unique(withinSP[[k]])
  withinSP[[k]] = withinSP[[k]][which(withinSP[[k]]$cDNA_ID != ""), ]
  withinSP[[k]][is.na(withinSP[[k]])] = 0
}
names(withinSP) = withinSP_names
print("DONE")
# Within SAAR
withinSAAR = list()
withinSAAR_names = c()
for (k in 1:length(myFiles_withinSAAR)) {
  withinSAAR_names[k] = paste(strsplit(myFiles_withinSAAR[k], "_")[[1]][1],
                              strsplit(myFiles_withinSAAR[k], "_")[[1]][2],
                              sep = "_")
  withinSAAR[[k]] = read.csv(myFiles_withinSAAR[k],
                             stringsAsFactors = F,
                             header = T)
  withinSAAR[[k]] = unique(withinSAAR[[k]])
  withinSAAR[[k]] = withinSAAR[[k]][which(withinSAAR[[k]]$cDNA_ID != ""), ]
  withinSAAR[[k]][is.na(withinSAAR[[k]])] = 0
}
names(withinSAAR) = withinSAAR_names
print("DONE")

# Analyze the data --------------------------------------------------------

unit_codons = names(which(GENETIC_CODE == repeat_unit))

# Protein excluding SP
codons_protein_fraction = matrix(NaN,
                                 nrow = length(overall_noSP) ,
                                 ncol = length(unit_codons))
for (org in seq_len(length(overall_noSP))) {
  mtx = overall_noSP[[org]][, -1]
  for (cod in seq_len(length(unit_codons))) {
    codons_protein_fraction[org, cod] = sum(mtx[, cod]) / sum(mtx)
  }
}
colnames(codons_protein_fraction) = unit_codons
rownames(codons_protein_fraction) = names(overall_noSP)


codons_SP_fraction = matrix(NaN, nrow = length(withinSP) , ncol = length(unit_codons))
for (org in seq_len(length(withinSP))) {
  mtx = withinSP[[org]][, -1]
  for (cod in seq_len(length(unit_codons))) {
    codons_SP_fraction[org, cod] = sum(mtx[, cod]) / sum(mtx)
  }
}
colnames(codons_SP_fraction) = unit_codons
rownames(codons_SP_fraction) = names(withinSP)


codons_SAAR_fraction = matrix(NaN,
                              nrow = length(withinSAAR) ,
                              ncol = length(unit_codons))
for (org in seq_len(length(withinSAAR))) {
  mtx = withinSAAR[[org]][, -1]
  for (cod in seq_len(length(unit_codons))) {
    codons_SAAR_fraction[org, cod] = sum(mtx[, cod]) / sum(mtx)
  }
}
colnames(codons_SAAR_fraction) = unit_codons
rownames(codons_SAAR_fraction) = names(withinSAAR)


# Visualize the results ---------------------------------------------------
print(paste(
  wd,
  "/codon_frequency/",
  OoI,
  "_orthologues",
  "/plots_",
  OoI,
  "_orthologues/",
  rownames(codons_SAAR_fraction)[1],
  "_codon_usage_region.svg",
  sep = ""
))
for (i in seq(1, dim(codons_SAAR_fraction)[1])) {
  svg(
    filename = paste(
      wd,
      "/codon_frequency/",
      OoI,
      "_orthologues",
      "/codon_frequency_plots_",
      OoI,
      "_orthologues/",
      rownames(codons_SAAR_fraction)[i],
      "_codon_usage_region.svg",
      sep = ""
    ),
    width = 40,
    height = 30,
    pointsize = 50
  )
  barplot(
    rbind(
      codons_protein_fraction[i, ] * 100,
      codons_SP_fraction[i, ] * 100,
      codons_SAAR_fraction[i, ] * 100
    ),
    beside = T,
    col = c("dodgerblue3", "firebrick2", "darkgreen"),
    xlab = "Codon",
    ylab = "Codon usage [%]",
    main = paste("L codon usage for", rownames(codons_SAAR_fraction)[i]),
    ylim = c(0, 100)
  )
  legend(
    "topleft",
    c("The rest of the protein", "Signal peptide", "L-SAAR"),
    col = c("dodgerblue3", "firebrick2", "darkgreen"),
    pch = 15,
    bty = "n"
  )
  dev.off()
}
