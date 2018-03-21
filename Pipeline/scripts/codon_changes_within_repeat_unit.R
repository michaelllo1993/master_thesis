

# Loading libraries -------------------------------------------------------

require(seqinr)
require(Biostrings)

# Getting the command line arguments ------------------------------------------

options(warn = -1)
wd = getwd()
args = commandArgs(trailingOnly = TRUE)

#read the arguments
revtrans_file = args[1]
codes_dict = t(read.csv(args[2], stringsAsFactors = F,header = F))
colnames(codes_dict) = codes_dict[1,]
codes_dict = codes_dict[-1,]
organism_of_interest_name = args[3]
repeat_unit = args[4]

#get all codons vector
all_codon_names = append("---", sort(names(GENETIC_CODE)))
#get the codons encoding the repeat forming AA
unit_codons = names(which(GENETIC_CODE == repeat_unit))
#get the ENSEMBL code of the organism of interest
organism_of_interest_code = as.character(codes_dict[organism_of_interest_name])
#get the organisms in the analysis
organisms = names(codes_dict)
#exclude the organism of interest
tmp_organisms = organisms[-which(organisms == organism_of_interest_name)]
#get the ENSEMBL code of the organism of all the organisms
ensembl_exclusive_codes = as.character(codes_dict)
#exclude the organism of interest
tmp_codes = ensembl_exclusive_codes[-which(ensembl_exclusive_codes == organism_of_interest_code)]

# Reading the data in --------------------------
#read the file
revtrans_data = as.matrix(read.csv(file = revtrans_file, header = F))
#create an empty list that will store the data in an easily accessbile way
organism_of_interest = lapply(1:length(tmp_codes), function(x)
  matrix(
    NaN,
    nrow = dim(revtrans_data)[1],
    ncol = dim(revtrans_data)[2]
  ))

for (org in seq(1, length(tmp_codes), by = 1)) {
  reg_ex = paste(tmp_codes[org], "[0-9]+", sep = "")
  
  #filter just the rows with organism in loop
  organism_of_interest[[org]] = revtrans_data[grep(pattern = reg_ex , revtrans_data[, 3], perl = T), ]
}
#name the list elements appropriately
names(organism_of_interest) = tmp_organisms

# Loop constructing thee lists with codons on leucine codons positions in organism of interest cDNA sequences --------
#create an empty list that will store the data in an easily accessbile way
OoI = lapply(1:length(organism_of_interest), function(x)
  lapply(1:length(unit_codons), function(x)
    NaN))
#create an empty list that will store the data in an easily accessbile way
OoI_SAAR = lapply(1:length(organism_of_interest), function(x)
  lapply(1:length(unit_codons), function(x)
    NaN))

other_codons = c()
other_codons_SAAR = c()

for (org in seq(1, length(tmp_codes), by = 1)) {
  for (l in seq(1, length(unit_codons), by = 1)) {
    for (i in seq(1, dim(organism_of_interest[[org]])[1], by = 1)) {
      #replace faulty "N" characters with "C"
      check_for_ns = s2c(organism_of_interest[[org]][i, 2])
      check_for_ns[which(check_for_ns == "N")] = "C"
      organism_of_interest[[org]][i, 2] = c2s(check_for_ns)
      #split the organism of interest cDNA sequence into codons (3 letter characters)
      organism_of_interest_codons = strsplit(organism_of_interest[[org]][i, 2], "(?<=.{3})", perl = TRUE)[[1]]
      #get the indices of codon forming repeats that is being analyzed in the iteration
      organism_of_interest_repeat_unit_indices = which(organism_of_interest_codons == unit_codons[l])
      #split the other organism cDNA sequence into codons (3 letter characters)
      other_codons_tmp = strsplit(organism_of_interest[[org]][i, 4], "(?<=.{3})", perl = TRUE)[[1]][organism_of_interest_repeat_unit_indices]
      #append to the resulting vector
      other_codons = append(other_codons, other_codons_tmp)
      #convert the gaps to "A" in order to use the translate function
      test = s2c(organism_of_interest[[org]][i, 2])
      test[s2c(organism_of_interest[[org]][i, 2]) == "-"] = "A"
      #translate the cDNA to AAs to search for the SAARs
      aas = s2c(as.vector(translate(DNAStringSet(c2s(
        test
      )))))
      index = grepRaw(c2s(rep(repeat_unit, 5)), (c2s(aas)))
      #execute the following block of code if there is any SAAR
      if (length(index) >= 1) {
        #search for consecutive AAs
        a = rle(aas)
        #get the length of the longest one
        run_length = a$lengths[which(a$values == repeat_unit)][which.max(a$lengths[which(a$values ==
                                                                                           repeat_unit)])]
        #get the indices of the AAs forming SAAR
        organism_of_interest_SAAR_indices = seq(
          from = index,
          to = (index + run_length - 1),
          by = 1
        )
        #get the codons encoding the SAAR
        other_codons_SAAR_tmp1 = strsplit(organism_of_interest[[org]][i, 2], "(?<=.{3})", perl = TRUE)[[1]][organism_of_interest_SAAR_indices]
        #get only the one that is analyzed in the current iteration
        which_lcodon_analyzed = which(other_codons_SAAR_tmp1 == unit_codons[l])
        proper_indices = (index + which_lcodon_analyzed)  - 1
        other_codons_SAAR_tmp = strsplit(organism_of_interest[[org]][i, 4], "(?<=.{3})", perl = TRUE)[[1]][proper_indices]
        #append to the resulting vector
        other_codons_SAAR = append(other_codons_SAAR, other_codons_SAAR_tmp)
      }
    }
    #copy the counted codons into the list
    OoI[[org]][[l]] = sort(table(other_codons), decreasing = T)
    other_codons = c()
    #copy the counted codons into the list
    OoI_SAAR[[org]][[l]] = sort(table(other_codons_SAAR), decreasing = T)
    other_codons_SAAR = c()
    
  }
}
#name the lists approprotely
names(OoI) = tmp_organisms
for (org in seq(1, length(OoI), by = 1)) {
  names(OoI[[org]]) = unit_codons
}

names(OoI_SAAR) = tmp_organisms
for (org in seq(1, length(OoI_SAAR), by = 1)) {
  names(OoI_SAAR[[org]]) = unit_codons
}

# results saving -------------------------------------------------------------

nms = names(OoI)
for (org in seq(1, length(OoI), by = 1)) {
  #create an empty matrix for the results
  output = matrix(data = NA,
                  nrow = length(all_codon_names),
                  ncol = 1)
  #set the rownames to all codons possible
  rownames(output) = all_codon_names
  for (l in seq(1, length(unit_codons), by = 1)) {
    #create a matrix for temporary results
    tmp = as.matrix(OoI[[org]][[l]])
    #search for faulty "N" nucleotides and delete the codons that consist of these
    if (any(grep("N", rownames(tmp), perl = T))) {
      rows_to_delete = grep("N", rownames(tmp), perl = T)
      tmp = tmp[-rows_to_delete, ]
      tmp = as.matrix(tmp)
    }
    #search for faulty lowercase nucleotides and delete the codons that consist of these
    if (any(grep('[a-z]', rownames(tmp), perl = T))) {
      rows_to_delete = grep('[a-z]', rownames(tmp), perl = T)
      tmp = tmp[-rows_to_delete, ]
      tmp = as.matrix(tmp)
    }
    #merge the temporary and overall results by rownames an keep all results
    output = merge(output, tmp, by = "row.names", all = T)
    output = as.matrix(output[, -1])
    rownames(output) = all_codon_names
  }
  output = output[, -1]
  #convert NAs to 0 - since that's what they mean
  output[which(is.na(output))] = 0
  #write the reslts into CSV file
  colnames(output) = unit_codons
  write.csv(
    x = output,
    file = paste(
      wd,
      "/",
      organism_of_interest_name,
      "_changes_within_repeatUnit/codon_changes_within_repeat_unit_",
      nms[org],
      ".csv",
      sep = ""
    ),
    row.names = T
  )
}

# results saving L-SAARs --------------------------------------------------

nms = names(OoI_SAAR)
for (org in seq(1, length(OoI_SAAR), by = 1)) {
  output = matrix(data = NA,
                  nrow = length(all_codon_names),
                  ncol = 1)
  rownames(output) = all_codon_names
  for (l in seq(1, length(unit_codons), by = 1)) {
    tmp = as.matrix(OoI_SAAR[[org]][[l]])
    if (any(grep("N", rownames(tmp), perl = T))) {
      rows_to_delete = grep("N", rownames(tmp), perl = T)
      tmp = tmp[-rows_to_delete, ]
    }
    output = merge(output, tmp, by = "row.names", all = T)
    output = as.matrix(output[, -1])
    rownames(output) = all_codon_names
  }
  output = output[, -1]
  output[which(is.na(output))] = 0
  
  colnames(output) = unit_codons
  write.csv(
    x = output,
    file = paste(
      wd,
      "/",
      organism_of_interest_name,
      "_changes_within_repeatUnit/codon_changes_within_SAAR_",
      nms[org],
      ".csv",
      sep = ""
    ),
    row.names = T
  )
}


# (repeat unit)\SAAR - calculate the differences and write files ----------------------------

setwd(paste(
  organism_of_interest_name,
  "_changes_within_repeatUnit",
  sep = ""
))
for (org in seq_len(length(nms))) {
  myFiles1 <- list.files(pattern = paste("_unit_",nms[org], ".csv", sep = ""))
  myFiles2 <- list.files(pattern = paste("_SAAR_",nms[org], ".csv", sep = ""))
  X = read.csv(myFiles1)[, 1]
  repeat_unit_data = read.csv(myFiles1)[, -1]
  SAAR_data = read.csv(myFiles2)[, -1]
  print(dim(repeat_unit_data))
  print(dim(SAAR_data))
  if (all(dim(repeat_unit_data) == dim(SAAR_data))) {
    repeat_unit_no_SAAR_data = repeat_unit_data - SAAR_data
    repeat_unit_no_SAAR_data = cbind(X, repeat_unit_no_SAAR_data)
    write.csv(
      x = repeat_unit_no_SAAR_data,
      file = paste(
        wd,
        "/",
        organism_of_interest_name,
        "_changes_within_repeatUnit/codon_changes_within_repeatUnitNoSAAR_",
        nms[org],
        ".csv",
        sep = ""
      ),
      row.names = T
    )
  } else{
    stop("The dimentions of the files are not equal!")
  }
}

