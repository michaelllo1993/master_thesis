wd = getwd()
require(ggplot2)
require(scales)
args = commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  print("Wrong number of arguments passed!")
  stop()
}

#USAGE: RScprit name input_file_name

SPL = read.csv(args[1], sep = ",", header = F)
#Read csv file
number_of_orgs = ((dim(SPL)[2]) / 3)
OoI_Ortho_LSAAR_length_diff = list()
OoI_Ortho_LSAAR_data = list()
OoI_Ortho_SP_length_diff = list()
OoI_Ortho_SP_data = list()
lengthened_indices = list()
correlations = list()
myNames = c()
mapper = read.csv(args[2],header = F)
fullNames = c()
for (i in seq(1, number_of_orgs, by = 1)) {
  test_index = which(as.vector(SPL[, (i * 3 - 2)]) != "NULL")[1]
  if (is.na(test_index)) {
    test_index = "UNKNOWN"
  }
  myNames[i] = strsplit(as.character(SPL[test_index, (i * 3 - 2)]), split = "0")[[1]][1]
  fullNames[i] = as.character(mapper[which(mapper[,2] == myNames[i]),1])
}
for (org in seq(2, number_of_orgs, by = 1)) {
  OoI_Ortho_LSAAR_length_diff[[org]] = SPL[, 3] - SPL[, org * 3]
  #length(homo sapiens LSAAR) - length(othologous LSAAR) -> positive number = LSAAR got longer -> negative number = LSAAR got shorter
  OoI_Ortho_SP_length_diff[[org]] = SPL[, 2] - SPL[, (org * 3) - 1]
  #length(homo sapiens SP) - length(othologous SP) -> positive number = SP got longer -> negative number = SP got shorter
  lengthened_indices[[org]] = which(OoI_Ortho_LSAAR_length_diff[[org]] >
                                      0)
  OoI_Ortho_LSAAR_data[[org]] = OoI_Ortho_LSAAR_length_diff[[org]][lengthened_indices[[org]]]
  OoI_Ortho_SP_data[[org]] = OoI_Ortho_SP_length_diff[[org]][lengthened_indices[[org]]]
  correlations[[org]] = cor.test(OoI_Ortho_LSAAR_data[[org]], OoI_Ortho_SP_data[[org]], method = "pearson",alternative = "two.sided")
  if (!dir.exists(path = paste(wd, sep = "/", "results/SP","length_analysis_plots"))) {
    dir.create(path = paste(wd, sep = "/", "results/SP","length_analysis_plots"))
  }
  to_plot = data.frame(SAAR = OoI_Ortho_LSAAR_data[[org]], SignalPeptide = OoI_Ortho_SP_data[[org]])
  sp = ggplot(data = to_plot, mapping = aes(x = SAAR, y = SignalPeptide))
  finale = sp + geom_point(alpha = .5,position = position_jitter(width=.05, height=.05)) + ggtitle(
    paste(
      "Scatterplot of SAAR length differences and SP length differences\n",
      fullNames[1],
      " - ",
      fullNames[org]
    )
  ) + theme(plot.title = element_text(hjust = 0.5)) + scale_x_continuous(breaks = seq(0, max(to_plot$SAAR)+1, by = 1)) + xlab("SAAR length difference") + ylab("Signal peptide length difference")
  ggsave(
    file = paste(
      sep = "",
      wd,
      "/results/SP/length_analysis_plots/",
      fullNames[1],
      "_",
      fullNames[org],
      "_scatterplot.svg"
    ),
    plot = finale,
    width = 5,
    height = 5
  )
}
correlations[[1]] = "Not applicable"
names(correlations) = fullNames
result = matrix(NA,ncol = 2,nrow = length(correlations))
rownames(result) = fullNames
colnames(result) = c("Pearson's R","p-value, two sided test")
for(i in seq(2,length(correlations))){
  result[i,] = c(round(correlations[[i]]$estimate,digits = 4),round(correlations[[i]]$p.value,digits=4))
}

write.csv(result,file = paste("results/SP/",fullNames[1],"_correlations_sigp_length_analysis.csv",sep = ""))
