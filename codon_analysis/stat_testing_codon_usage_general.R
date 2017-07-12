list_test=list()
shapiro_results = c()

for(i in 1:dim(codons_count_protein_fraction)[1]){
  list_test[[i]] = codons_count_protein_fraction[i,]
}
names(list_test) = rownames(codons_count_protein_fraction)

# Testing for normality
for(i in 1:length(list_test)){
  shapiro_results[i]=shapiro.test(list_test[[i]])$p.value
}

percentages = c()
names = percentages
df_codon_usage = as.data.frame(t(codons_count_protein_fraction))
for(i in 1:dim(df_codon_usage)[2]){
  percentages = append(percentages,df_codon_usage[,i])
  names = append(names,rep(colnames(df_codon_usage)[i],6))
}
df_test = data.frame(percentages, names)
# Testing for homoscedasticity
levenes_test_result = leveneTest(
  df_test$percentages,group = as.factor(df_test$names),
  location = "mean"
)

# Testing for differences in means (ANOVA)
a=aov(df_test$percentages ~ df_test$names)
summary(a)

