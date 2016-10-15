library(seqinr)
SPL=read.csv('SPL_length_Hsap_Ortho.csv',sep=",");#Read csv file
SP=read.csv('SP_length_Hsap_Ortho.csv',sep=",")

#ANALYSIS OF LENGTH OF SIGNAL PEPTIDES WITH AT LEAST 5 LEUCINE REPEATS
Hsap_Ptro_SPL= SPL[,2] - SPL[,5]; #length(homo sapiens sp) - length(othologous sp) -> positive number = sp got longer -> negative number = sp got shorter
Hsap_Mmus_SPL = SPL[,2] - SPL[,8];
Hsap_Btau_SPL = SPL[,2] - SPL[,11];
Hsap_Ggal_SPL = SPL[,2] - SPL[,14];
Hsap_Xtro_SPL = SPL[,2] - SPL[,17];

Hs_Pt_SPL=SPL[,3] - SPL[,6];
Hs_Pt.i_SPL=which(Hs_Pt_SPL>0)

Hs_Mm_SPL=SPL[,3] - SPL[,9]; #length of leucine runs in homo sapiens' sp - length of leucine runs in orthologous' sp
Hs_Mm.i_SPL=which(Hs_Mm_SPL>0) #Which l-runs got longer 

Hs_Bt_SPL=SPL[,3] - SPL[,12];
Hs_Bt.i_SPL=which(Hs_Bt_SPL>0)

Hs_Gg_SPL=SPL[,3] - SPL[,15];
Hs_Gg.i_SPL=which(Hs_Gg_SPL>0)

Hs_Xt_SPL=SPL[,3] - SPL[,18];
Hs_Xt.i_SPL=which(Hs_Xt_SPL>0)

Hsap_Mmus1_SPL=Hsap_Mmus_SPL[Hs_Mm.i_SPL];#Signal pepdtide length differences only where l-run got longer
Hsap_Ptro1_SPL=Hsap_Ptro_SPL[Hs_Pt.i_SPL];
Hsap_Btau1_SPL=Hsap_Btau_SPL[Hs_Bt.i_SPL];
Hsap_Ggal1_SPL=Hsap_Ggal_SPL[Hs_Gg.i_SPL];
Hsap_Xtro1_SPL=Hsap_Xtro_SPL[Hs_Xt.i_SPL];

Hsap_Mmus2_SPL=Hsap_Mmus_SPL[which(!is.na(Hsap_Mmus_SPL))];
Hsap_Ptro2_SPL=Hsap_Ptro_SPL[which(!is.na(Hsap_Ptro_SPL))];
Hsap_Btau2_SPL=Hsap_Btau_SPL[which(!is.na(Hsap_Btau_SPL))];
Hsap_Ggal2_SPL=Hsap_Ggal_SPL[which(!is.na(Hsap_Ggal_SPL))];
Hsap_Xtro2_SPL=Hsap_Xtro_SPL[which(!is.na(Hsap_Xtro_SPL))];

#Correlation coefficients (both Spearman & Pearson methods)
ce_spearman_human_mouse=cor(Hsap_Mmus_SPL[which(!is.na(Hsap_Mmus_SPL))],Hs_Mm_SPL[which(!is.na(Hs_Mm_SPL))],method = c("spearman"))
ce_pearson_human_mouse=cor(Hsap_Mmus_SPL[which(!is.na(Hsap_Mmus_SPL))],Hs_Mm_SPL[which(!is.na(Hs_Mm_SPL))],method = c("pearson"))
ce_spearman_human_chimp=cor(Hsap_Ptro_SPL[which(!is.na(Hsap_Ptro_SPL))],Hs_Pt_SPL[which(!is.na(Hs_Pt_SPL))],method = c("spearman"));
ce_pearson_human_chimp=cor(Hsap_Ptro_SPL[which(!is.na(Hsap_Ptro_SPL))],Hs_Pt_SPL[which(!is.na(Hs_Pt_SPL))],method = c("pearson"));
ce_spearman_human_cow=cor(Hsap_Btau_SPL[which(!is.na(Hsap_Btau_SPL))],Hs_Bt_SPL[which(!is.na(Hs_Bt_SPL))],method = c("spearman"));
ce_pearson_human_cow=cor(Hsap_Btau_SPL[which(!is.na(Hsap_Btau_SPL))],Hs_Bt_SPL[which(!is.na(Hs_Bt_SPL))],method = c("pearson"));
ce_spearman_human_chicken=cor(Hsap_Ggal_SPL[which(!is.na(Hsap_Ggal_SPL))],Hs_Gg_SPL[which(!is.na(Hs_Gg_SPL))],method = c("spearman"));
ce_pearson_human_chicken=cor(Hsap_Ggal_SPL[which(!is.na(Hsap_Ggal_SPL))],Hs_Gg_SPL[which(!is.na(Hs_Gg_SPL))],method = c("pearson"));
ce_spearman_human_frog=cor(Hsap_Xtro_SPL[which(!is.na(Hsap_Xtro_SPL))],Hs_Xt_SPL[which(!is.na(Hs_Xt_SPL))],method = c("spearman"));
ce_pearson_human_frog=cor(Hsap_Xtro_SPL[which(!is.na(Hsap_Xtro_SPL))],Hs_Xt_SPL[which(!is.na(Hs_Xt_SPL))],method = c("pearson"))



#ANALYSIS OF LENGTH OF ALL SIGNAL PEPTIDES
Hsap_Ptro_SP= SP[,2] - SP[,4]; #length(homo sapiens sp) - length(othologous sp) -> positive number = sp got longer -> negative number = sp got shorter
Hsap_Mmus_SP = SP[,2] - SP[,6];
Hsap_Btau_SP = SP[,2] - SP[,8];
Hsap_Ggal_SP = SP[,2] - SP[,10];
Hsap_Xtro_SP = SP[,2] - SP[,12];

Hsap_Mmus_SP = Hsap_Mmus_SP[which(!is.na(Hsap_Mmus_SP))];
Hsap_Btau_SP = Hsap_Btau_SP[which(!is.na(Hsap_Btau_SP))];
Hsap_Ptro_SP = Hsap_Ptro_SP[which(!is.na(Hsap_Ptro_SP))];
Hsap_Ggal_SP = Hsap_Ggal_SP[which(!is.na(Hsap_Ggal_SP))];
Hsap_Xtro_SP = Hsap_Xtro_SP[which(!is.na(Hsap_Xtro_SP))];

ALL_splengthsdiffs_SPL = c(Hsap_Ptro_SPL,Hsap_Mmus_SPL,Hsap_Btau_SPL,Hsap_Ggal_SPL,Hsap_Xtro_SPL);
ALL_llengthsdiffs_SPL = c(Hs_Pt_SPL,Hs_Mm_SPL,Hs_Bt_SPL,Hs_Gg_SPL,Hs_Xt_SPL);
ALL_SPL=c(Hsap_Mmus1_SPL,Hsap_Ptro1_SPL,Hsap_Btau1_SPL,Hsap_Ggal1_SPL,Hsap_Xtro1_SPL);
ALL2_SPL=c(Hsap_Mmus2_SPL,Hsap_Ptro2_SPL,Hsap_Btau2_SPL,Hsap_Ggal2_SPL,Hsap_Xtro2_SPL);
ALL_SP = c(Hsap_Mmus_SP,Hsap_Btau_SP,Hsap_Ptro_SP,Hsap_Ggal_SP,Hsap_Xtro_SP);

#Set of plots Human - ALL
par(mfrow=c(2,2))
plot(density(ALL_SP), main="Signal peptide length difference (Human - orhologues); all cases (SP)");
plot(density(ALL2_SPL,bw=0.5),main="Signal peptide length difference (Human - orhologues); all cases (SPL)");
plot(density(ALL_SPL),main="Signal peptide length difference (Human - orhologues); leucine gain cases (SPL)");
#plot(ALL_splengthsdiffs_SPL,ALL_llengthsdiffs_SPL,main="Scatterplot ALL",pch=20);
smoothScatter(ALL_splengthsdiffs_SPL,ALL_llengthsdiffs_SPL,main="Scatterplot Human - ALL",xlab="SP length differenece",ylab="leucine run length difference",nrpoints=0)


#Set of plots Human-Mouse
par(mfrow=c(2,2))
plot(density(Hsap_Mmus_SP,bw=0.5),main="Signal peptide length difference (Human - Mouse); all cases (SP)")
plot(density(Hsap_Mmus2_SPL,bw=0.5),main="Signal peptide length difference (Human - Mouse); all cases (SPL)")
plot(density(Hsap_Mmus1_SPL),main="Signal peptide length difference (Human - Mouse); leucine gain cases (SPL)")
#plot(Hsap_Mmus_SPL[which(!is.na(Hsap_Mmus_SPL))],Hs_Mm_SPL[which(!is.na(Hs_Mm_SPL))],main="Scatterplot Homo sapiens - Mus musculus",pch=20)
smoothScatter(Hsap_Mmus_SPL[which(!is.na(Hsap_Mmus_SPL))],Hs_Mm_SPL[which(!is.na(Hs_Mm_SPL))],main=bquote(atop(Scatterplot~Human~-~Mouse,CE~pearson==.(ce_pearson_human_mouse)~~CE~Spearman==.(ce_spearman_human_mouse))),xlab="SP length differenece",ylab="leucine run length difference",nrpoints=0)

#Set of plots Human-Chimp
par(mfrow=c(2,2))
plot(density(Hsap_Ptro_SP),main="Signal peptide length difference (Human - Chimp); all cases (SP)")
plot(density(Hsap_Ptro2_SPL),main="Signal peptide length difference (Human - Chimp); all cases (SPL)")
plot(density(Hsap_Ptro1_SPL),main="Signal peptide length difference (Human - Chimp); leucine gain cases (SPL)")
#plot(Hsap_Ptro_SPL[which(!is.na(Hsap_Ptro_SPL))],Hs_Pt_SPL[which(!is.na(Hs_Pt_SPL))],main="Scatterplot Homo sapiens - Pan troglodytes",pch=20)
smoothScatter(Hsap_Ptro1_SPL[which(!is.na(Hsap_Ptro_SPL))], Hs_Pt_SPL[which(!is.na(Hs_Pt_SPL))],main=bquote(atop(Scatterplot~Human~-~Chimp,CE~pearson==.(ce_pearson_human_chimp)~~CE~Spearman==.(ce_spearman_human_chimp))),xlab="SP length differenece",ylab="leucine run length difference",nrpoints=0,xlim=c(-3,5), ylim=c(-3,5))

#Set of plots Human-Cow
par(mfrow=c(2,2))
plot(density(Hsap_Btau_SP),main="Signal peptide length difference (Human - Cow); all cases (SP)")
plot(density(Hsap_Btau2_SPL),main="Signal peptide length difference (Human - Cow); all cases (SPL)")
plot(density(Hsap_Btau1_SPL),main="Signal peptide length difference (Human - Cow); leucine gain cases (SPL)")
#plot(Hsap_Btau_SPL[which(!is.na(Hsap_Btau_SPL))],Hs_Bt_SPL[which(!is.na(Hs_Bt_SPL))],main="Scatterplot Homo sapiens - Bos taurus",pch=20)
smoothScatter(Hsap_Btau_SPL[which(!is.na(Hsap_Btau_SPL))],Hs_Bt_SPL[which(!is.na(Hs_Bt_SPL))],main=bquote(atop(Scatterplot~Human~-~Cow,CE~pearson==.(ce_pearson_human_cow)~~CE~Spearman==.(ce_spearman_human_cow))),xlab="SP length differenece",ylab="leucine run length difference",nrpoints=0)

#Set of plots Human-Chicken
par(mfrow=c(2,2))
plot(density(Hsap_Ggal_SP),main="Signal peptide length difference (Human - Chicken); all cases (SP)")
plot(density(Hsap_Ggal2_SPL),main="Signal peptide length difference (Human - Chicken); all cases (SPL)")
plot(density(Hsap_Ggal1_SPL),main="Signal peptide length difference (Human - Chicken); leucine gain cases (SPL)");
#plot(Hsap_Ggal_SPL[which(!is.na(Hsap_Ggal_SPL))],Hs_Gg_SPL[which(!is.na(Hs_Gg_SPL))],main="Scatterplot Homo sapiens - Gallus gallus",pch=20)
smoothScatter(Hsap_Ggal_SPL[which(!is.na(Hsap_Ggal_SPL))],Hs_Gg_SPL[which(!is.na(Hs_Gg_SPL))],main=bquote(atop(Scatterplot~Human~-~Chicken,CE~pearson==.(ce_pearson_human_chicken)~~CE~Spearman==.(ce_spearman_human_chicken))),xlab="SP length differenece",ylab="leucine run length difference")

#Set of plots Human-Frog
par(mfrow=c(2,2))
plot(density(Hsap_Xtro_SP),main="Signal peptide length difference (Human - Frog); all cases (SP)")
plot(density(Hsap_Xtro2_SPL),main="Signal peptide length difference (Human - Frog); all cases (SPL)")
plot(density(Hsap_Xtro1_SPL),main="Signal peptide length difference (Human - Frog); leucine gain cases (SPL)")
#plot(Hsap_Xtro_SPL[which(!is.na(Hsap_Xtro_SPL))],Hs_Xt_SPL[which(!is.na(Hs_Xt_SPL))],main="Scatterplot Homo sapiens - Xenopus tropicalis",pch=20)
smoothScatter(Hsap_Xtro_SPL[which(!is.na(Hsap_Xtro_SPL))],Hs_Xt_SPL[which(!is.na(Hs_Xt_SPL))],main=bquote(atop(Scatterplot~Human~-~Frog,CE~pearson==.(ce_pearson_human_frog)~~CE~Spearman==.(ce_spearman_human_frog))),xlab="SP length differenece",ylab="leucine run length difference",nrpoints=0)


#TABLES DATA
#H.sapiens - P.troglodytes
sum(Hsap_Ptro_SPL[which(!is.na(Hsap_Ptro_SPL))]==0 & Hs_Pt_SPL[which(!is.na(Hs_Pt_SPL))]==0)#no of aa with no changes in sp length & leucine run length
sum(Hsap_Ptro_SPL[which(!is.na(Hsap_Ptro_SPL))]!=0 & Hs_Pt_SPL[which(!is.na(Hs_Pt_SPL))]!=0)#no of aa with changes in sp length & leucine run length
sum(Hsap_Ptro_SPL[which(!is.na(Hsap_Ptro_SPL))]==0 & Hs_Pt_SPL[which(!is.na(Hs_Pt_SPL))]!=0)#no of aa with changes in leucine runs length but with no changes in sp length
sum(Hsap_Ptro_SPL[which(!is.na(Hsap_Ptro_SPL))]!=0 & Hs_Pt_SPL[which(!is.na(Hs_Pt_SPL))]==0)#no of aa with changes in sp length but with no changes in leucine runs lengths
length(Hsap_Ptro_SPL[which(!is.na(Hsap_Ptro_SPL))])#no of all sp changes

#H.sapiens - M.musculus
sum(Hsap_Mmus_SPL[which(!is.na(Hsap_Mmus_SPL))]==0 & Hs_Mm_SPL[which(!is.na(Hs_Mm_SPL))]==0)#no of aa with no changes in sp length & leucine run length
sum(Hsap_Mmus_SPL[which(!is.na(Hsap_Mmus_SPL))]!=0 & Hs_Mm_SPL[which(!is.na(Hs_Mm_SPL))]!=0)#no of aa with changes in sp length & leucine run length
sum(Hsap_Mmus_SPL[which(!is.na(Hsap_Mmus_SPL))]==0 & Hs_Mm_SPL[which(!is.na(Hs_Mm_SPL))]!=0)#no of aa with changes in leucine runs length but with no changes in sp length
sum(Hsap_Mmus_SPL[which(!is.na(Hsap_Mmus_SPL))]!=0 & Hs_Mm_SPL[which(!is.na(Hs_Mm_SPL))]==0)#no of aa with changes in sp length but with no changes in leucine runs lengths
length(Hsap_Mmus_SPL[which(!is.na(Hsap_Mmus_SPL))])#no of all sp changes

#H.sapiens - B.taurus
sum(Hsap_Btau_SPL[which(!is.na(Hsap_Btau_SPL))]==0 & Hs_Bt_SPL[which(!is.na(Hs_Bt_SPL))]==0)#no of aa with no changes in sp length & leucine run length
sum(Hsap_Btau_SPL[which(!is.na(Hsap_Btau_SPL))]!=0 & Hs_Bt_SPL[which(!is.na(Hs_Bt_SPL))]!=0)#no of aa with changes in sp length & leucine run length
sum(Hsap_Btau_SPL[which(!is.na(Hsap_Btau_SPL))]==0 & Hs_Bt_SPL[which(!is.na(Hs_Bt_SPL))]!=0)#no of aa with changes in leucine runs length but with no changes in sp length
sum(Hsap_Btau_SPL[which(!is.na(Hsap_Btau_SPL))]!=0 & Hs_Bt_SPL[which(!is.na(Hs_Bt_SPL))]==0)#no of aa with changes in sp length but with no changes in leucine runs lengths
length(Hsap_Btau_SPL[which(!is.na(Hsap_Btau_SPL))])#no of all sp changes

#H.sapiens - G.gallus
sum(Hsap_Ggal_SPL[which(!is.na(Hsap_Ggal_SPL))]==0 & Hs_Gg_SPL[which(!is.na(Hs_Gg_SPL))]==0)#no of aa with no changes in sp length & leucine run length
sum(Hsap_Ggal_SPL[which(!is.na(Hsap_Ggal_SPL))]!=0 & Hs_Gg_SPL[which(!is.na(Hs_Gg_SPL))]!=0)#no of aa with changes in sp length & leucine run length
sum(Hsap_Ggal_SPL[which(!is.na(Hsap_Ggal_SPL))]==0 & Hs_Gg_SPL[which(!is.na(Hs_Gg_SPL))]!=0)#no of aa with changes in leucine runs length but with no changes in sp length
sum(Hsap_Ggal_SPL[which(!is.na(Hsap_Ggal_SPL))]!=0 & Hs_Gg_SPL[which(!is.na(Hs_Gg_SPL))]==0)#no of aa with changes in sp length but with no changes in leucine runs lengths
length(Hsap_Ggal_SPL[which(!is.na(Hsap_Ggal_SPL))])#no of all sp changes

#H.sapiens - X.tropicalis
sum(Hsap_Xtro_SPL[which(!is.na(Hsap_Xtro_SPL))]==0 & Hs_Xt_SPL[which(!is.na(Hs_Xt_SPL))]==0)#no of aa with no changes in sp length & leucine run length
sum(Hsap_Xtro_SPL[which(!is.na(Hsap_Xtro_SPL))]!=0 & Hs_Xt_SPL[which(!is.na(Hs_Xt_SPL))]!=0)#no of aa with changes in sp length & leucine run length
sum(Hsap_Xtro_SPL[which(!is.na(Hsap_Xtro_SPL))]==0 & Hs_Xt_SPL[which(!is.na(Hs_Xt_SPL))]!=0)#no of aa with changes in leucine runs length but with no changes in sp length
sum(Hsap_Xtro_SPL[which(!is.na(Hsap_Xtro_SPL))]!=0 & Hs_Xt_SPL[which(!is.na(Hs_Xt_SPL))]==0)#no of aa with changes in sp length but with no changes in leucine runs lengths
length(Hsap_Xtro_SPL[which(!is.na(Hsap_Xtro_SPL))])#no of all sp changes
