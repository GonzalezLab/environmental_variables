### Created on 31st Januar 2020
## Same xtx_square.R but only for the TE batch


library(qvalue)

nfiles=50 # Only when the input data was splitted in different files
npops=20 # To be changed depending on the number of populations

# freq=matrix(read.table("/homes/users/mbogaerts/scratch/Baypass/DrosEU/DrosEU_2014_V3/all30/final/core/noinvAll30_1_summary_yij_pij.out",h=T)$M_Pstd,ncol=npops,byrow=TRUE)
# for(i in 2:nfiles){
#   freq=rbind(freq,matrix(read.table(paste0("/homes/users/mbogaerts/scratch/Baypass/DrosEU/DrosEU_2014_V3/all30/final/core/noinvAll30_",i,"_summary_yij_pij.out",sep =""),h=T)$M_Pstd,ncol=npops, byrow=TRUE))
# }
# mu=mean(as.numeric(freq))
# sd=sd(as.numeric(freq))
# freq=(freq-mu)/sd
# xtx.star=rowSums(freq**2)


### All30 autosomes only TEs
all30.freq=matrix(read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/only_TEs/all30/all30coreTEs_summary_yij_pij.out",h=T)$M_Pstd,ncol=npops,byrow=TRUE)
all30.mu=mean(as.numeric(all30.freq))
all30.sd=sd(as.numeric(all30.freq))
all30.freq=(all30.freq-all30.mu)/all30.sd
all30.xtx.star=rowSums(all30.freq**2)
all30.pi = read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/only_TEs/all30/all30coreTEs_summary_pi_xtx.out", h=T)
all30.TE.names = read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/all30/corrected_TEs/all30_droseu14_noinv_auto_only_TEs_names.txt", h=F)
all30.pi$TE = all30.TE.names$V1
all30.pi$POS = all30.TE.names$V2
# Get MAF:
all30.TE.MAF = (0.5 - abs(0.5 -all30.pi$M_P))
all30.pi$MAF = all30.TE.MAF
all30.pi$star = all30.xtx.star
all30.star.xtx.MAF = all30.pi[all30.pi$MAF>=0.01,]
# Get p-values
pvalues.bilateral.all30.TEs=1-(2*(abs(pchisq(all30.pi$M_XtX, df=20)-0.5)))
pvalues.bilateral.all30.TEs.MAF=1-(2*(abs(pchisq(all30.star.xtx.MAF$star, df=20)-0.5)))
# Get q-values
qvalues.all30.TEs = qvalue(pvalues.bilateral.all30.TEs)
qvalues.all30.TEs.MAF = qvalue(pvalues.bilateral.all30.TEs.MAF)
# Create final matrix
all30.TE.def = all30.star.xtx.MAF
all30.TE.def$qvalue = qvalues.all30.TEs.MAF$qvalues
significant.all30.xtx.TEs.MAF=all30.TE.def[all30.TE.def$qvalue<=0.05,]
# Plot
par(mfrow=c(2,2))
hist(pvalues.bilateral.all30.TEs, breaks=100, freq=F, main="p-values without MAF filtering")
abline(h=1, col="red")
hist(pvalues.bilateral.all30.TEs.MAF, breaks=100, freq=F, main="p-values with MAF filtering")
abline(h=1, col="red")
plot(all30.pi$M_XtX, all30.pi$star)
plot(all30.star.xtx.MAF$M_XtX, all30.star.xtx.MAF$star)

### All30 Xchr only TEs
npops=20
all30.X.freq=matrix(read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/only_TEs/all30/all30XnoinvsizeTEs_summary_yij_pij.out",h=T)$M_Pstd,ncol=npops,byrow=TRUE)
all30.X.mu=mean(as.numeric(all30.X.freq))
all30.X.sd=sd(as.numeric(all30.X.freq))
all30.X.freq=(all30.X.freq-all30.X.mu)/all30.X.sd
all30.X.xtx.star=rowSums(all30.X.freq**2)
all30.X.pi = read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/only_TEs/all30/all30XnoinvsizeTEs_summary_pi_xtx.out", h=T)
all30.X.TE.names = read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/all30/corrected_TEs/all30_droseu14_noinv_X_only_TEs_names.txt", h=F)
all30.X.pi$TE = all30.X.TE.names$V1
all30.X.pi$POS = all30.X.TE.names$V2
# # Get MAF:
# all30.TE.MAF = (0.5 - abs(0.5 -all30.pi$M_P))
# all30.pi$MAF = all30.TE.MAF
# all30.pi$star = all30.xtx.star
# all30.star.xtx.MAF = all30.pi[all30.pi$MAF>=0.01,]
# # Get p-values
# pvalues.bilateral.all30.TEs=1-(2*(abs(pchisq(all30.pi$M_XtX, df=20)-0.5)))
# pvalues.bilateral.all30.TEs.MAF=1-(2*(abs(pchisq(all30.star.xtx.MAF$star, df=20)-0.5)))
# # Get q-values
# qvalues.all30.TEs = qvalue(pvalues.bilateral.all30.TEs)
# qvalues.all30.TEs.MAF = qvalue(pvalues.bilateral.all30.TEs.MAF)
# # Create final matrix
# all30.TE.def = all30.star.xtx.MAF
# all30.TE.def$qvalue = qvalues.all30.TEs.MAF$qvalues
# significant.all30.xtx.TEs.MAF=all30.TE.def[all30.TE.def$qvalue<=0.05,]
# # Plot
# par(mfrow=c(2,2))
# hist(pvalues.bilateral.all30.TEs, breaks=100, freq=F, main="p-values without MAF filtering")
# abline(h=1, col="red")
# hist(pvalues.bilateral.all30.TEs.MAF, breaks=100, freq=F, main="p-values with MAF filtering")
# abline(h=1, col="red")
# plot(all30.pi$M_XtX, all30.pi$star)
# plot(all30.star.xtx.MAF$M_XtX, all30.star.xtx.MAF$star)


dev.off()