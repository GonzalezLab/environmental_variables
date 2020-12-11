### XtX analysis for 11 NA populations (Machado 2019)
## Created 1st September 2020

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("qvalue")
library(qvalue)


####### NORTH AMERICA ####### 
## Autosomes ##
xtx.pi.NA = read.table("/Users/pogo/Documents/Maria/Baypass/NA/analysis/core/xtx/machadoAUTO_clean_concatenated_summary_pi_xtx.out", h=F)
NA.names=read.table("/Users/pogo/Documents/Maria/Baypass/NA/baypass_data_def/split_files/map/all_concatenated_machadoTEs.map.names")
xtx.pi.NA$CHROM = NA.names$V1
xtx.pi.NA$POS = NA.names$V2
xtx.pi.NA$TE = NA.names$V3
NA.xtx.star = read.table("/Users/pogo/Documents/Maria/Baypass/NA/analysis/core/xtx/xtx.star.NA.11pop.txt",h=T)
# Get MAF
NA.MAF = (0.5 - abs(0.5 -xtx.pi.NA$V2))
NA.xtx.star$MAF = NA.MAF
NA.xtx.star.MAF = NA.xtx.star[NA.xtx.star$MAF>=0.01,]
# Get p-values
pvalues.bilateral.NA=1-(2*(abs(pchisq(NA.xtx.star$x, df=11)-0.5)))
pvalues.bilateral.NA.MAF=1-(2*(abs(pchisq(NA.xtx.star.MAF$x, df=11)-0.5)))
# Get q-values
qvalues.bilateral.NA = qvalue(pvalues.bilateral.NA)
# Create final matrix
xtx.pi.NA$MAF = NA.MAF
xtx.pi.NA$xtxstar = NA.xtx.star$x
xtx.pi.NA$pvalues = pvalues.bilateral.NA
xtx.pi.NA$qvalues = qvalues.bilateral.NA$qvalues
xtx.pi.NA.MAF = xtx.pi.NA[xtx.pi.NA$MAF>=0.01,]

write.table(xtx.pi.NA.MAF, "/Users/pogo/Documents/Maria/Baypass/NA/analysis/core/xtx/NA.xtx.MAF.allSNPs.txt", sep = "\t", row.names= F, quote=F)

num_sig_NA_auto = length(xtx.pi.NA.MAF$V1)*0.0005
significant.xtx.pi.NA.MAF = tail(xtx.pi.NA.MAF[order(xtx.pi.NA.MAF$xtxstar),], 1025L)
write.table(significant.xtx.pi.NA.MAF, "/Users/pogo/Documents/Maria/Baypass/NA/analysis/core/xtx/NA.xtx.MAF.significant0.005.txt", sep = "\t", row.names= F, quote=F)

# Write the bedfile:
bedfile.NA.MAF = significant.xtx.pi.NA.MAF[,c(8,9,9,6,12)]
write.table(bedfile.NA.MAF, "/Users/pogo/Documents/Maria/Baypass/NA/analysis/core/xtx/NA.xtx.MAF.significant0.005.bed", sep = "\t", row.names= F, quote=F)


# P-value histograms:
png("/Users/pogo/Documents/Maria/Baypass/NA/analysis/core/xtx/NA.auto.plots.png")
par(mfrow=c(2,2))
hist(pvalues.bilateral.NA, breaks=100, freq=F, main="p-values without MAF filtering")
abline(h=1, col="red")
hist(pvalues.bilateral.NA.MAF, breaks=100, freq=F, main="p-values with MAF filtering")
abline(h=1, col="red")
plot(xtx.pi.NA$V6, xtx.pi.NA$xtxstar)
plot(xtx.pi.NA.MAF$V6, xtx.pi.NA.MAF$xtxstar)
dev.off()


## X chromosome ##
xtx.pi.NA.Xchrom = read.table("/Users/pogo/Documents/Maria/Baypass/NA/analysis/core/xtx/machadoX_summary_pi_xtx.out", h=T)
NA.X.names=read.table("/Users/pogo/Documents/Maria/Baypass/NA/baypass_data_def/machadoX.NAMES.txt")
xtx.pi.NA.Xchrom$CHROM = NA.X.names$V1
xtx.pi.NA.Xchrom$POS = NA.X.names$V2
xtx.pi.NA.Xchrom$TE = NA.X.names$V3
NA.X.xtx.star = read.table("/Users/pogo/Documents/Maria/Baypass/NA/analysis/core/xtx/xtx.star.NA.Xchrom.11pop.txt",h=T)
# Get MAF
NA.X.MAF = (0.5 - abs(0.5 -xtx.pi.NA.Xchrom$M_P))
NA.X.xtx.star$MAF = NA.X.MAF
NA.X.xtx.star.MAF = NA.X.xtx.star[NA.X.xtx.star$MAF>=0.01,]
# Get p-values
pvalues.bilateral.NA.X=1-(2*(abs(pchisq(NA.X.xtx.star$x, df=11)-0.5)))
pvalues.bilateral.NA.X.MAF=1-(2*(abs(pchisq(NA.X.xtx.star.MAF$x, df=11)-0.5)))
# Get q-values
qvalues.bilateral.NA.X = qvalue(pvalues.bilateral.NA.X)
# Create final matrix
xtx.pi.NA.Xchrom$MAF = NA.X.MAF
xtx.pi.NA.Xchrom$xtxstar = NA.X.xtx.star$x
xtx.pi.NA.Xchrom$pvalues = pvalues.bilateral.NA.X
xtx.pi.NA.Xchrom$qvalues = qvalues.bilateral.NA.X$qvalues
xtx.pi.NA.Xchrom.MAF = xtx.pi.NA.Xchrom[xtx.pi.NA.Xchrom$MAF>=0.01,]

write.table(xtx.pi.NA.Xchrom.MAF, "/Users/pogo/Documents/Maria/Baypass/NA/analysis/core/xtx/NA.Xchrom.xtx.MAF.allSNPs.txt", sep = "\t", row.names= F, quote=F)

num_sig_NA_X = length(xtx.pi.NA.Xchrom.MAF$MRK)*0.0005
significant.xtx.pi.NA.X.MAF = tail(xtx.pi.NA.Xchrom.MAF[order(xtx.pi.NA.Xchrom.MAF$xtxstar),], 139L)
write.table(significant.xtx.pi.NA.X.MAF, "/Users/pogo/Documents/Maria/Baypass/NA/analysis/core/xtx/NA.Xchrom.xtx.MAF.significant0.005.txt", sep = "\t", row.names= F, quote=F)

# Write the bedfile:
bedfile.NA.X.MAF = significant.xtx.pi.NA.X.MAF[,c(8,9,9,6,12)]
write.table(bedfile.NA.X.MAF, "/Users/pogo/Documents/Maria/Baypass/NA/analysis/core/xtx/NA.Xchrom.xtx.MAF.significant0.005.bed", sep = "\t", row.names= F, quote=F)


# P-value histograms:
png("/Users/pogo/Documents/Maria/Baypass/NA/analysis/core/xtx/NA.Xchrom.plots.png")
par(mfrow=c(2,2))
hist(pvalues.bilateral.NA.X, breaks=100, freq=F, main="p-values without MAF filtering")
abline(h=1, col="red")
hist(pvalues.bilateral.NA.X.MAF, breaks=100, freq=F, main="p-values with MAF filtering")
abline(h=1, col="red")
plot(xtx.pi.NA.Xchrom$M_XtX, xtx.pi.NA.Xchrom$xtxstar)
plot(xtx.pi.NA.Xchrom.MAF$M_XtX, xtx.pi.NA.Xchrom.MAF$xtxstar)
dev.off()


# #write.table(spring30.original.xtx.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/spring_auto_MAF.txt", sep = "\t", row.names= F, quote=F)
# spring30.original.xtx.MAF$qvalue = qvalues.spring30.MAF$qvalue
# write.table(spring30.original.xtx.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/without_inversions/spring_auto_MAF_qvalue.txt", sep = "\t", row.names= F, quote=F)
# 
# significant.spring30.xtx.MAF = spring30.original.xtx.MAF[spring30.original.xtx.MAF$qvalue<=0.0001,]
# # Write a table with all data:
# write.table(significant.spring30.xtx.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/spring30_auto_xtx_complete_table.txt", sep = "\t", row.names= F, quote=F)
# # Write the bedfile:
# bedfile.spring30.MAF = significant.spring30.xtx.MAF[,c(8,9,9,6,12)]
# write.table(bedfile.spring30.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/spring30_auto_xtx.bed", sep = "\t", row.names= F, quote=F)
# 
# # New threshold (17-01-20):
# num_sig_sp_auto = length(spring30.original.xtx.MAF$V1)*0.0005
# significant.spring30.xtx.MAF.new = tail(spring30.original.xtx.MAF[order(spring30.original.xtx.MAF$star),], 700L)
# #significant.spring30.xtx.MAF.new = spring30.original.xtx.MAF[spring30.original.xtx.MAF$V6>=27,]
# #significant.spring30.xtx.MAF.new = spring30.original.xtx.MAF[spring30.original.xtx.MAF$qvalue<=0.000001,]
# # Write a table with all data:
# write.table(significant.spring30.xtx.MAF.new, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/spring30/spring30_auto_xtx_complete_table_top0.05_def.txt", sep = "\t", row.names= F, quote=F)
# # Write the bedfile:
# bedfile.spring30.MAF.new = significant.spring30.xtx.MAF.new[,c(8,9,9,6,11,12)]
# write.table(bedfile.spring30.MAF.new, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/spring30_auto_xtx_top0.05.bed", sep = "\t", row.names= F, quote=F)
# # For TEs with q-value <= 0-05:
# # For TEs with top 0.1%:
# # Finally we are keeping the 5% 70000 lines.
# # Finally we are keeping the 2.5% 35000 lines.
# significant.TE.spring30.xtx.MAF = tail(spring30.original.xtx.MAF[order(spring30.original.xtx.MAF$star),], 35000L)
# write.table(significant.TE.spring30.xtx.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/TE_thresholds_top_0.1/spring30_auto_xtx_TEs_top2.5%.txt", sep = "\t", row.names= F, quote=F)
# 
# 
# # P-value histograms:
# png("/Users/pogo/Documents/Maria/Baypass/paper/figures/spring30.auto.plots.png")
# par(mfrow=c(2,2))
# hist(pvalues.bilateral.spring30, breaks=100, freq=F, main="p-values without MAF filtering")
# abline(h=1, col="red")
# hist(pvalues.bilateral.spring30.MAF, breaks=100, freq=F, main="p-values with MAF filtering")
# abline(h=1, col="red")
# plot(spring30.original.xtx$V6, spring30.original.xtx$star)
# plot(spring30.original.xtx.MAF$V6, spring30.original.xtx.MAF$star)
# dev.off()
# 
# ####### EUROPE SPRING X CHROM ####### 
# spring30.original.xtx.X = read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/spring30/XtX/star.xtx/spring30Xsoze_summary_pi_xtx.out", h=T)
# spring30.xtx.X.names = read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/spring/30samples_batch/spring30_droseu14_X_names_SNPs.txt")
# spring30.original.xtx.X$CHR = spring30.xtx.X.names$V1
# spring30.original.xtx.X$POS = spring30.xtx.X.names$V2
# spring30.star.xtx.X = read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/spring30/XtX/star.xtx/spring.star.xtx.X.txt", h=T)
# #plot(spring30.original.xtx.X$M_XtX, spring30.star.xtx.X$x)
# # Get MAF
# MAF.spring.X = (0.5 - abs(0.5 -spring30.original.xtx.X$M_P))
# spring30.star.xtx.X$MAF = MAF.spring.X
# spring30.star.xtx.X.MAF = spring30.star.xtx.X[spring30.star.xtx.X$MAF>=0.01,]
# # Get p-value
# pvalues.bilateral.spring.X=1-(2*(abs(pchisq(spring30.star.xtx.X$x, df=14)-0.5)))
# pvalues.bilateral.spring.X.MAF=1-(2*(abs(pchisq(spring30.star.xtx.X.MAF$x, df=14)-0.5)))
# # Get q-value
# qvalues.spring.X.MAF = qvalue(pvalues.bilateral.spring.X.MAF)
# summary(qvalues.spring.X.MAF)
# # Create final matrix)
# spring30.original.xtx.X$MAF = MAF.spring.X
# spring30.original.xtx.X$star = spring30.star.xtx.X$x
# spring30.original.xtx.X.MAF = spring30.original.xtx.X[spring30.original.xtx.X$MAF>=0.01,]
# #write.table(spring30.original.xtx.X.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/spring_Xchr_MAF.txt", sep = "\t", row.names= F, quote=F)
# spring30.original.xtx.X.MAF$qvalues = qvalues.spring.X.MAF$qvalues
# write.table(spring30.original.xtx.X.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/without_inversions/spring_Xchr_MAF_qvalue.txt", sep = "\t", row.names= F, quote=F)
# 
# significant.spring30.xtx.X = spring30.original.xtx.X.MAF[spring30.original.xtx.X.MAF$qvalues<=0.05,]
# # Threshold: 18.210040
# def.spring30.X.MAF=significant.spring30.xtx.X[significant.spring30.xtx.X$M_XtX>=18.21,]
# # Write a table with all data:
# write.table(def.spring30.X.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/spring30_Xchrom_xtx_complete_table.txt", sep = "\t", row.names= F, quote=F)
# # Write the bedfile:
# bedfile.spring30.X.MAF = def.spring30.X.MAF[,c(8,9,9,6,12)]
# write.table(bedfile.spring30.X.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/spring30_Xchrom_xtx.bed", sep = "\t", row.names= F, quote=F)
# 
# # New threshold (17-02-20):
# num_sig_sp_X = length(spring30.original.xtx.X.MAF$MRK)*0.0005
# significant.spring30.xtx.X.new = tail(spring30.original.xtx.X.MAF[order(spring30.original.xtx.X.MAF$star),], 52L)
# #significant.spring30.xtx.X.new = spring30.original.xtx.X.MAF[spring30.original.xtx.X.MAF$qvalues<=0.001,]
# #significant.spring30.xtx.X.new=significant.spring30.xtx.X[significant.spring30.xtx.X$M_XtX>=23,]
# # Write a table with all data:
# write.table(significant.spring30.xtx.X.new, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/spring30/spring30_Xchrom_xtx_complete_table_top0.05_def.txt", sep = "\t", row.names= F, quote=F)
# # Write the bedfile:
# bedfile.spring30.X.MAF.new = significant.spring30.xtx.X.new[,c(8,9,9,6,11,12)]
# write.table(bedfile.spring30.X.MAF.new, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/spring30_Xchrom_xtx_top0.05.bed", sep = "\t", row.names= F, quote=F)
# # For TEs with threshold q-value <= 0.05:
# # For TEs with threshold of top 0.1%:
# # Finally we are keeping the 5% 5200 lines.
# # Finally we are keeping the 2.5% 2600 lines.
# significant.TE.spring30.xtx.X = tail(spring30.original.xtx.X.MAF[order(spring30.original.xtx.X.MAF$star),], 2600L)
# write.table(significant.TE.spring30.xtx.X, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/TE_thresholds_top_0.1/spring30_Xchrom_TE_top2.5%.txt", sep = "\t", row.names= F, quote=F)
# 
# 
# # P-value histograms:
# png("/Users/pogo/Documents/Maria/Baypass/paper/figures/spring30.X.plots.png")
# par(mfrow=c(2,2))
# hist(pvalues.bilateral.spring.X, breaks=100, freq=F, main="p-values without MAF filtering")
# abline(h=1, col="red")
# hist(pvalues.bilateral.spring.X.MAF, breaks=100, freq=F, main="p-values with MAF filtering")
# abline(h=1, col="red")
# plot(spring30.original.xtx.X$M_XtX, spring30.original.xtx.X$star)
# plot(spring30.original.xtx.X.MAF$M_XtX, spring30.original.xtx.X.MAF$star)
# dev.off()


