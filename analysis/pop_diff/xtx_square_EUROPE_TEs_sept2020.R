### Created on 29th November 2019
## By Mathieu Gautier
## For the for loop, use it in the cluster with interactive session due to the large amount of data
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("qvalue")
library(qvalue)

nfiles=50 # Only when the input data was splitted in different files
npops=20 # To be changed depending on the number of populations

freq=matrix(read.table("/homes/users/mbogaerts/scratch/Baypass/DrosEU/DrosEU_2014_V3/all30/final/core/noinvAll30_1_summary_yij_pij.out",h=T)$M_Pstd,ncol=npops,byrow=TRUE)
for(i in 2:nfiles){
  freq=rbind(freq,matrix(read.table(paste0("/homes/users/mbogaerts/scratch/Baypass/DrosEU/DrosEU_2014_V3/all30/final/core/noinvAll30_",i,"_summary_yij_pij.out",sep =""),h=T)$M_Pstd,ncol=npops, byrow=TRUE))
}
mu=mean(as.numeric(freq))
sd=sd(as.numeric(freq))
freq=(freq-mu)/sd
xtx.star=rowSums(freq**2)


####### EUROPE ALL AUTOSOMES ####### 
all30.xtx=read.table("/Volumes/Maxtor/Maria/Baypass/13-02-20/all30/new/XtX/xtx_star/all30_concatenated_summary_pi_xtx.out", h=F)
names=read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/all30/corrected_TEs/splitted_geno/map/all30_droseu14_auto.geno.names.perfiles.txt")
all30.xtx$V8=names$V1
all30.xtx$V9=names$V2
xtx.star=read.table("/Volumes/Maxtor/Maria/Baypass/13-02-20/all30/new/XtX/xtx_star/xtx.star.txt",h=T)
# Get MAF
MAF = (0.5 - abs(0.5 -all30.xtx$V2))
xtx.star$MAF = MAF
xtx.star.MAF = xtx.star[xtx.star$MAF>=0.01,]
# Compute p-value
pvalues.bilateral=1-(2*(abs(pchisq(xtx.star$x, df=20)-0.5)))
pvalues.bilateral.MAF=1-(2*(abs(pchisq(xtx.star.MAF$x, df=20)-0.5)))
# Compute q-value
qvalues.xtx.MAF= qvalue(pvalues.bilateral.MAF)
summary(qvalues.xtx.MAF)
# Create final matrix
all30.xtx$V10=MAF
all30.xtx$V11=xtx.star$x
all30.xtx.MAF=all30.xtx[all30.xtx$V10>=0.01,]
#write.table(all30.xtx.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/all_auto_MAF.txt", sep = "\t", row.names= F, quote=F)
all30.xtx.MAF$qvalue=qvalues.xtx.MAF$qvalues
#write.table(all30.xtx.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/all_auto_MAF_qvalue.txt", sep = "\t", row.names= F, quote=F)

## Without inversions
#write.table(all30.xtx.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/without_inversions/all_auto_MAF.txt", sep = "\t", row.names= F, quote=F)
#significant.all30.xtx.MAF=all30.xtx.MAF[all30.xtx.MAF$qvalue<=0.000001,]
# Write a table with all data:
#write.table(significant.all30.xtx.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/all30_autosomes_xtx_complete_table.txt", sep = "\t", row.names= F, quote=F)
# Write the bedfile:
#bedfile.all30.MAF = significant.all30.xtx.MAF[,c(8,9,9,6,12)]
#write.table(bedfile.all30.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/all30_autosomes_xtx.bed", sep = "\t", row.names= F, quote=F)

# New threshold (17-01-20):
# Threshold for the top 0.05% (0.0005)
#num_sig = length(all30.xtx.MAF$V1)*0.0005
#significant.all30.xtx.MAF.new=tail(all30.xtx.MAF[order(all30.xtx.MAF$V11),], 665L)
#significant.all30.xtx.MAF.new=all30.xtx.MAF[all30.xtx.MAF$V6>=35,]
#significant.all30.xtx.MAF.new=all30.xtx.MAF[all30.xtx.MAF$qvalue<=2.935491e-08,]
# Write a table with all data:
#write.table(significant.all30.xtx.MAF.new, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/all30_autosomes_xtx_complete_table_top0.005_prueba2.txt", sep = "\t", row.names= F, quote=F)
# Write the bedfile:
#bedfile.all30.MAF.new = significant.all30.xtx.MAF.new[,c(8,9,9,6,11,12)]
#write.table(bedfile.all30.MAF.new, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/all30_autosomes_xtx_top0.05_prueba.bed", sep = "\t", row.names= F, quote=F)
# For TEs using a threshold of q-value <= 0.05:
# For TEs using a threshold of top 2.5%:
# Finally we are keeping the 5% 66500 lines.
significant.TE.all30.xtx.MAF=tail(all30.xtx.MAF[order(all30.xtx.MAF$V11),], 33275L)
write.table(significant.TE.all30.xtx.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_sept_2020_std_corrected_variables/all30/xtx/TEs/all30_autosomes_xtx_TEs_top_2.5%.txt", sep = "\t", row.names= F, quote=F)
significant.TE.all30.xtx.MAF=tail(all30.xtx.MAF[order(all30.xtx.MAF$V11),], 13310L)
write.table(significant.TE.all30.xtx.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_sept_2020_std_corrected_variables/all30/xtx/TEs/all30_autosomes_xtx_TEs_top_1%.txt", sep = "\t", row.names= F, quote=F)





####### EUROPE ALL CHROMOSOME X ####### 
# Done in the cluster
all30.freq.X=matrix(read.table("/homes/users/mbogaerts/scratch/Baypass/DrosEU/DrosEU_2014_V3/all30/final/core/all30Xnoinvsize_summary_yij_pij.out",h=T)$M_Pstd,ncol=npops,byrow=TRUE)
all30.mu.X=mean(as.numeric(all30.freq.X))
all30.sd.X=sd(as.numeric(all30.freq.X))
all30.freq.X=(all30.freq.X-all30.mu.X)/all30.sd.X
xtx.star.X=rowSums(all30.freq.X**2)
# In local
all30.xtx.X = read.table("/Volumes/Maxtor/Maria/Baypass/13-02-20/all30/new/XtX/xtx_star/all30Xnoinvsize_summary_pi_xtx.out", h=T)
all30.X.names=read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/all30/corrected_TEs/all30_droseu14_X_SNPs_names.txt", h=F)
all30.xtx.X$CHROM = all30.X.names$V1
all30.xtx.X$POS = all30.X.names$V2
all30.star.xtx.X = read.table("/Volumes/Maxtor/Maria/Baypass/13-02-20/all30/new/XtX/xtx_star/xtx.star.X.txt", h=F)
#plot(all30.xtx.X$M_XtX, all30.star.xtx.X$V2)
# Get MAF
all30.MAF.X = (0.5 - abs(0.5 -all30.xtx.X$M_P))
all30.star.xtx.X$V3 = all30.MAF.X
all30.star.xtx.X.MAF = all30.star.xtx.X[all30.star.xtx.X$V3>=0.01,]
# Get p-values
pvalues.bilateral.all30.X=1-(2*(abs(pchisq(all30.star.xtx.X$V2, df=20)-0.5)))
pvalues.bilateral.all30.X.MAF=1-(2*(abs(pchisq(all30.star.xtx.X.MAF$V2, df=20)-0.5)))
# Get q-values
qvalues.all30.X = qvalue(pvalues.bilateral.all30.X)
qvalues.all30.X.MAF = qvalue(pvalues.bilateral.all30.X.MAF)
# Create final matrix
all30.xtx.X$MAF = all30.MAF.X
all30.xtx.X$star = all30.star.xtx.X$V2
all30.xtx.X.MAF = all30.xtx.X[all30.xtx.X$MAF>0.01,]
#write.table(all30.xtx.X.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/all_xchr_MAF.txt", sep = "\t", row.names= F, quote=F)
all30.xtx.X.MAF$qvalue = qvalues.all30.X.MAF$qvalue
#write.table(all30.xtx.X.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/all_xchr_MAF_qvalue.txt", sep = "\t", row.names= F, quote=F)

#significant.all30.xtx.X.MAF=all30.xtx.X.MAF[all30.xtx.X.MAF$qvalue<=0.05,]
# Threshold: 22.98934
#def.all30.X.MAF=significant.all30.xtx.X.MAF[significant.all30.xtx.X.MAF$M_XtX>=22.98,]
# Write a table with all data:
#write.table(def.all30.X.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/all30_Xchr_xtx_complete_table.txt", sep = "\t", row.names= F, quote=F)
# Write the bedfile:
#bedfile.all30.X.MAF = def.all30.X.MAF[,c(8,9,9,6,12)]
#write.table(bedfile.all30.X.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/all30_Xchr_xtx.bed", sep = "\t", row.names= F, quote=F)

# New threshold (17-01-20):
num_sig_allX = length(all30.xtx.X.MAF$MRK)*0.0005
significant.all30.xtx.X.MAF.new=tail(all30.xtx.X.MAF[order(all30.xtx.X.MAF$star),], 54L)
#significant.all30.xtx.X.MAF.new=all30.xtx.X.MAF[all30.xtx.X.MAF$qvalue<=0.0005687158,]
#significant.all30.xtx.X.MAF.new=all30.xtx.X.MAF[all30.xtx.X.MAF$M_XtX>=30,]
# Write a table with all data:
#write.table(significant.all30.xtx.X.MAF.new, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/all30_Xchr_xtx_complete_table_top0.05_def.txt", sep = "\t", row.names= F, quote=F)
# Write the bedfile:
#bedfile.all30.X.MAF.new = significant.all30.xtx.X.MAF.new[,c(8,9,9,6,11,12)]
#write.table(bedfile.all30.X.MAF.new, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/all30_Xchr_xtx_top0.05.bed", sep = "\t", row.names= F, quote=F)
# For TEs using a threshold of q-value <= 0.05:
# For TEs using a threshold of top 0.1%:
# Finally we are keeping the 5% 5400 lines.
# Finally we are keeping the 2.5% 2700 lines.
significant.TE.all30.xtx.X.MAF=tail(all30.xtx.X.MAF[order(all30.xtx.X.MAF$star),], 2700L)
write.table(significant.TE.all30.xtx.X.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_sept_2020_std_corrected_variables/all30/xtx/TEs/all30_Xchr_xtx_TEs_top2.5%.txt", sep = "\t", row.names= F, quote=F)
significant.TE.all30.xtx.X.MAF=tail(all30.xtx.X.MAF[order(all30.xtx.X.MAF$star),], 1072L)
write.table(significant.TE.all30.xtx.X.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_sept_2020_std_corrected_variables/all30/xtx/TEs/all30_Xchr_xtx_TEs_top1%.txt", sep = "\t", row.names= F, quote=F)


####### EUROPE SPRING AUTOSOMES ####### 
spring30.original.xtx=read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/spring30/XtX/star.xtx/spring30_all_concatenated_summary_pi_xtx.txt", h=F)
spring30.names=read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/spring/30samples_batch/splitted_autosomes/spring30_droseu14_auto.map.NAMES.ALL.txt")
spring30.original.xtx$CHROM = spring30.names$V1
spring30.original.xtx$POS = spring30.names$V2
spring30.xtx.star=read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/spring30/XtX/star.xtx/spring30.xtx.star.txt",h=T)
#plot(spring30.original.xtx$V6,spring30.xtx.star$x)
# Get MAF
spring30.MAF = (0.5 - abs(0.5 -spring30.original.xtx$V2))
spring30.xtx.star$MAF = spring30.MAF
spring30.xtx.star.MAF = spring30.xtx.star[spring30.xtx.star$MAF>=0.01,]
# Get p-values
pvalues.bilateral.spring30=1-(2*(abs(pchisq(spring30.xtx.star$x, df=14)-0.5)))
pvalues.bilateral.spring30.MAF=1-(2*(abs(pchisq(spring30.xtx.star.MAF$x, df=14)-0.5)))
# Get q-values
qvalues.spring30.MAF = qvalue(pvalues.bilateral.spring30.MAF)
summary(qvalues.spring30.MAF)
# Create final matrix
spring30.original.xtx$MAF = spring30.MAF
spring30.original.xtx$star = spring30.xtx.star$x
spring30.original.xtx.MAF = spring30.original.xtx[spring30.original.xtx$MAF>=0.01,]
#write.table(spring30.original.xtx.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/spring_auto_MAF.txt", sep = "\t", row.names= F, quote=F)
spring30.original.xtx.MAF$qvalue = qvalues.spring30.MAF$qvalue
#write.table(spring30.original.xtx.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/without_inversions/spring_auto_MAF_qvalue.txt", sep = "\t", row.names= F, quote=F)

#significant.spring30.xtx.MAF = spring30.original.xtx.MAF[spring30.original.xtx.MAF$qvalue<=0.0001,]
# Write a table with all data:
#write.table(significant.spring30.xtx.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/spring30_auto_xtx_complete_table.txt", sep = "\t", row.names= F, quote=F)
# Write the bedfile:
#bedfile.spring30.MAF = significant.spring30.xtx.MAF[,c(8,9,9,6,12)]
#write.table(bedfile.spring30.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/spring30_auto_xtx.bed", sep = "\t", row.names= F, quote=F)

# New threshold (17-01-20):
num_sig_sp_auto = length(spring30.original.xtx.MAF$V1)*0.0005
significant.spring30.xtx.MAF.new = tail(spring30.original.xtx.MAF[order(spring30.original.xtx.MAF$star),], 700L)
#significant.spring30.xtx.MAF.new = spring30.original.xtx.MAF[spring30.original.xtx.MAF$V6>=27,]
#significant.spring30.xtx.MAF.new = spring30.original.xtx.MAF[spring30.original.xtx.MAF$qvalue<=0.000001,]
# Write a table with all data:
write.table(significant.spring30.xtx.MAF.new, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/spring30/spring30_auto_xtx_complete_table_top0.05_def.txt", sep = "\t", row.names= F, quote=F)
# Write the bedfile:
#bedfile.spring30.MAF.new = significant.spring30.xtx.MAF.new[,c(8,9,9,6,11,12)]
#write.table(bedfile.spring30.MAF.new, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/spring30_auto_xtx_top0.05.bed", sep = "\t", row.names= F, quote=F)
# For TEs with q-value <= 0-05:
# For TEs with top 0.1%:
# Finally we are keeping the 5% 70000 lines.
# Finally we are keeping the 2.5% 35000 lines.
significant.TE.spring30.xtx.MAF = tail(spring30.original.xtx.MAF[order(spring30.original.xtx.MAF$star),], 35000L)
write.table(significant.TE.spring30.xtx.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_sept_2020_std_corrected_variables/spring30/xtx/TEs/spring30_auto_xtx_TEs_top2.5%.txt", sep = "\t", row.names= F, quote=F)
significant.TE.spring30.xtx.MAF = tail(spring30.original.xtx.MAF[order(spring30.original.xtx.MAF$star),], 13998L)
write.table(significant.TE.spring30.xtx.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_sept_2020_std_corrected_variables/spring30/xtx/TEs/spring30_auto_xtx_TEs_top1%.txt", sep = "\t", row.names= F, quote=F)



####### EUROPE SPRING X CHROM ####### 
spring30.original.xtx.X = read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/spring30/XtX/star.xtx/spring30Xsoze_summary_pi_xtx.out", h=T)
spring30.xtx.X.names = read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/spring/30samples_batch/spring30_droseu14_X_names_SNPs.txt")
spring30.original.xtx.X$CHR = spring30.xtx.X.names$V1
spring30.original.xtx.X$POS = spring30.xtx.X.names$V2
spring30.star.xtx.X = read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/spring30/XtX/star.xtx/spring.star.xtx.X.txt", h=T)
#plot(spring30.original.xtx.X$M_XtX, spring30.star.xtx.X$x)
# Get MAF
MAF.spring.X = (0.5 - abs(0.5 -spring30.original.xtx.X$M_P))
spring30.star.xtx.X$MAF = MAF.spring.X
spring30.star.xtx.X.MAF = spring30.star.xtx.X[spring30.star.xtx.X$MAF>=0.01,]
# Get p-value
pvalues.bilateral.spring.X=1-(2*(abs(pchisq(spring30.star.xtx.X$x, df=14)-0.5)))
pvalues.bilateral.spring.X.MAF=1-(2*(abs(pchisq(spring30.star.xtx.X.MAF$x, df=14)-0.5)))
# Get q-value
qvalues.spring.X.MAF = qvalue(pvalues.bilateral.spring.X.MAF)
summary(qvalues.spring.X.MAF)
# Create final matrix)
spring30.original.xtx.X$MAF = MAF.spring.X
spring30.original.xtx.X$star = spring30.star.xtx.X$x
spring30.original.xtx.X.MAF = spring30.original.xtx.X[spring30.original.xtx.X$MAF>=0.01,]
#write.table(spring30.original.xtx.X.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/spring_Xchr_MAF.txt", sep = "\t", row.names= F, quote=F)
spring30.original.xtx.X.MAF$qvalues = qvalues.spring.X.MAF$qvalues
#write.table(spring30.original.xtx.X.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/without_inversions/spring_Xchr_MAF_qvalue.txt", sep = "\t", row.names= F, quote=F)

#significant.spring30.xtx.X = spring30.original.xtx.X.MAF[spring30.original.xtx.X.MAF$qvalues<=0.05,]
# Threshold: 18.210040
#def.spring30.X.MAF=significant.spring30.xtx.X[significant.spring30.xtx.X$M_XtX>=18.21,]
# Write a table with all data:
#write.table(def.spring30.X.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/spring30_Xchrom_xtx_complete_table.txt", sep = "\t", row.names= F, quote=F)
# Write the bedfile:
#bedfile.spring30.X.MAF = def.spring30.X.MAF[,c(8,9,9,6,12)]
#write.table(bedfile.spring30.X.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/spring30_Xchrom_xtx.bed", sep = "\t", row.names= F, quote=F)

# New threshold (17-02-20):
num_sig_sp_X = length(spring30.original.xtx.X.MAF$MRK)*0.0005
significant.spring30.xtx.X.new = tail(spring30.original.xtx.X.MAF[order(spring30.original.xtx.X.MAF$star),], 52L)
#significant.spring30.xtx.X.new = spring30.original.xtx.X.MAF[spring30.original.xtx.X.MAF$qvalues<=0.001,]
#significant.spring30.xtx.X.new=significant.spring30.xtx.X[significant.spring30.xtx.X$M_XtX>=23,]
# Write a table with all data:
write.table(significant.spring30.xtx.X.new, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/spring30/spring30_Xchrom_xtx_complete_table_top0.05_def.txt", sep = "\t", row.names= F, quote=F)
# Write the bedfile:
#bedfile.spring30.X.MAF.new = significant.spring30.xtx.X.new[,c(8,9,9,6,11,12)]
#write.table(bedfile.spring30.X.MAF.new, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/spring30_Xchrom_xtx_top0.05.bed", sep = "\t", row.names= F, quote=F)
# For TEs with threshold q-value <= 0.05:
# For TEs with threshold of top 0.1%:
# Finally we are keeping the 5% 5200 lines.
# Finally we are keeping the 2.5% 2600 lines.
significant.TE.spring30.xtx.X = tail(spring30.original.xtx.X.MAF[order(spring30.original.xtx.X.MAF$star),], 2613L)
write.table(significant.TE.spring30.xtx.X, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_sept_2020_std_corrected_variables/spring30/xtx/TEs/spring30_Xchrom_TE_top2.5%.txt", sep = "\t", row.names= F, quote=F)
significant.TE.spring30.xtx.X = tail(spring30.original.xtx.X.MAF[order(spring30.original.xtx.X.MAF$star),], 1045L)
write.table(significant.TE.spring30.xtx.X, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_sept_2020_std_corrected_variables/spring30/xtx/TEs/spring30_Xchrom_TE_top1%.txt", sep = "\t", row.names= F, quote=F)




####### EUROPE FALL AUTOSOMES #######
fall30.original.xtx=read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/fall30/XtX/star.xtx/fall30_all_concatenated_summary_pi_xtx.txt", h=F)
fall30.names=read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/fall/30samples_batch/splitted_autosomes/fall30_droseu14_auto.map.NAMES.ALL.txt")
fall30.original.xtx$CHROM = fall30.names$V1
fall30.original.xtx$POS = fall30.names$V2
fall30.xtx.star=read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/fall30/XtX/star.xtx/fall.star.xtx.txt",h=T)
#plot(fall30.original.xtx$V6, fall30.xtx.star$x)
# Get MAF
MAF.fall = (0.5 - abs(0.5 -fall30.original.xtx$V2))
fall30.xtx.star$MAF = MAF.fall
fall30.xtx.star.MAF = fall30.xtx.star[fall30.xtx.star$MAF>=0.01,]
# Get p-value
pvalues.bilateral.fall=1-(2*(abs(pchisq(fall30.xtx.star$x, df=10)-0.5)))
pvalues.bilateral.fall.MAF=1-(2*(abs(pchisq(fall30.xtx.star.MAF$x, df=10)-0.5)))
# Get q-value
qvalues.fall30.MAF = qvalue(pvalues.bilateral.fall.MAF)
summary(qvalues.fall30.MAF)
# Final matrix
fall30.original.xtx$MAF = MAF.fall
fall30.original.xtx$star = fall30.xtx.star$x
fall30.original.xtx.MAF = fall30.original.xtx[fall30.original.xtx$MAF>=0.01,]
#write.table(fall30.original.xtx.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/fall_auto_MAF.txt", sep = "\t", row.names= F, quote=F)
fall30.original.xtx.MAF$qvalues = qvalues.fall30.MAF$qvalues
#write.table(fall30.original.xtx.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/fall_auto_MAF_qvalue.txt", sep = "\t", row.names= F, quote=F)

#significant.fall30.xtx.MAF = fall30.original.xtx.MAF[fall30.original.xtx.MAF$qvalues<=0.001,]
# Write a table with all data:
#rite.table(significant.fall30.xtx.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/fall30_auto_xtx_complete_table.txt", sep = "\t", row.names= F, quote=F)
# Write the bedfile:
#bedfile.fall30.MAF = significant.fall30.xtx.MAF[,c(8,9,9,6,12)]
#write.table(bedfile.fall30.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/fall30_auto_xtx.bed", sep = "\t", row.names= F, quote=F)

# New threshold (17-01-20):
num_sig_fl_auto = length(fall30.original.xtx.MAF$V1)*0.0005
significant.fall30.xtx.MAF.new = tail(fall30.original.xtx.MAF[order(fall30.original.xtx.MAF$star),], 768L)
#significant.fall30.xtx.MAF.new = fall30.original.xtx.MAF[fall30.original.xtx.MAF$qvalues<=0.000005,]
#significant.fall30.xtx.MAF.new = fall30.original.xtx.MAF[fall30.original.xtx.MAF$V6>=20,]
# Write a table with all data:
write.table(significant.fall30.xtx.MAF.new, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/fall30_auto_xtx_complete_table_top0.05_prueba.txt", sep = "\t", row.names= F, quote=F)
# Write the bedfile:
#bedfile.fall30.MAF.new = significant.fall30.xtx.MAF.new[,c(8,9,9,6,11,12)]
#write.table(bedfile.fall30.MAF.new, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/fall30_auto_xtx_new_top0.05.bed", sep = "\t", row.names= F, quote=F)
# For TEs with q-value <= 0.05:
# For TEs with of top 0.1%:
# Finally we are keeping the 5% 76800 lines.
# Finally we are keeping the 2.5% 38400 lines.
significant.TE.fall30.xtx.MAF = tail(fall30.original.xtx.MAF[order(fall30.original.xtx.MAF$star),], 38400L)
write.table(significant.TE.fall30.xtx.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_sept_2020_std_corrected_variables/fall30/results/xtx/TEs/fall30_auto_TE_top2.5%.txt", sep = "\t", row.names= F, quote=F)
significant.TE.fall30.xtx.MAF = tail(fall30.original.xtx.MAF[order(fall30.original.xtx.MAF$star),], 15356L)
write.table(significant.TE.fall30.xtx.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_sept_2020_std_corrected_variables/fall30/results/xtx/TEs/fall30_auto_TE_top1%.txt", sep = "\t", row.names= F, quote=F)



####### EUROPE FALL X CHROM ####### 
fall30.original.xtx.X=read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/fall30/XtX/star.xtx/fall30Xsize_summary_pi_xtx.out", h=T)
fall30.X.names=read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/fall/30samples_batch/fall30_droseu14_X_names_SNPs.txt")
fall30.original.xtx.X$CHR = fall30.X.names$V1
fall30.original.xtx.X$POS = fall30.X.names$V2
fall30.xtx.star.X=read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/fall30/XtX/star.xtx/fall30.xtx.star.X.txt",h=T)
#plot(fall30.original.xtx.X$M_XtX, fall30.xtx.star.X$x)
# Get MAF
MAF.fall.X = (0.5 - abs(0.5 -fall30.original.xtx.X$M_P))
fall30.xtx.star.X$MAF = MAF.fall.X
fall30.xtx.star.X.MAF = fall30.xtx.star.X[fall30.xtx.star.X$MAF>=0.01,]
# Get p-value
pvalues.bilateral.fall.X=1-(2*(abs(pchisq(fall30.xtx.star.X$x, df=10)-0.5)))
pvalues.bilateral.fall.X.MAF=1-(2*(abs(pchisq(fall30.xtx.star.X.MAF$x, df=10)-0.5)))

# Get q-value
qvalues.fall.X.MAF = qvalue(pvalues.bilateral.fall.X.MAF)
summary(qvalues.fall.X.MAF)
# Final matrix
fall30.original.xtx.X$MAF = MAF.fall.X
fall30.original.xtx.X$star = fall30.xtx.star.X$x
fall30.original.xtx.X.MAF = fall30.original.xtx.X[fall30.original.xtx.X$MAF>=0.01,]
#write.table(fall30.original.xtx.X.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/fall_xchr_MAF.txt", sep = "\t", row.names= F, quote=F)
fall30.original.xtx.X.MAF$qvalues = qvalues.fall.X.MAF$qvalues
#write.table(fall30.original.xtx.X.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/fall_xchr_MAF_qvalue.txt", sep = "\t", row.names= F, quote=F)

#significant.fall30.xtx.X = fall30.original.xtx.X.MAF[fall30.original.xtx.X.MAF$qvalues<=0.05,]
# Threshold: 13.908057
#def.fall30.X.MAF=significant.fall30.xtx.X[significant.fall30.xtx.X$M_XtX>=13.9,]
# Write a table with all data:
#write.table(def.fall30.X.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/fall30_Xchrom_xtx_complete_table.txt", sep = "\t", row.names= F, quote=F)
# Write the bedfile:
#bedfile.fall30.X.MAF = def.fall30.X.MAF[,c(8,9,9,6,12)]
#write.table(bedfile.fall30.X.MAF, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/fall30_Xchrom_xtx.bed", sep = "\t", row.names= F, quote=F)

# New threshold (17-01-20):
num_sig_fl_X = length(fall30.original.xtx.X.MAF$MRK)*0.0005
significant.fall30.xtx.X.new = tail(fall30.original.xtx.X.MAF[order(fall30.original.xtx.X.MAF$star),], 53L)
#significant.fall30.xtx.X.new = fall30.original.xtx.X.MAF[fall30.original.xtx.X.MAF$qvalues<=0.001,]
#significant.fall30.xtx.X.new=significant.fall30.xtx.X[significant.fall30.xtx.X$M_XtX>=18,]
# Write a table with all data:
#write.table(significant.fall30.xtx.X.new, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/fall30_Xchrom_xtx_complete_table_top0.05_prueba.txt", sep = "\t", row.names= F, quote=F)
# Write the bedfile:
#bedfile.fall30.X.MAF.new = significant.fall30.xtx.X.new[,c(8,9,9,6,11,12)]
#write.table(bedfile.fall30.X.MAF.new, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/fall30_Xchrom_xtx_top0.05.bed", sep = "\t", row.names= F, quote=F)
# TE threshold of q-value <= 0.05:
# TE threshold of of top 0.1%:
# Finally we are keeping the 5% 5300 lines.
# Finally we are keeping the 5% 2650 lines.
significant.TE.fall30.xtx.X = tail(fall30.original.xtx.X.MAF[order(fall30.original.xtx.X.MAF$star),], 2650L)
write.table(significant.TE.fall30.xtx.X, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_sept_2020_std_corrected_variables/fall30/results/xtx/TEs/fall30_Xchrom_TE_top2.5%.txt", sep = "\t", row.names= F, quote=F)
significant.TE.fall30.xtx.X = tail(fall30.original.xtx.X.MAF[order(fall30.original.xtx.X.MAF$star),], 1060L)
write.table(significant.TE.fall30.xtx.X, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_sept_2020_std_corrected_variables/fall30/results/xtx/TEs/fall30_Xchrom_TE_top1%.txt", sep = "\t", row.names= F, quote=F)


