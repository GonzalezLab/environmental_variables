## Created 3rd September 2020

# Analysis for 11 populations in Machado 2019

set.seed(61)

install.packages("DataCombine")
library(DataCombine)

### Loop for variables in autosomes:

for (i in 1:71) {
  print(i)
  stat = read.table(paste(paste("/Users/pogo/Documents/Maria/Baypass/NA/analysis/standard/auto/",i,sep=""),"_machadoAUTOstd_observed_complete.out", sep=""))
  stat.pi = read.table("/Users/pogo/Documents/Maria/Baypass/NA/analysis/standard/auto/all_machadoAUTOstd_summary_pi_xtx_names_sorted.out", h=F)
  stat0 = read.table(paste(paste("/Users/pogo/Documents/Maria/Baypass/NA/analysis/standard/simulated/auto/",i,sep=""),"_NAautoSIM_median_simulated.out", sep=""))
  stat0.pi = read.table("/Users/pogo/Documents/Maria/Baypass/NA/analysis/standard/simulated/auto/NAautoSIM_summary_pi_xtx.out",h=F)
  # Filter by MAF in observed data
  machadoNA.MAF = (0.5 - abs(0.5 -stat.pi$V5))
  stat$MAF = machadoNA.MAF
  stat.MAF = stat[stat$MAF>=0.01,]
  significant.30BF.machadoNA.MAF = stat.MAF[stat.MAF$V2>=30,]
  # Filter by MAF in simulated data
  machadoNA.MAF.sim = (0.5 - abs(0.5 -stat0.pi$V2))
  stat0$MAF = machadoNA.MAF.sim
  stat0.MAF = stat0[stat0$MAF>=0.01,]
  total.sim = length(stat0.MAF$V1)
  significant.30BF.machadoNA.MAF.sim = stat0.MAF[stat0.MAF$V2>=30,]
  FDR.calculation=(length(significant.30BF.machadoNA.MAF.sim$V1)*100)/total.sim
  pod.thresh=quantile(stat0.MAF$V2, probs=0.999)
  write.table(significant.30BF.machadoNA.MAF, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/NA/analysis/standard/results/auto/",i,sep=""),"_machadoNAauto_MAF_BF30.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  write.table(FDR.calculation, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/NA/analysis/standard/results/auto/",i,sep=""),"_machadoNAauto_MAF_FDRsim.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  write.table(pod.thresh, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/NA/analysis/standard/results/auto/",i,sep=""),"_machadoNAauto_threshold_sim_9.99%.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
}

### Loop for variables in X chromsoome:
for (i in 1:71) {
  print(i)
  stat = read.table(paste(paste("/Users/pogo/Documents/Maria/Baypass/NA/analysis/standard/X/",i,sep=""),"_machadoXstd_median_simulated_names.out", sep=""))
  stat.pi = read.table("/Users/pogo/Documents/Maria/Baypass/NA/analysis/standard/X/machadoXstd_summary_pi_xtx.out", h=F)
  stat0 = read.table(paste(paste("/Users/pogo/Documents/Maria/Baypass/NA/analysis/standard/simulated/X/",i,sep=""),"_NAsimX_median_simulated.out", sep=""))
  stat0.pi = read.table("/Users/pogo/Documents/Maria/Baypass/NA/analysis/standard/simulated/X/machadoXstdSIM_summary_pi_xtx.out",h=F)
  # Filter by MAF in observed data
  machadoNA.X.MAF = (0.5 - abs(0.5 -stat.pi$V2))
  stat$MAF = machadoNA.X.MAF
  stat.MAF = stat[stat$MAF>=0.01,]
  significant.30BF.machadoNA.X.MAF = stat.MAF[stat.MAF$V2>=30,]
  # Filter by MAF in simulated data
  all30.MAF.sim = (0.5 - abs(0.5 -stat0.pi$V2))
  stat0$MAF = all30.MAF.sim
  stat0.MAF = stat0[stat0$MAF>=0.01,]
  total.sim = length(stat0.MAF$V1)
  significant.30BF.all30.MAF.sim = stat0.MAF[stat0.MAF$V2>=30,]
  FDR.calculation=(length(significant.30BF.all30.MAF.sim$V1)*100)/total.sim
  pod.thresh=quantile(stat0.MAF$V2, probs=0.999)
  write.table(significant.30BF.machadoNA.X.MAF, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/NA/analysis/standard/results/X/",i,sep=""),"_machadoNA_X_MAF_BF30.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  write.table(FDR.calculation, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/NA/analysis/standard/results/X/",i,sep=""),"_machadoNA_X_MAF_FDRsim_BF30.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  write.table(pod.thresh, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/NA/analysis/standard/results/X/",i,sep=""),"_machadoNA_X_threshold_sim_9.99%.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
}


## Lighthours
# Autosomes
for (i in 1:7) {
  print(i)
  stat = read.table(paste(paste("/Users/pogo/Documents/Maria/Baypass/NA/analysis/standard/auto/lighthours/",i,sep=""),"_machadoAUTOstdLH_observed_complete.out", sep=""))
  stat.pi = read.table("/Users/pogo/Documents/Maria/Baypass/NA/analysis/standard/auto/all_machadoAUTOstd_summary_pi_xtx_names_sorted.out", h=F)
  stat0 = read.table(paste(paste("/Users/pogo/Documents/Maria/Baypass/NA/analysis/standard/simulated/auto/lighthours/",i,sep=""),"_NAautoSIMlh_median_simulated.out", sep=""))
  stat0.pi = read.table("/Users/pogo/Documents/Maria/Baypass/NA/analysis/standard/simulated/auto/lighthours/NAautoSIMlh_summary_pi_xtx.out",h=F)
  # Filter by MAF in observed data
  machadoNA.MAF = (0.5 - abs(0.5 -stat.pi$V5))
  stat$MAF = machadoNA.MAF
  stat.MAF = stat[stat$MAF>=0.01,]
  significant.30BF.machadoNA.MAF = stat.MAF[stat.MAF$V2>=30,]
  # Filter by MAF in simulated data
  machadoNA.MAF.sim = (0.5 - abs(0.5 -stat0.pi$V2))
  stat0$MAF = machadoNA.MAF.sim
  stat0.MAF = stat0[stat0$MAF>=0.01,]
  total.sim = length(stat0.MAF$V1)
  significant.30BF.machadoNA.MAF.sim = stat0.MAF[stat0.MAF$V2>=30,]
  FDR.calculation=(length(significant.30BF.machadoNA.MAF.sim$V1)*100)/total.sim
  pod.thresh=quantile(stat0.MAF$V2, probs=0.999)
  write.table(significant.30BF.machadoNA.MAF, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/NA/analysis/standard/results/auto/lighthours/",i,sep=""),"_machadoNAautoLH_MAF_BF30.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  write.table(FDR.calculation, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/NA/analysis/standard/results/auto/lighthours/",i,sep=""),"_machadoNAautoLH_MAF_FDRsim.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  write.table(pod.thresh, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/NA/analysis/standard/results/auto/lighthours/",i,sep=""),"_machadoNAautoLH_threshold_sim_9.99%.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
}

### Loop for variables in X chromsoome:
for (i in 1:7) {
  print(i)
  stat = read.table(paste(paste("/Users/pogo/Documents/Maria/Baypass/NA/analysis/standard/X/lighthours/",i,sep=""),"_machadoXstdLH_median_observed_names.out", sep=""))
  stat.pi = read.table("/Users/pogo/Documents/Maria/Baypass/NA/analysis/standard/X/lighthours/machadoXstdLH_summary_pi_xtx.out", h=F)
  stat0 = read.table(paste(paste("/Users/pogo/Documents/Maria/Baypass/NA/analysis/standard/simulated/X/lighthours/",i,sep=""),"_machadoXstdSIMlh_median_simulated.out", sep=""))
  stat0.pi = read.table("/Users/pogo/Documents/Maria/Baypass/NA/analysis/standard/simulated/X/lighthours/machadoXstdSIMlh_summary_pi_xtx.out",h=F)
  # Filter by MAF in observed data
  machadoNA.X.MAF = (0.5 - abs(0.5 -stat.pi$V2))
  stat$MAF = machadoNA.X.MAF
  stat.MAF = stat[stat$MAF>=0.01,]
  significant.30BF.machadoNA.X.MAF = stat.MAF[stat.MAF$V2>=30,]
  # Filter by MAF in simulated data
  all30.MAF.sim = (0.5 - abs(0.5 -stat0.pi$V2))
  stat0$MAF = all30.MAF.sim
  stat0.MAF = stat0[stat0$MAF>=0.01,]
  total.sim = length(stat0.MAF$V1)
  significant.30BF.all30.MAF.sim = stat0.MAF[stat0.MAF$V2>=30,]
  FDR.calculation=(length(significant.30BF.all30.MAF.sim$V1)*100)/total.sim
  pod.thresh=quantile(stat0.MAF$V2, probs=0.999)
  write.table(significant.30BF.machadoNA.X.MAF, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/NA/analysis/standard/results/X/lighthours/",i,sep=""),"_machadoNAlh_X_MAF_BF30.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  write.table(FDR.calculation, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/NA/analysis/standard/results/X/lighthours/",i,sep=""),"_machadoNAlh_MAF_FDRsim_BF30.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  write.table(pod.thresh, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/NA/analysis/standard/results/X/lighthours/",i,sep=""),"_machadoNAlh_threshold_sim_9.99%.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
}

