### 28th November 2019
## Created by Matthieu Gautier
## Adjust BFs and its PODs
## First filtering by MAF (this can be set as an argument)

obs.data <- read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/all30/new/standard/auto/observed/separated_seeds/22_all30IS_all_noseed_betai_names_head.txt", h=T)
obs.data.pi <- read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/all30/new/standard/auto/observed/separated_seeds/all30IS_all_noseed_pi_xtx.txt", h=F)

pods.data <- read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/all30/new/standard/auto/simulated/22_all30autoSTDsim_summary_betai_reg.out", h=F)
pods.data.pi <- read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/all30/new/standard/auto/simulated/all30autoSTDsim_summary_pi_xtx.out", h=T)

pods.data.500 <- read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/all30/new/XtX/calibration/var22_all30autoSTDsim_500K_summary_betai_reg.out", h=F)
pods.data.pi.500 <- read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/all30/new/XtX/calibration/all30autoSTDsim_500k_summary_pi_xtx.out", h=F)

#var22_pvalues_100 = compute.empPvals(obs.data$BF.dB.,obs.data.pi$V2,pods.data$V5,pods.data.pi$M_P)
var22_pvalues_500 = compute.empPvals(obs.data$BF.dB.,obs.data.pi$V2,pods.data.500$V5,pods.data.pi.500$V2)
var22_qvalues_500 = qvalue(var22_pvalues_500)
summary(var22_qvalues_500)


stat.obs=obs.data$BF.dB
pi.obs=obs.data.pi$V2
stat.pods=pods.data.500$V5
pi.pods=pods.data.pi.500$V2
maf.thr=0.05
tol=0.01
plot=TRUE


var22_pvalues_500_2 = compute.empPvals.2(obs.data$BF.dB.,obs.data.pi$V2,pods.data.500$V5,pods.data.pi.500$V2)

compute.empPvals<-function(stat.obs,pi.obs,stat.pods,pi.pods,maf.thr=0.05,tol=0.01,plot=TRUE){
  #tol=maximum difference in allele frequency bin porportion difference
  require(qvalue)
  maf.obs=0.5-abs(0.5-pi.obs)
  maf.pods=0.5-abs(0.5-pi.pods)
  sel.obs=maf.obs>maf.thr
  sel.pods=maf.pods>maf.thr
  breaks.vector=seq(max(0.005,maf.thr-0.005),min(0.995,(1-maf.thr)+0.005),0.01)
  
  freq.bins.obs=hist(pi.obs[sel.obs],breaks=breaks.vector,plot=F)
  prop.snp.per.bin=freq.bins.obs$counts/sum(freq.bins.obs$counts)
  names(prop.snp.per.bin)=freq.bins.obs$mids
  freq.bins.pods=hist(pi.pods[sel.pods],breaks=breaks.vector,plot=F)
  count.pods=freq.bins.pods$counts
  
  nsnps.to.sample=length(stat.pods)
  crit=1
  while(crit>tol & nsnps.to.sample>0){
    nsnps.to.sample=nsnps.to.sample-1
    counts.to.sample=round(prop.snp.per.bin*nsnps.to.sample)
    bins.with.not.enough.snps=prop.snp.per.bin*nsnps.to.sample>count.pods
    counts.to.sample[bins.with.not.enough.snps]=count.pods[bins.with.not.enough.snps]
    realized.prop=counts.to.sample/sum(counts.to.sample)
    crit=max(abs(realized.prop-prop.snp.per.bin))
    #  crit=max(prop.snp.per.bin*nsnps.to.sample-count.pods) #maximal nber of missing SNPs per category
  }
  
  if(nsnps.to.sample<1000){
    stop("Less than 1,000 eligible SNPs: You may increase the number of SNPs in the PODs\n")
  }
  
  sampled.snp=c()
  for(i in names(counts.to.sample)[counts.to.sample>0]){
    lim.inf=as.numeric(i)-0.005
    lim.sup=lim.inf+0.01
    tmp.num=which(pi.pods>=lim.inf & pi.pods<lim.sup)
    sampled.snp=c(sampled.snp,sample(tmp.num,counts.to.sample[i]))
  }
  
  res=matrix(NA,length(stat.obs),2)
  colnames(res)=c("Empirical Pvalue (adjusted PODS)","Empirical Pvalue (non-adjusted PODS)")
  res[sel.obs,1]=empPvals(stat = stat.obs[sel.obs],stat0 = stat.pods[sampled.snp])
  res[sel.obs,2]=empPvals(stat = stat.obs[sel.obs],stat0 = stat.pods)
  
  if(plot){
    layout(matrix(1:2,2,1))
    hist(res[sel.obs,1],breaks=50,freq=F,xlab="P-value",main="Empirical P-value distribution (adjusted PODs distribution)")
    abline(h=1,lty=2,col="red")
    hist(res[sel.obs,2],breaks=50,freq=F,xlab="P-value",main="Empirical P-value distribution (non-adjusted PODs distribution)")
    abline(h=1,lty=2,col="red")
  }
  
  return(res)
  
}


compute.empPvals.2 <-function(stat.obs,pi.obs,stat.pods,pi.pods,maf.thr=0.05,tol=0.01,plot=TRUE){
  #tol=maximum difference in allele frequency bin porportion difference
  require(qvalue)
  maf.obs=0.5-abs(0.5-pi.obs)
  maf.pods=0.5-abs(0.5-pi.pods)
  sel.obs=maf.obs>maf.thr
  sel.pods=maf.pods>maf.thr
  breaks.vector=seq(max(0.005,maf.thr-0.005),min(0.995,(1-maf.thr)+0.005),0.01)
  
  freq.bins.obs=hist(pi.obs[sel.obs],breaks=breaks.vector,plot=F)
  prop.snp.per.bin=freq.bins.obs$counts/sum(freq.bins.obs$counts)
  names(prop.snp.per.bin)=freq.bins.obs$mids
  freq.bins.pods=hist(pi.pods[sel.pods],breaks=breaks.vector,plot=F)
  count.pods=freq.bins.pods$counts
  
  nsnps.to.sample=length(stat.pods[sel.pods])
  crit=1
  while(crit>tol & nsnps.to.sample>0){
    nsnps.to.sample=nsnps.to.sample-1
    counts.to.sample=round(prop.snp.per.bin*nsnps.to.sample)
    bins.with.not.enough.snps=prop.snp.per.bin*nsnps.to.sample>count.pods
    counts.to.sample[bins.with.not.enough.snps]=count.pods[bins.with.not.enough.snps]
    realized.prop=counts.to.sample/sum(counts.to.sample)
    crit=max(abs(realized.prop-prop.snp.per.bin))
    #  crit=max(prop.snp.per.bin*nsnps.to.sample-count.pods) #maximal nber of missing SNPs per category
  }
  
  if(nsnps.to.sample<1000){
    stop("Less than 1,000 eligible SNPs: You may increase the number of SNPs in the PODs\n")
  }
  
  sampled.snp=c()
  for(i in names(counts.to.sample)[counts.to.sample>0]){
    lim.inf=as.numeric(i)-0.005
    lim.sup=lim.inf+0.01
    tmp.num=which(pi.pods>=lim.inf & pi.pods<lim.sup)
    sampled.snp=c(sampled.snp,sample(tmp.num,counts.to.sample[i]))
  }
  
  stat.pods.maf=stat.pods[sel.pods]
  
  
  res=matrix(NA,length(stat.obs),2)
  colnames(res)=c("Empirical Pvalue (adjusted PODS)","Empirical Pvalue (non-adjusted PODS)")
  res[sel.obs,1]=empPvals(stat = stat.obs[sel.obs],stat0 = stat.pods[sampled.snp])
  res[sel.obs,2]=empPvals(stat = stat.obs[sel.obs],stat0 = stat.pods.maf)
  
  if(plot){
    layout(matrix(1:2,2,1))
    hist(res[sel.obs,1],breaks=50,freq=F,xlab="P-value",main="Empirical P-value distribution (adjusted PODs distribution)")
    abline(h=1,lty=2,col="red")
    hist(res[sel.obs,2],breaks=50,freq=F,xlab="P-value",main="Empirical P-value distribution (non-adjusted PODs distribution)")
    abline(h=1,lty=2,col="red")
  }
  
  return(res)
  
}
compute.empPvals.3 <-function(stat.obs,pi.obs,stat.pods,pi.pods,maf.thr=0.1,tol=0.01,plot=TRUE){
  #tol=maximum difference in allele frequency bin porportion difference
  require(qvalue)
  maf.obs=0.5-abs(0.5-pi.obs)
  maf.pods=0.5-abs(0.5-pi.pods)
  sel.obs=maf.obs>maf.thr
  sel.pods=maf.pods>maf.thr
  breaks.vector=seq(max(0.005,maf.thr-0.005),min(0.995,(1-maf.thr)+0.005),0.01)
  
  freq.bins.obs=hist(pi.obs[sel.obs],breaks=breaks.vector,plot=F)
  prop.snp.per.bin=freq.bins.obs$counts/sum(freq.bins.obs$counts)
  names(prop.snp.per.bin)=freq.bins.obs$mids
  freq.bins.pods=hist(pi.pods[sel.pods],breaks=breaks.vector,plot=F)
  count.pods=freq.bins.pods$counts
  
  nsnps.to.sample=length(stat.pods[sel.pods])
  crit=1
  while(crit>tol & nsnps.to.sample>0){
    nsnps.to.sample=nsnps.to.sample-1
    counts.to.sample=round(prop.snp.per.bin*nsnps.to.sample)
    bins.with.not.enough.snps=prop.snp.per.bin*nsnps.to.sample>count.pods
    counts.to.sample[bins.with.not.enough.snps]=count.pods[bins.with.not.enough.snps]
    realized.prop=counts.to.sample/sum(counts.to.sample)
    crit=max(abs(realized.prop-prop.snp.per.bin))
    #  crit=max(prop.snp.per.bin*nsnps.to.sample-count.pods) #maximal nber of missing SNPs per category
  }
  
  if(nsnps.to.sample<1000){
    stop("Less than 1,000 eligible SNPs: You may increase the number of SNPs in the PODs\n")
  }
  
  sampled.snp=c()
  for(i in names(counts.to.sample)[counts.to.sample>0]){
    lim.inf=as.numeric(i)-0.005
    lim.sup=lim.inf+0.01
    tmp.num=which(pi.pods>=lim.inf & pi.pods<lim.sup)
    sampled.snp=c(sampled.snp,sample(tmp.num,counts.to.sample[i]))
  }
  
  stat.pods.maf=stat.pods[sel.pods]
  
  
  res=matrix(NA,length(stat.obs),2)
  colnames(res)=c("Empirical Pvalue (adjusted PODS)","Empirical Pvalue (non-adjusted PODS)")
  res[sel.obs,1]=empPvals(stat = stat.obs[sel.obs],stat0 = stat.pods[sampled.snp])
  res[sel.obs,2]=empPvals(stat = stat.obs[sel.obs],stat0 = stat.pods.maf)
  
  if(plot){
    layout(matrix(1:2,2,1))
    hist(res[sel.obs,1],breaks=50,freq=F,xlab="P-value",main="Empirical P-value distribution (adjusted PODs distribution)")
    abline(h=1,lty=2,col="red")
    hist(res[sel.obs,2],breaks=50,freq=F,xlab="P-value",main="Empirical P-value distribution (non-adjusted PODs distribution)")
    abline(h=1,lty=2,col="red")
  }
  
  return(res)
  
}
compute.empPvals.4 <-function(stat.obs,pi.obs,stat.pods,pi.pods,maf.thr=0.01,tol=0.01,plot=TRUE){
  #tol=maximum difference in allele frequency bin porportion difference
  require(qvalue)
  maf.obs=0.5-abs(0.5-pi.obs)
  maf.pods=0.5-abs(0.5-pi.pods)
  sel.obs=maf.obs>maf.thr
  sel.pods=maf.pods>maf.thr
  breaks.vector=seq(max(0.005,maf.thr-0.005),min(0.995,(1-maf.thr)+0.005),0.01)
  
  freq.bins.obs=hist(pi.obs[sel.obs],breaks=breaks.vector,plot=F)
  prop.snp.per.bin=freq.bins.obs$counts/sum(freq.bins.obs$counts)
  names(prop.snp.per.bin)=freq.bins.obs$mids
  freq.bins.pods=hist(pi.pods[sel.pods],breaks=breaks.vector,plot=F)
  count.pods=freq.bins.pods$counts
  
  nsnps.to.sample=length(stat.pods[sel.pods])
  crit=1
  while(crit>tol & nsnps.to.sample>0){
    nsnps.to.sample=nsnps.to.sample-1
    counts.to.sample=round(prop.snp.per.bin*nsnps.to.sample)
    bins.with.not.enough.snps=prop.snp.per.bin*nsnps.to.sample>count.pods
    counts.to.sample[bins.with.not.enough.snps]=count.pods[bins.with.not.enough.snps]
    realized.prop=counts.to.sample/sum(counts.to.sample)
    crit=max(abs(realized.prop-prop.snp.per.bin))
    #  crit=max(prop.snp.per.bin*nsnps.to.sample-count.pods) #maximal nber of missing SNPs per category
  }
  
  if(nsnps.to.sample<1000){
    stop("Less than 1,000 eligible SNPs: You may increase the number of SNPs in the PODs\n")
  }
  
  sampled.snp=c()
  for(i in names(counts.to.sample)[counts.to.sample>0]){
    lim.inf=as.numeric(i)-0.005
    lim.sup=lim.inf+0.01
    tmp.num=which(pi.pods>=lim.inf & pi.pods<lim.sup)
    sampled.snp=c(sampled.snp,sample(tmp.num,counts.to.sample[i]))
  }
  
  stat.pods.maf=stat.pods[sel.pods]
  
  
  res=matrix(NA,length(stat.obs),2)
  colnames(res)=c("Empirical Pvalue (adjusted PODS)","Empirical Pvalue (non-adjusted PODS)")
  res[sel.obs,1]=empPvals(stat = stat.obs[sel.obs],stat0 = stat.pods[sampled.snp])
  res[sel.obs,2]=empPvals(stat = stat.obs[sel.obs],stat0 = stat.pods.maf)
  
  if(plot){
    layout(matrix(1:2,2,1))
    hist(res[sel.obs,1],breaks=50,freq=F,xlab="P-value",main="Empirical P-value distribution (adjusted PODs distribution)")
    abline(h=1,lty=2,col="red")
    hist(res[sel.obs,2],breaks=50,freq=F,xlab="P-value",main="Empirical P-value distribution (non-adjusted PODs distribution)")
    abline(h=1,lty=2,col="red")
  }
  
  return(res)
  
}

set.seed(61)

install.packages("DataCombine")
library(DataCombine)

### Loop for variables in All30 (autosomes):
for (i in 1:69) {
  print(i)
  stat = read.table(paste(paste("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/all30/new/standard/auto/observed/",i,sep=""),"_all30IS_observed_complete.out", sep=""))
  stat.pi = read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/all30/new/standard/auto/observed/all_TE_all30IS_all_SORTED_summary_pi_xtx_names.out", h=F)
  stat0 = read.table(paste(paste("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/all30/new/standard/auto/simulated/",i,sep=""),"_all30autoSTDsim_median_simulated.out", sep=""))
  stat0.pi = read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/all30/new/standard/auto/simulated/all30autoSTDsim_summary_pi_xtx.out",h=F)
  # Filter by MAF in observed data
  all30.MAF = (0.5 - abs(0.5 -stat.pi$V4))
  stat$MAF = all30.MAF
  stat.MAF = stat[stat$MAF>=0.01,]
  significant.30BF.all30.MAF = stat.MAF[stat.MAF$V2>=30,]
  # Filter by MAF in simulated data
  all30.MAF.sim = (0.5 - abs(0.5 -stat0.pi$V2))
  stat0$MAF = all30.MAF.sim
  stat0.MAF = stat0[stat0$MAF>=0.01,]
  total.sim = length(stat0.MAF$V1)
  significant.30BF.all30.MAF.sim = stat0.MAF[stat0.MAF$V2>=30,]
  FDR.calculation=(length(significant.30BF.all30.MAF.sim$V1)*100)/total.sim
  pod.thresh=quantile(stat0.MAF$V2, probs=0.999)
  # write.table(significant.30BF.all30.MAF, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/all30/auto/",i,sep=""),"_all30ISauto_MAF_BF20.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  # #write.table(FDR.calculation, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/all30/auto/",i,sep=""),"_all30ISauto_MAF_FDRsim.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  #write.table(pod.thresh, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/all30/auto/",i,sep=""),"_all30ISauto_threshold_sim_9.99%.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  TE_eu_auto <- grepl.sub(data = stat.MAF, pattern = "FBti", Var = "V3")
  write.table(TE_eu_auto, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/all30/auto/",i,sep=""),"_all30ISauto_MAF_all_TEs.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
}

### Loop for variables in All30 (X chromosome):
for (i in 1:69) {
  print(i)
  stat = read.table(paste(paste("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/all30/new/standard/X/observed/",i,sep=""),"_all30XnoinvSTDsize_median_observed_names.out", sep=""))
  stat.pi = read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/all30/new/standard/X/observed/all30XnoinvSTDsize_summary_pi_xtx.out", h=F)
  stat0 = read.table(paste(paste("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/all30/new/standard/X/simulated/",i,sep=""),"_all30XnoinvSTDsimsize_median_simulated.out", sep=""))
  stat0.pi = read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/all30/new/standard/X/simulated/all30XnoinvSTDsimsize_summary_pi_xtx.out",h=F)
  # Filter by MAF in observed data
  all30.MAF = (0.5 - abs(0.5 -stat.pi$V2))
  stat$MAF = all30.MAF
  stat.MAF = stat[stat$MAF>=0.01,]
  significant.30BF.all30.MAF = stat.MAF[stat.MAF$V2>=30,]
  # Filter by MAF in simulated data
  all30.MAF.sim = (0.5 - abs(0.5 -stat0.pi$V2))
  stat0$MAF = all30.MAF.sim
  stat0.MAF = stat0[stat0$MAF>=0.01,]
  total.sim = length(stat0.MAF$V1)
  significant.30BF.all30.MAF.sim = stat0.MAF[stat0.MAF$V2>=30,]
  FDR.calculation=(length(significant.30BF.all30.MAF.sim$V1)*100)/total.sim
  pod.thresh=quantile(stat0.MAF$V2, probs=0.999)
  # write.table(significant.30BF.all30.MAF, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/all30/X/",i,sep=""),"_all30ISX_MAF_BF20.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  # write.table(FDR.calculation, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/all30/X/",i,sep=""),"_all30ISX_MAF_FDRsim_BF20.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  write.table(pod.thresh, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/all30/X/",i,sep=""),"_all30ISX_threshold_sim_9.99%.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  #TE_eu_X <- grepl.sub(data = stat.MAF, pattern = "FBti", Var = "V3")
  #write.table(TE_eu_X, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/all30/X/",i,sep=""),"_all30ISX_MAF_all_TEs.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
}


### Loop for variables in Bergland 2008 (autosomes):
for (i in 1:71) {
  print(i)
  stat = read.table(paste(paste("/Users/pogo/Documents/Maria/Baypass/Bergland/analysis/2008/auto/observed/",i,sep=""),"_bergland08ISAUTOstd_median_observed_names.out", sep=""))
  stat.pi = read.table("/Users/pogo/Documents/Maria/Baypass/Bergland/analysis/2008/auto/observed/bergland08ISAUTOstd_summary_pi_xtx.out", h=F)
  stat0 = read.table(paste(paste("/Users/pogo/Documents/Maria/Baypass/Bergland/analysis/2008/auto/simulated/",i,sep=""),"_bergland08ISAUTOsim_median_simulated.out", sep=""))
  stat0.pi = read.table("/Users/pogo/Documents/Maria/Baypass/Bergland/analysis/2008/auto/simulated/bergland08ISAUTOsim_summary_pi_xtx.out",h=F)
  # Filter by MAF in observed data
  all30.MAF = (0.5 - abs(0.5 -stat.pi$V2))
  stat$MAF = all30.MAF
  stat.MAF = stat[stat$MAF>=0.01,]
  significant.30BF.all30.MAF = stat.MAF[stat.MAF$V2>=30,]
  # Filter by MAF in simulated data
  # all30.MAF.sim = (0.5 - abs(0.5 -stat0.pi$V2))
  # stat0$MAF = all30.MAF.sim
  # stat0.MAF = stat0[stat0$MAF>=0.01,]
  # total.sim = length(stat0.MAF$V1)
  # significant.30BF.all30.MAF.sim = stat0.MAF[stat0.MAF$V2>=30,]
  # FDR.calculation=(length(significant.30BF.all30.MAF.sim$V1)*100)/total.sim
  # pod.thresh=quantile(stat0.MAF$V2, probs=0.999)
  #write.table(significant.30BF.all30.MAF, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/NA08/auto/",i,sep=""),"_NA08_MAF_BF30.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  #write.table(FDR.calculation, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/NA08/auto/",i,sep=""),"_NA08_MAF_FDRsim_BF30.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  #write.table(pod.thresh, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/NA/NA08/auto/",i,sep=""),"_threshold_sim_9.99%.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  TE_NA08_auto <- grepl.sub(data = stat.MAF, pattern = "FBti", Var = "V3")
  write.table(TE_NA08_auto, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/NA/NA08/auto/",i,sep=""),"_NA08_MAF_all_TEs.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
}

### Loop for variables in Bergland 2008 (X chromosome):
for (i in 1:71) {
  print(i)
  stat = read.table(paste(paste("/Users/pogo/Documents/Maria/Baypass/Bergland/analysis/2008/X/observed/",i,sep=""),"_bergland08ISXstdsize_median_observed_names.out", sep=""))
  stat.pi = read.table("/Users/pogo/Documents/Maria/Baypass/Bergland/analysis/2008/X/observed/bergland08ISXstdsize_summary_pi_xtx.out", h=F)
  stat0 = read.table(paste(paste("/Users/pogo/Documents/Maria/Baypass/Bergland/analysis/2008/X/simulated/",i,sep=""),"_bergland08ISXsimsize_median_simulated.out", sep=""))
  stat0.pi = read.table("/Users/pogo/Documents/Maria/Baypass/Bergland/analysis/2008/X/simulated/bergland08ISXsimsize_summary_pi_xtx.out",h=F)
  # Filter by MAF in observed data
  all30.MAF = (0.5 - abs(0.5 -stat.pi$V2))
  stat$MAF = all30.MAF
  stat.MAF = stat[stat$MAF>=0.01,]
  significant.30BF.all30.MAF = stat.MAF[stat.MAF$V2>=30,]
  # Filter by MAF in simulated data
  all30.MAF.sim = (0.5 - abs(0.5 -stat0.pi$V2))
  stat0$MAF = all30.MAF.sim
  stat0.MAF = stat0[stat0$MAF>=0.01,]
  total.sim = length(stat0.MAF$V1)
  significant.30BF.all30.MAF.sim = stat0.MAF[stat0.MAF$V2>=30,]
  FDR.calculation=(length(significant.30BF.all30.MAF.sim$V1)*100)/total.sim
  pod.thresh=quantile(stat0.MAF$V2, probs=0.999)
  # write.table(significant.30BF.all30.MAF, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/NA08/X/",i,sep=""),"_NA08_X_MAF_BF30.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  # write.table(FDR.calculation, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/NA08/X/",i,sep=""),"_NA08_X_MAF_FDRsim_BF30.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  write.table(pod.thresh, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/NA/NA08/X/",i,sep=""),"_NA08_X_threshold_sim_9.99%.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  #TE_NA08_X <- grepl.sub(data = stat.MAF, pattern = "FBti", Var = "V3")
  #write.table(TE_NA08_X, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/NA/NA08/X/",i,sep=""),"_NA08_X_MAF_all_TEs.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
}

### Loop for variables in Bergland 2010 (autosomes):
for (i in 1:72) {
  print(i)
  stat = read.table(paste(paste("/Users/pogo/Documents/Maria/Baypass/Bergland/analysis/2010/auto/observed/",i,sep=""),"_bergland10IS1AUTOstd_median_observed_names.out", sep=""))
  stat.pi = read.table("/Users/pogo/Documents/Maria/Baypass/Bergland/analysis/2010/auto/observed/bergland10IS1AUTOstd_summary_pi_xtx.out", h=F)
  stat0 = read.table(paste(paste("/Users/pogo/Documents/Maria/Baypass/Bergland/analysis/2010/auto/simulated/",i,sep=""),"_bergland10ISAUTOsim_median_simulated.out", sep=""))
  stat0.pi = read.table("/Users/pogo/Documents/Maria/Baypass/Bergland/analysis/2010/auto/simulated/bergland10ISAUTOsim_summary_pi_xtx.out",h=F)
  # Filter by MAF in observed data
  all30.MAF = (0.5 - abs(0.5 -stat.pi$V2))
  stat$MAF = all30.MAF
  stat.MAF = stat[stat$MAF>=0.01,]
  significant.30BF.all30.MAF = stat.MAF[stat.MAF$V2>=30,]
  # Filter by MAF in simulated data
  # all30.MAF.sim = (0.5 - abs(0.5 -stat0.pi$V2))
  # stat0$MAF = all30.MAF.sim
  # stat0.MAF = stat0[stat0$MAF>=0.01,]
  # total.sim = length(stat0.MAF$V1)
  # significant.30BF.all30.MAF.sim = stat0.MAF[stat0.MAF$V2>=30,]
  # FDR.calculation=(length(significant.30BF.all30.MAF.sim$V1)*100)/total.sim
  # pod.thresh=quantile(stat0.MAF$V2, probs=0.999)
  #write.table(significant.30BF.all30.MAF, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/NA10/auto/",i,sep=""),"_NA10_MAF_BF30.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  #write.table(FDR.calculation, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/NA10/auto/",i,sep=""),"_NA10_MAF_FDRsim_BF30.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  #write.table(pod.thresh, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/NA/NA10/auto/",i,sep=""),"_NA10_threshold_sim_9.99%.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")  
  TE_NA10_auto <- grepl.sub(data = stat.MAF, pattern = "FBti", Var = "V3")
  write.table(TE_NA10_auto, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/NA/NA10/auto/",i,sep=""),"_NA10_MAF_all_TEs.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  
}

### Loop for variables in Bergland 2010 (X chrom):
for (i in 1:72) {
  print(i)
  stat = read.table(paste(paste("/Users/pogo/Documents/Maria/Baypass/Bergland/analysis/2010/X/observed/",i,sep=""),"_bergland10IS1Xstdsize_median_observed_names.out", sep=""))
  stat.pi = read.table("/Users/pogo/Documents/Maria/Baypass/Bergland/analysis/2010/X/observed/bergland10IS1Xstdsize_summary_pi_xtx.out", h=F)
  stat0 = read.table(paste(paste("/Users/pogo/Documents/Maria/Baypass/Bergland/analysis/2010/X/simulated/",i,sep=""),"_bergland10ISXsimsize_median_simulated.out", sep=""))
  stat0.pi = read.table("/Users/pogo/Documents/Maria/Baypass/Bergland/analysis/2010/X/simulated/bergland10ISXsimsize_summary_pi_xtx.out",h=F)
  # Filter by MAF in observed data
  all30.MAF = (0.5 - abs(0.5 -stat.pi$V2))
  stat$MAF = all30.MAF
  stat.MAF = stat[stat$MAF>=0.01,]
  significant.30BF.all30.MAF = stat.MAF[stat.MAF$V2>=30,]
  # Filter by MAF in simulated data
  all30.MAF.sim = (0.5 - abs(0.5 -stat0.pi$V2))
  stat0$MAF = all30.MAF.sim
  stat0.MAF = stat0[stat0$MAF>=0.01,]
  total.sim = length(stat0.MAF$V1)
  significant.30BF.all30.MAF.sim = stat0.MAF[stat0.MAF$V2>=30,]
  FDR.calculation=(length(significant.30BF.all30.MAF.sim$V1)*100)/total.sim
  pod.thresh=quantile(stat0.MAF$V2, probs=0.999)
  # write.table(significant.30BF.all30.MAF, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/NA10/X/",i,sep=""),"_NA10_X_MAF_BF30.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  # write.table(FDR.calculation, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/NA10/X/",i,sep=""),"_NA10_X_MAF_FDRsim_BF30.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  write.table(pod.thresh, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/NA/NA10/X/",i,sep=""),"_NA10_X_threshold_sim_9.99%.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  #TE_NA10_X <- grepl.sub(data = stat.MAF, pattern = "FBti", Var = "V3")
  #write.table(TE_NA10_X, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/NA/NA10/X/",i,sep=""),"_NA10_X_MAF_all_TEs.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
}

### Loop for variables in Spring30 (autosomes):
for (i in 1:68) {
  print(i)
  stat = read.table(paste(paste("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/spring30/standard/auto/observed/",i,sep=""),"_spring30IS_observed_complete.out", sep=""))
  stat.pi = read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/spring30/standard/auto/observed/all_TE_spring30IS_all_SORTED_summary_pi_xtx_names.out", h=F)
  stat0 = read.table(paste(paste("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/spring30/standard/auto/simulated/",i,sep=""),"_spring30autoSTDsim_median_simulated.out", sep=""))
  stat0.pi = read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/spring30/standard/auto/simulated/spring30autoSTDsim_summary_pi_xtx.out",h=F)
  # Filter by MAF in observed data
  all30.MAF = (0.5 - abs(0.5 -stat.pi$V4))
  stat$MAF = all30.MAF
  stat.MAF = stat[stat$MAF>=0.01,]
  significant.30BF.all30.MAF = stat.MAF[stat.MAF$V2>=30,]
  # Filter by MAF in simulated data
  # all30.MAF.sim = (0.5 - abs(0.5 -stat0.pi$V2))
  # stat0$MAF = all30.MAF.sim
  # stat0.MAF = stat0[stat0$MAF>=0.01,]
  # total.sim = length(stat0.MAF$V1)
  # significant.30BF.all30.MAF.sim = stat0.MAF[stat0.MAF$V2>=30,]
  # FDR.calculation=(length(significant.30BF.all30.MAF.sim$V1)*100)/total.sim
  # pod.thresh=quantile(stat0.MAF$V2, probs=0.999)
  #write.table(significant.30BF.all30.MAF, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/spring30/auto/",i,sep=""),"_spring30ISauto_MAF_BF30.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  #write.table(FDR.calculation, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/spring30/auto/",i,sep=""),"_spring30ISauto_MAF_FDRsim.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  #write.table(pod.thresh, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/spring30/auto/",i,sep=""),"_spring30ISauto_threshold_sim_9.99%.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  TE_SP_auto <- grepl.sub(data = stat.MAF, pattern = "FBti", Var = "V3")
  write.table(TE_SP_auto, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/spring30/auto/",i,sep=""),"_spring30ISauto_MAF_all_TEs.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
}

### Loop for variables in Spring30 (X chromosome):
for (i in 1:68) {
  print(i)
  stat = read.table(paste(paste("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/spring30/standard/X/observed/",i,sep=""),"_spring30Xstds567size_median_observed_names.out", sep=""))
  stat.pi = read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/spring30/standard/X/observed/spring30Xstdsize_summary_pi_xtx.out", h=F)
  stat0 = read.table(paste(paste("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/spring30/standard/X/simulated/",i,sep=""),"_spring30XstdSIMsize_median_simulated.out", sep=""))
  stat0.pi = read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/spring30/standard/X/simulated/spring30XstdSIMsize_summary_pi_xtx.out",h=F)
  # Filter by MAF in observed data
  all30.MAF = (0.5 - abs(0.5 -stat.pi$V2))
  stat$MAF = all30.MAF
  stat.MAF = stat[stat$MAF>=0.01,]
  significant.30BF.all30.MAF = stat.MAF[stat.MAF$V2>=30,]
  # Filter by MAF in simulated data
  all30.MAF.sim = (0.5 - abs(0.5 -stat0.pi$V2))
  stat0$MAF = all30.MAF.sim
  stat0.MAF = stat0[stat0$MAF>=0.01,]
  total.sim = length(stat0.MAF$V1)
  significant.30BF.all30.MAF.sim = stat0.MAF[stat0.MAF$V2>=30,]
  FDR.calculation=(length(significant.30BF.all30.MAF.sim$V1)*100)/total.sim
  pod.thresh=quantile(stat0.MAF$V2, probs=0.999)
  # write.table(significant.30BF.all30.MAF, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/spring30/X/",i,sep=""),"_spring30ISX_MAF_BF30.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  # write.table(FDR.calculation, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/spring30/X/",i,sep=""),"_spring30ISX_MAF_FDRsim_BF30.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  write.table(pod.thresh, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/spring30/X/",i,sep=""),"_spring30ISX_threshold_sim_9.99%.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  #TE_SP_X <- grepl.sub(data = stat.MAF, pattern = "FBti", Var = "V3")
  #write.table(TE_SP_X, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/spring30/X/",i,sep=""),"_spring30ISX_MAF_all_TEs.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  
}

### Loop for variables in Fall30 (autosomes):
for (i in 1:68) {
  print(i)
  stat = read.table(paste(paste("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/fall30/standard/auto/observed/",i,sep=""),"_fall30IS_observed_complete.out", sep=""))
  stat.pi = read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/fall30/standard/auto/observed/all_TE_fall30IS_all_SORTED_summary_pi_xtx_names.out", h=F)
  stat0 = read.table(paste(paste("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/fall30/standard/auto/simulated/",i,sep=""),"_fall30autoSTDsim_median_simulated.out", sep=""))
  stat0.pi = read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/fall30/standard/auto/simulated/fall30AUTOstdSIM_summary_pi_xtx.out",h=F)
  # Filter by MAF in observed data
  all30.MAF = (0.5 - abs(0.5 -stat.pi$V4))
  stat$MAF = all30.MAF
  stat.MAF = stat[stat$MAF>=0.01,]
  significant.30BF.all30.MAF = stat.MAF[stat.MAF$V2>=30,]
  # Filter by MAF in simulated data
  all30.MAF.sim = (0.5 - abs(0.5 -stat0.pi$V2))
  stat0$MAF = all30.MAF.sim
  stat0.MAF = stat0[stat0$MAF>=0.01,]
  total.sim = length(stat0.MAF$V1)
  significant.30BF.all30.MAF.sim = stat0.MAF[stat0.MAF$V2>=30,]
  FDR.calculation=(length(significant.30BF.all30.MAF.sim$V1)*100)/total.sim
  pod.thresh=quantile(stat0.MAF$V2, probs=0.999)
  # write.table(significant.30BF.all30.MAF, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/fall30/auto/",i,sep=""),"_fall30ISauto_MAF_BF30.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  # write.table(FDR.calculation, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/fall30/auto/",i,sep=""),"_fall30ISauto_MAF_FDRsim.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  write.table(pod.thresh, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/fall30/auto/",i,sep=""),"_fall30ISauto_threshold_sim_9.99.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  #TE_FL_auto <- grepl.sub(data = stat.MAF, pattern = "FBti", Var = "V3")
  #write.table(TE_FL_auto, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/fall30/auto/",i,sep=""),"_fall30ISauto_MAF_all_TEs.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
}

### Loop for variables in Fall30 (X chromosome):
for (i in 1:68) {
  print(i)
  stat = read.table(paste(paste("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/fall30/standard/X/observed/",i,sep=""),"_fall30Xstdsize_median_observed_names.out", sep=""))
  stat.pi = read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/fall30/standard/X/observed/fall30Xstdsize_summary_pi_xtx.out", h=F)
  stat0 = read.table(paste(paste("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/fall30/standard/X/simulated/",i,sep=""),"_fall30XstdSIMsize_median_simulated.out", sep=""))
  stat0.pi = read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/fall30/standard/X/simulated/fall30XstdSIMsize_summary_pi_xtx.out",h=F)
  # Filter by MAF in observed data
  all30.MAF = (0.5 - abs(0.5 -stat.pi$V2))
  stat$MAF = all30.MAF
  stat.MAF = stat[stat$MAF>=0.01,]
  significant.30BF.all30.MAF = stat.MAF[stat.MAF$V2>=30,]
  # Filter by MAF in simulated data
  all30.MAF.sim = (0.5 - abs(0.5 -stat0.pi$V2))
  stat0$MAF = all30.MAF.sim
  stat0.MAF = stat0[stat0$MAF>=0.01,]
  total.sim = length(stat0.MAF$V1)
  significant.30BF.all30.MAF.sim = stat0.MAF[stat0.MAF$V2>=30,]
  FDR.calculation=(length(significant.30BF.all30.MAF.sim$V1)*100)/total.sim
  pod.thresh=quantile(stat0.MAF$V2, probs=0.999)
  # write.table(significant.30BF.all30.MAF, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/fall30/X/",i,sep=""),"_fall30ISX_MAF_BF30.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  # write.table(FDR.calculation, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/fall30/X/",i,sep=""),"_fall30ISX_MAF_FDRsim_BF30.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  write.table(pod.thresh, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/fall30/X/",i,sep=""),"_fall30ISX_threshold_sim_9.99%",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  #TE_FL_X <- grepl.sub(data = stat.MAF, pattern = "FBti", Var = "V3")
  #write.table(TE_FL_X, file=paste(paste("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/fall30/X/",i,sep=""),"_fall30ISX_MAF_all_TEs.txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
}


############################################# PRUEBAS ############################################# 
### Test if graphs are similar using different seeds and median of BF:
obs.data.seed <- read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/all30/new/standard/auto/observed/separated_seeds/22_all30ISs1500_all_betai_names_head.out", h=T)
obs.data.seed.pi <- read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/all30/new/standard/auto/observed/separated_seeds/all30ISs1500_all_summary_pi_xtx.out", h=F)
#pods.data <- read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/all30/new/standard/auto/simulated/22_all30autoSTDsim_summary_betai_reg.out", h=F)
#pods.data.pi <- read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/all30/new/standard/auto/simulated/all30autoSTDsim_summary_pi_xtx.out", h=T)
pods.data.500 <- read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/all30/new/XtX/calibration/var22_all30autoSTDsim_500K_summary_betai_reg.out", h=F)
pods.data.pi.500 <- read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/all30/new/XtX/calibration/all30autoSTDsim_500k_summary_pi_xtx.out", h=F)

var22_pvalues_500_seed1 = compute.empPvals.2(obs.data.seed$BF.dB.,obs.data.seed.pi$V2,pods.data.500$V5,pods.data.pi.500$V2)

# Test different variables:
obs.data.var11 <- read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/all30/new/standard/auto/observed/separated_seeds/11_all30IS_all_betai_names_head.out", h=T)
obs.data.seed.pi <- read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/all30/new/standard/auto/observed/separated_seeds/all30ISs1500_all_summary_pi_xtx.out", h=F)
pods.data.500.var11 <- read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/all30/new/standard/auto/simulated_500k/11_all30autoSTDsim_500K_summary_betai_reg.out", h=F)
pods.data.pi.500 <- read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/all30/new/XtX/calibration/all30autoSTDsim_500k_summary_pi_xtx.out", h=F)
var22_pvalues_500_var11 = compute.empPvals.2(obs.data.var11$BF.dB.,obs.data.seed.pi$V2,pods.data.500.var11$V5,pods.data.pi.500$V2)

obs.data.var51 <- read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/all30/new/standard/auto/observed/separated_seeds/51_all30IS_all_betai_names_head.out", h=T)
obs.data.seed.pi <- read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/all30/new/standard/auto/observed/separated_seeds/all30ISs1500_all_summary_pi_xtx.out", h=F)
pods.data.500.var51 <- read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/all30/new/standard/auto/simulated_500k/51_all30autoSTDsim_500K_summary_betai_reg.out", h=F)
pods.data.pi.500 <- read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/all30/new/XtX/calibration/all30autoSTDsim_500k_summary_pi_xtx.out", h=F)
var51_pvalues_500 = compute.empPvals.2(obs.data.var51$BF.dB.,obs.data.seed.pi$V2,pods.data.500.var51$V5,pods.data.pi.500$V2)


# With median values (what pi value to use???):

