## Created 10th July 2019
## Simulations for the standard IS output
source("/Users/pogo/Documents/Maria/Baypass/baypass_2.1/utils/baypass_utils.R")
library(corrplot)
library(ape)
library(mvtnorm)
library(data.table)

#### Repeated on 19th September 2019

#### BERGLAND ####
### 2008
# Autosomes
# For the several simulations modified on 28th November 2019
omega=as.matrix(read.table("/Users/pogo/Documents/Maria/Baypass/Bergland/analysis/2008/auto/observed/omegas/bergland08ISAUTOstd_mat_omega.out"))
geno.data <- geno2YN("/Users/pogo/Documents/Maria/Baypass/Bergland/definitive_files_2019/bergland_auto_nofix.geno")
pi.beta.coef=read.table("/Users/pogo/Documents/Maria/Baypass/Bergland/analysis/2008/auto/observed/omegas/bergland08ISAUTOstd_summary_beta_params.out",h=T)$Mean
coverage= read.table("/Users/pogo/Documents/Maria/Baypass/Bergland/definitive_files_2019/bergland_auto_coverage.txt")
samplesize <- read.table("/Users/pogo/Documents/Maria/Baypass/Bergland/analysis/bergland_corrected.poolsize", header = FALSE)
simu.geno <- simulate.baypass(omega.mat=omega, nsnp=100000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0, coverage=coverage, suffix="bergland08stdAUTOsim")
simu.geno1 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland08stdAUTOsim1")
simu.geno2 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland08stdAUTOsim2")
simu.geno3 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland08stdAUTOsim3")
simu.geno4 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland08stdAUTOsim4")
simu.geno5 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland08stdAUTOsim5")
simu.geno6 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland08stdAUTOsim6")
simu.geno7 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland08stdAUTOsim7")
simu.geno8 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland08stdAUTOsim8")
simu.geno9 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland08stdAUTOsim9")
simu.geno10 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland08stdAUTOsim10")



# Chr. X (modified 21st October 2019)
# For the several simulations modified on 28th November 2019
omega=as.matrix(read.table("/Users/pogo/Documents/Maria/Baypass/Bergland/analysis/2008/X/for_sim/bergland08ISXstdsize_mat_omega.out"))
geno.data <- geno2YN("/Users/pogo/Documents/Maria/Baypass/Bergland/definitive_files_2019/bergland_X_nofix.geno")
pi.beta.coef=read.table("/Users/pogo/Documents/Maria/Baypass/Bergland/analysis/2008/X/for_sim/bergland08ISXstdsize_summary_beta_params.out",h=T)$Mean
coverage= read.table("/Users/pogo/Documents/Maria/Baypass/Bergland/definitive_files_2019/bergland_X_coverage.txt")
samplesize <- read.table("/Users/pogo/Documents/Maria/Baypass/Bergland/analysis/2008/X/for_sim/bergland_corrected_X.poolsize", header = FALSE)
simu.geno <- simulate.baypass(omega.mat=omega, nsnp=100000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0, coverage=coverage, suffix="bergland08stdXsimsize")
simu.geno1 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland08stdXsimsize1")
simu.geno2 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland08stdXsimsize2")
simu.geno3 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland08stdXsimsize3")
simu.geno4 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland08stdXsimsize4")
simu.geno5 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland08stdXsimsize5")
simu.geno6 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland08stdXsimsize6")
simu.geno7 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland08stdXsimsize7")
simu.geno8 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland08stdXsimsize8")
simu.geno9 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland08stdXsimsize9")
simu.geno10 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland08stdXsimsize10")


### 2010
# Autosomes
# For the several simulations modified on 28th November 2019
omega=as.matrix(read.table("/Users/pogo/Documents/Maria/Baypass/Bergland/analysis/2010/auto/observed/omegas/bergland10IS1AUTOstd_mat_omega.out"))
geno.data <- geno2YN("/Users/pogo/Documents/Maria/Baypass/Bergland/definitive_files_2019/bergland_auto_nofix.geno")
pi.beta.coef=read.table("/Users/pogo/Documents/Maria/Baypass/Bergland/analysis/2010/auto/observed/omegas/bergland10IS1AUTOstd_summary_beta_params.out",h=T)$Mean
coverage= read.table("/Users/pogo/Documents/Maria/Baypass/Bergland/definitive_files_2019/bergland_auto_coverage.txt")
samplesize <- read.table("/Users/pogo/Documents/Maria/Baypass/Bergland/analysis/bergland_corrected.poolsize", header = FALSE)
simu.geno <- simulate.baypass(omega.mat=omega, nsnp=100000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0, coverage=coverage, suffix="bergland10stdAUTOsim")
simu.geno1 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland10stdAUTOsim1")
simu.geno2 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland10stdAUTOsim2")
simu.geno3 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland10stdAUTOsim3")
simu.geno4 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland10stdAUTOsim4")
simu.geno5 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland10stdAUTOsim5")
simu.geno6 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland10stdAUTOsim6")
simu.geno7 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland10stdAUTOsim7")
simu.geno8 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland10stdAUTOsim8")
simu.geno9 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland10stdAUTOsim9")
simu.geno10 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland10stdAUTOsim10")


# 
# # Chr. X (modified 21st October 2019)
# For the several simulations modified on 28th November 2019
omega=as.matrix(read.table("/Users/pogo/Documents/Maria/Baypass/Bergland/analysis/2010/X/for_sim/bergland10IS1Xstdsize_mat_omega.out"))
geno.data <- geno2YN("/Users/pogo/Documents/Maria/Baypass/Bergland/definitive_files_2019/bergland_X_nofix.geno")
pi.beta.coef=read.table("/Users/pogo/Documents/Maria/Baypass/Bergland/analysis/2010/X/for_sim/bergland10IS1Xstdsize_summary_beta_params.out",h=T)$Mean
coverage= read.table("/Users/pogo/Documents/Maria/Baypass/Bergland/definitive_files_2019/bergland_X_coverage.txt")
samplesize <- read.table("/Users/pogo/Documents/Maria/Baypass/Bergland/analysis/2010/X/for_sim/bergland_corrected_X.poolsize", header = FALSE)
simu.geno <- simulate.baypass(omega.mat=omega, nsnp=100000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0, coverage=coverage, suffix="bergland10stdXsimsize")
simu.geno1 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland10stdXsimsize1")
simu.geno2 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland10stdXsimsize2")
simu.geno3 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland10stdXsimsize3")
simu.geno4 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland10stdXsimsize4")
simu.geno5 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland10stdXsimsize5")
simu.geno6 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland10stdXsimsize6")
simu.geno7 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland10stdXsimsize7")
simu.geno8 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland10stdXsimsize8")
simu.geno9 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland10stdXsimsize9")
simu.geno10 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="bergland10stdXsimsize10")


#### DrosEU 2014 ####

### FALL30
## Fall autosomes
# For the several simulations modified on 28th November 2019
omega=as.matrix(read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/fall30/standard/for_sim/fall30IS_25_mat_omega.out"))
geno.data <- geno2YN("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/fall/30samples_batch/fall30_droseu14_auto.geno")
pi.beta.coef=read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/fall30/standard/for_sim/fall30IS_25_summary_beta_params.out",h=T)$Mean
coverage= read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/fall/30samples_batch/fall30_droseu14_auto_coverage.txt")
samplesize <- read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/fall/30samples_batch/fall30.poolsize", header = FALSE)
simu.geno <- simulate.baypass(omega.mat=omega, nsnp=100000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0, coverage=coverage, suffix="fall30autoSTDsim")
simu.geno1 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="fall30autoSTDsim1")
simu.geno2 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="fall30autoSTDsim2")
simu.geno3 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="fall30autoSTDsim3")
simu.geno4 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="fall30autoSTDsim4")
simu.geno5 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="fall30autoSTDsim5")
simu.geno6 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="fall30autoSTDsim6")
simu.geno7 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="fall30autoSTDsim7")
simu.geno8 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="fall30autoSTDsim8")
simu.geno9 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="fall30autoSTDsim9")
simu.geno10 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="fall30autoSTDsim10")


## Fall X (modified on 21st October 2019)
# For the several simulations modified on 28th November 2019
omega=as.matrix(read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/fall30/standard/for_sim/fall30Xstdsize_mat_omega.out"))
geno.data <- geno2YN("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/fall/30samples_batch/fall30_droseu14_X.geno")
pi.beta.coef=read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/fall30/standard/for_sim/fall30Xstdsize_summary_beta_params.out",h=T)$Mean
coverage= read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/fall/30samples_batch/fall30_droseu14_X_coverage.txt")
samplesize <- read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/fall30/standard/for_sim/fall30_X.poolsize", header = FALSE)
simu.geno <- simulate.baypass(omega.mat=omega, nsnp=100000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0, coverage=coverage, suffix="fall30XSTDsimsize")
simu.geno1 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="fall30XSTDsimsize1")
simu.geno2 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="fall30XSTDsimsize2")
simu.geno3 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="fall30XSTDsimsize3")
simu.geno4 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="fall30XSTDsimsize4")
simu.geno5 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="fall30XSTDsimsize5")
simu.geno6 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="fall30XSTDsimsize6")
simu.geno7 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="fall30XSTDsimsize7")
simu.geno8 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="fall30XSTDsimsize8")
simu.geno9 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="fall30XSTDsimsize9")
simu.geno10 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="fall30XSTDsimsize10")



### SPRING30
# For the several simulations modified on 28th November 2019
omega=as.matrix(read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/spring30/standard/for_sim/spring30IS_25_mat_omega.out"))
geno.data <- geno2YN("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/spring/30samples_batch/spring30_droseu14_auto.geno")
pi.beta.coef=read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/spring30/standard/for_sim/spring30IS_25_summary_beta_params.out",h=T)$Mean
coverage= read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/spring/30samples_batch/spring30_droseu14_auto_coverage.txt")
samplesize <- read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/spring/30samples_batch/spring30.poolsize", header = FALSE)
simu.geno <- simulate.baypass(omega.mat=omega, nsnp=100000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0, coverage=coverage, suffix="spring30autoSTDsim")
simu.geno1 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="spring30autoSTDsim1")
simu.geno2 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="spring30autoSTDsim2")
simu.geno3 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="spring30autoSTDsim3")
simu.geno4 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="spring30autoSTDsim4")
simu.geno5 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="spring30autoSTDsim5")
simu.geno6 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="spring30autoSTDsim6")
simu.geno7 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="spring30autoSTDsim7")
simu.geno8 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="spring30autoSTDsim8")
simu.geno9 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="spring30autoSTDsim9")
simu.geno10 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="spring30autoSTDsim10")



## Spring X (modified on 21st October 2019)
# For the several simulations modified on 28th November 2019
omega=as.matrix(read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/spring30/standard/for_sim/spring30Xstdsize_mat_omega.out"))
geno.data <- geno2YN("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/spring/30samples_batch/spring30_droseu14_X.geno")
pi.beta.coef=read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/spring30/standard/for_sim/spring30Xstdsize_summary_beta_params.out",h=T)$Mean
coverage= read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/spring/30samples_batch/spring30_droseu14_X_coverage.txt")
samplesize <- read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/spring30/standard/for_sim/spring30_X.poolsize", header = FALSE)
simu.geno <- simulate.baypass(omega.mat=omega, nsnp=100000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0, coverage=coverage, suffix="spring30XSTDsimsize")
simu.geno1 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="spring30XSTDsimsize1")
simu.geno2 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="spring30XSTDsimsize2")
simu.geno3 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="spring30XSTDsimsize3")
simu.geno4 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="spring30XSTDsimsize4")
simu.geno5 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="spring30XSTDsimsize5")
simu.geno6 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="spring30XSTDsimsize6")
simu.geno7 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="spring30XSTDsimsize7")
simu.geno8 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="spring30XSTDsimsize8")
simu.geno9 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="spring30XSTDsimsize9")
simu.geno10 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="spring30XSTDsimsize10")



### ALL30
# AUTOSOMES (with splitted files)
# (Modified 28th November 2019)
omega=as.matrix(read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/all30/new/standard/auto/observed/omegas/all30IS_25_mat_omega.out"))
geno.data <- geno2YN("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/all30/corrected_TEs/all30_droseu14_noinv_auto.geno")
pi.beta.coef=read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/all30/new/standard/for_sim/observed/all30IS_25_summary_beta_params.out",h=T)$Mean
# In: /Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/all30/corrected_TEs
# python /Users/pogo/Documents/Maria/Scripts/Baypass/create_coverage_file_multiple.py all30_droseu14_noinv_auto.geno all30_droseu14_noinv_auto_coverage.txt 40
coverage= read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/all30/corrected_TEs/all30_droseu14_noinv_auto_coverage.txt")
samplesize <- read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/all30/corrected_TEs/all_30pop.poolsize", header = FALSE)
simu.geno <- simulate.baypass(omega.mat=omega, nsnp=100000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0, coverage=coverage, suffix="all30autoSTDsim")
simu.geno1 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="all30autoSTDsim1")
simu.geno2 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="all30autoSTDsim2")
simu.geno3 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="all30autoSTDsim3")
simu.geno4 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="all30autoSTDsim4")
simu.geno5 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="all30autoSTDsim5")
simu.geno6 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="all30autoSTDsim6")
simu.geno7 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="all30autoSTDsim7")
simu.geno8 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="all30autoSTDsim8")
simu.geno9 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="all30autoSTDsim9")
simu.geno10 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="all30autoSTDsim10")


# Chr. X (mofidied on 21st October 2019)
# (Multiple simulations modified on 28th November 2019)
omega=as.matrix(read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/all30/new/standard/X/for_sim/all30XnoinvSTDsize_mat_omega.out"))
geno.data <- geno2YN("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/all30/corrected_TEs/all30_droseu14_noinv_X.geno")
pi.beta.coef=read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/all30/new/standard/X/for_sim/all30XnoinvSTDsize_summary_beta_params.out",h=T)$Mean
# In: /Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/all30/corrected_TEs
# python /Users/pogo/Documents/Maria/Scripts/Baypass/create_coverage_file_multiple.py all30_droseu14_noinv_X.geno all30_droseu14_noinv_X_coverage.txt 40
coverage= read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/all30/corrected_TEs/all30_droseu14_noinv_X_coverage.txt")
samplesize <- read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/all30/new/standard/X/for_sim/all_30pop_X.poolsize", header = FALSE)
simu.geno <- simulate.baypass(omega.mat=omega, nsnp=100000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0, coverage=coverage, suffix="all30XNOINVsimsize")
simu.geno1 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="all30XNOINVsimsize1")
simu.geno2 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="all30XNOINVsimsize2")
simu.geno3 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="all30XNOINVsimsize3")
simu.geno4 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="all30XNOINVsimsize4")
simu.geno5 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="all30XNOINVsimsize5")
simu.geno6 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="all30XNOINVsimsize6")
simu.geno7 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="all30XNOINVsimsize7")
simu.geno8 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="all30XNOINVsimsize8")
simu.geno9 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="all30XNOINVsimsize9")
simu.geno10 <- simulate.baypass(omega.mat=omega, nsnp=50000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0.01, coverage=coverage, suffix="all30XNOINVsimsize10")





########## THIS WAY IS NOT CORRECT DUE TO WE ARE SIMULATING BY USING THE REAL DATA AND WE WANT NEUTRAL SIMULATIONS ############
## This loop is needed for each seed and each batch
## In bergland two batches (2010 and 2008) and four seeds (but we will use most probably only 3)
# pi.beta.coefis=read.table("/Users/pogo/Documents/Maria/Baypass/results_v2/Bergland/standard_IS/auto/2008/bergland08ISAUTO_summary_beta_params.out",h=T)$Mean
# berglandAUTO.data<-geno2YN("/Users/pogo/Documents/Maria/Baypass/results_v2/Bergland/files/bergland_autosomes.geno")
# beta.coefis=fread("/Users/pogo/Documents/Maria/Baypass/results_v2/Bergland/standard_IS/auto/2008/bergland08ISAUTO_summary_betai_reg.out",h=T)
# omega.is=as.matrix(read.table(file="/Users/pogo/Documents/Maria/Baypass/results_v2/Bergland/standard_IS/auto/2008/bergland08ISAUTO_mat_omega.out", header=F))
# cov=read.table("/Users/pogo/Documents/Maria/Baypass/results_v2/Bergland/files/variables/bergland_2008_standard_model.txt",h=F) # Environmental variables file
# samplesize=read.table("/Users/pogo/Documents/Maria/Baypass/results_v2/Bergland/files/bergland.poolsize")
# coveragefile=read.table("/Users/pogo/Documents/Maria/Baypass/results_v2/Bergland/files/bergland_auto_coverage.txt")
# for (i in 1:71) {
#   cov1=as.numeric(cov[i,])
#   beta.coefisvar=beta.coefis[beta.coefis$COVARIABLE == i, ]$Beta_is
#   simulate.baypass(omega.mat=omega.is, nsnp=100000, coverage=coveragefile, beta.coef=beta.coefisvar, beta.pi=pi.beta.coefis, pop.trait=cov1, sample.size=samplesize, pi.maf=0, suffix=paste("berglandAUTO_",i, sep=""))
# }

