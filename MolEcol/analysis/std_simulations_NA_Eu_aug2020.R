## Created 28th August 2020
## Simulations for the standard IS output
## NA: for new dataset using 11 populations from Machado et al. 2019
## Europe: only for variables that are now corrected (mostly precipitation and evaporation)

source("/Users/pogo/Documents/Maria/Baypass/baypass_2.1/utils/baypass_utils.R")
library(corrplot)
library(ape)
library(mvtnorm)
library(data.table)


#### DrosEU 2014 ####

### FALL30
## Fall autosomes
omega=as.matrix(read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/corrected_var/fall30/data/auto/fall30ISvar_25_mat_omega.out"))
#geno.data <- geno2YN("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/fall/30samples_batch/fall30_droseu14_auto.geno")
pi.beta.coef=read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/corrected_var/fall30/data/auto/fall30ISvar_25_summary_beta_params.out",h=T)$Mean
coverage= read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/fall/30samples_batch/fall30_droseu14_auto_coverage.txt")
samplesize <- read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/fall/30samples_batch/fall30.poolsize", header = FALSE)
simu.geno <- simulate.baypass(omega.mat=omega, nsnp=100000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0, coverage=coverage, suffix="fall30autoSTDsimVAR")

## Fall X 
omega=as.matrix(read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/corrected_var/fall30/data/Xchrom/fall30Xstdvar_mat_omega.out"))
#geno.data <- geno2YN("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/fall/30samples_batch/fall30_droseu14_X.geno")
pi.beta.coef=read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/corrected_var/fall30/data/Xchrom/fall30Xstdvar_summary_beta_params.out",h=T)$Mean
coverage= read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/fall/30samples_batch/fall30_droseu14_X_coverage.txt")
samplesize <- read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/fall30/standard/for_sim/fall30_X.poolsize", header = FALSE)
simu.geno <- simulate.baypass(omega.mat=omega, nsnp=100000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0, coverage=coverage, suffix="fall30XSTDsimVAR")

### SPRING30
# AUTOSOMES
omega=as.matrix(read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/corrected_var/spring30/data/auto/spring30ISvar_25_mat_omega.out"))
#geno.data <- geno2YN("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/spring/30samples_batch/spring30_droseu14_auto.geno")
pi.beta.coef=read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/corrected_var/spring30/data/auto/spring30ISvar_25_summary_beta_params.out",h=T)$Mean
coverage= read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/spring/30samples_batch/spring30_droseu14_auto_coverage.txt")
samplesize <- read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/spring/30samples_batch/spring30.poolsize", header = FALSE)
simu.geno <- simulate.baypass(omega.mat=omega, nsnp=100000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0, coverage=coverage, suffix="spring30autoSTDsimVAR")

## Spring X 
omega=as.matrix(read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/corrected_var/spring30/data/Xchrom/spring30XstdVAR_mat_omega.out"))
#geno.data <- geno2YN("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/spring/30samples_batch/spring30_droseu14_X.geno")
pi.beta.coef=read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/corrected_var/spring30/data/Xchrom/spring30XstdVAR_summary_beta_params.out",h=T)$Mean
coverage= read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/spring/30samples_batch/spring30_droseu14_X_coverage.txt")
samplesize <- read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/spring30/standard/for_sim/spring30_X.poolsize", header = FALSE)
simu.geno <- simulate.baypass(omega.mat=omega, nsnp=100000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0, coverage=coverage, suffix="spring30XSTDsimVAR")


### ALL30
# AUTOSOMES (with splitted files)
omega=as.matrix(read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/corrected_var/all30/data/auto/all30ISvar_25_mat_omega.out"))
#geno.data <- geno2YN("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/all30/corrected_TEs/all30_droseu14_noinv_auto.geno")
pi.beta.coef=read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/corrected_var/all30/data/auto/all30ISvar_25_summary_beta_params.out",h=T)$Mean
# In: /Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/all30/corrected_TEs
# python /Users/pogo/Documents/Maria/Scripts/Baypass/create_coverage_file_multiple.py all30_droseu14_noinv_auto.geno all30_droseu14_noinv_auto_coverage.txt 40
coverage= read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/all30/corrected_TEs/all30_droseu14_noinv_auto_coverage.txt")
samplesize <- read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/all30/corrected_TEs/all_30pop.poolsize", header = FALSE)
simu.geno <- simulate.baypass(omega.mat=omega, nsnp=100000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0, coverage=coverage, suffix="all30autoSTDsimVAR")


# Chr. X 
omega=as.matrix(read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/corrected_var/all30/data/Xchrom/all30XnoinvSTDvar_mat_omega.out"))
#geno.data <- geno2YN("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/all30/corrected_TEs/all30_droseu14_noinv_X.geno")
pi.beta.coef=read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/corrected_var/all30/data/Xchrom/all30XnoinvSTDvar_summary_beta_params.out",h=T)$Mean
# In: /Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/all30/corrected_TEs
# python /Users/pogo/Documents/Maria/Scripts/Baypass/create_coverage_file_multiple.py all30_droseu14_noinv_X.geno all30_droseu14_noinv_X_coverage.txt 40
coverage= read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/all30/corrected_TEs/all30_droseu14_noinv_X_coverage.txt")
coverage2 = as.matrix(read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/all30/corrected_TEs/all30_droseu14_noinv_X_coverage.txt"))
samplesize <- read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/corrected_var/all30/data/all_30pop_X.poolsize", header = FALSE)
simu.geno <- simulate.baypass(omega.mat=omega, nsnp=100000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0, coverage=coverage, suffix="all30XNOINVsimVAR")


#### NA data (Machado 2019) ####
## AUTOSOMES
omega=as.matrix(read.table("/Users/pogo/Documents/Maria/Baypass/NA/simulations/data/machadoAUTOstd_25_mat_omega.out"))
#geno.data <- geno2YN("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/all30/corrected_TEs/all30_droseu14_noinv_auto.geno")
pi.beta.coef=read.table("/Users/pogo/Documents/Maria/Baypass/NA/simulations/data/machadoAUTOstd_25_summary_beta_params.out",h=T)$Mean
# In: /Users/pogo/Documents/Maria/Baypass/NA/baypass_data_def
#  python /Users/pogo/Documents/Maria/Scripts/Baypass/create_coverage_file_multiple_NA_2020.py machado.autosomes.TEs.sorted.geno machado.autosomes.TEs.sorted.coverage.txt 22
coverage= read.table("/Users/pogo/Documents/Maria/Baypass/NA/baypass_data_def/machado.autosomes.TEs.sorted.coverage.txt")
samplesize <- read.table("/Users/pogo/Documents/Maria/Baypass/NA/baypass_data_def/machado2L.subset.pruebaNA.poolsize", header = FALSE)
simu.geno <- simulate.baypass(omega.mat=omega, nsnp=100000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0, coverage=coverage, suffix="NAstdAUTOsim")

## CHROM X
omega=as.matrix(read.table("/Users/pogo/Documents/Maria/Baypass/NA/simulations/data/Xchrom/machadoXstd_mat_omega.out"))
pi.beta.coef=read.table("/Users/pogo/Documents/Maria/Baypass/NA/simulations/data/Xchrom/machadoXstd_summary_beta_params.out",h=T)$Mean
# /Users/pogo/Documents/Maria/Baypass/NA/baypass_data_def
#python /Users/pogo/Documents/Maria/Scripts/Baypass/create_coverage_file_multiple_NA_2020.py machadoX.subset.TEs.sorted.geno machadoX.subset.TEs.sorted.coverage.txt 22
coverage= read.table("/Users/pogo/Documents/Maria/Baypass/NA/baypass_data_def/machadoX.subset.TEs.sorted.coverage.txt")
samplesize <- read.table("/Users/pogo/Documents/Maria/Baypass/NA/simulations/data/Xchrom/machadoX.subset.pruebaNA.poolsize.good", header = FALSE)
simu.geno <- simulate.baypass(omega.mat=omega, nsnp=100000, sample.size=samplesize, beta.pi=pi.beta.coef,pi.maf=0, coverage=coverage, suffix="NAstdXsim")





