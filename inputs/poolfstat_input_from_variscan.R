library(poolfstat)


# Machado 2L: "/Users/pogo/Documents/Maria/Baypass/Machado_DrosEU/mel_2L_mapping2016.varscanSNP_noindel_dpfilter_biallelic_repeatmask.recode.vcf"

poolsizes_machado <- read.table("/Users/pogo/Documents/Maria/Baypass/Machado_DrosEU/machado_number_individuals.txt")
#### Repasar los pool sizes porque deben ser iguales que en baypass (doble numero de individuos en autosomas y numero de individuos en X si machos)
poolnames_machado <- read.table("/Users/pogo/Documents/Maria/Baypass/Machado_DrosEU/machado_names.txt")

machado_info <- read.table("/Users/pogo/Documents/Maria/Baypass/Machado_DrosEU/final_order_poolsize.txt")

pooldata_machado2L <- vcf2pooldata(vcf.file = "/Users/pogo/Documents/Maria/Baypass/Machado_DrosEU/mel_2L_mapping2016.varscanSNP_noindel_dpfilter_biallelic_repeatmask.recode.vcf", 
                                   poolsizes = machado_info[,3], poolnames = machado_info[,1],
                                   min.cov.per.pool = -1, min.rc = 1, max.cov.per.pool = 1e+06,
                                   min.maf = 0.01, nlines.per.readblock = 1e+06, nthreads = 1)

pooldata_machado2L_subset <- pooldata.subset(pooldata_machado2L, pool.index = c(21, 22, 23, 25, 26, 27, 28, 34, 39, 50, 64), min.cov.per.pool = -1,
                                             max.cov.per.pool = 1e+06, min.maf = -1)
#write.table(pooldata_machado2L_subset@snp.info, file = "/Users/pogo/Documents/Maria/Baypass/NA/baypass/machado2L.subset.names")

pooldata_machado2L_subset_baypass <- pooldata2genobaypass(pooldata_machado2L_subset, writing.dir = getwd(), prefix = "machado2L.subset.pruebaRAL",
                     subsamplesize = -1, subsamplingmethod = "thinning")

pooldata_machado2R <- vcf2pooldata(vcf.file ="/Users/pogo/Documents/Maria/Baypass/NA/mel_2R_mapping2016.varscanSNP_noindel_dpfilter_biallelic_repeatmask.recode.vcf.gz",
                                  poolsizes = machado_info[,3], poolnames = machado_info[,1],
                                  min.cov.per.pool = -1, min.rc = 1, max.cov.per.pool = 1e+06,
                                  min.maf = 0.01, nlines.per.readblock = 1e+06, nthreads = 1)

pooldata_machado2R_subset <- pooldata.subset(pooldata_machado2R, pool.index = c(21, 22, 23, 25, 26, 27, 28, 34, 39, 50, 64), min.cov.per.pool = -1,
                                             max.cov.per.pool = 1e+06, min.maf = -1)

pooldata_machado3L <- vcf2pooldata(vcf.file ="/Users/pogo/Documents/Maria/Baypass/NA/mel_3L_mapping2016.varscanSNP_noindel_dpfilter_biallelic_repeatmask.recode.vcf.gz",
                                   poolsizes = machado_info[,3], poolnames = machado_info[,1],
                                   min.cov.per.pool = -1, min.rc = 1, max.cov.per.pool = 1e+06,
                                   min.maf = 0.01, nlines.per.readblock = 1e+06, nthreads = 1)

pooldata_machado3L_subset <- pooldata.subset(pooldata_machado3L, pool.index = c(21, 22, 23, 25, 26, 27, 28, 34, 39, 50, 64), min.cov.per.pool = -1,
                                             max.cov.per.pool = 1e+06, min.maf = -1)

pooldata_machado3R <- vcf2pooldata(vcf.file ="/Users/pogo/Documents/Maria/Baypass/NA/mel_3R_mapping2016.varscanSNP_noindel_dpfilter_biallelic_repeatmask.recode.vcf.gz",
                                   poolsizes = machado_info[,3], poolnames = machado_info[,1],
                                   min.cov.per.pool = -1, min.rc = 1, max.cov.per.pool = 1e+06,
                                   min.maf = 0.01, nlines.per.readblock = 1e+06, nthreads = 1)

pooldata_machado3R_subset <- pooldata.subset(pooldata_machado3R, pool.index = c(21, 22, 23, 25, 26, 27, 28, 34, 39, 50, 64), min.cov.per.pool = -1,
                                             max.cov.per.pool = 1e+06, min.maf = -1)

pooldata_machadoX <- vcf2pooldata(vcf.file ="/Users/pogo/Documents/Maria/Baypass/NA/mel_X_mapping2016.varscanSNP_noindel_dpfilter_biallelic_repeatmask.recode.vcf.gz",
                                   poolsizes = machado_info[,2], poolnames = machado_info[,1],
                                   min.cov.per.pool = -1, min.rc = 1, max.cov.per.pool = 1e+06,
                                   min.maf = 0.01, nlines.per.readblock = 1e+06, nthreads = 1)

pooldata_machadoX_subset <- pooldata.subset(pooldata_machadoX, pool.index = c(21, 22, 23, 25, 26, 27, 28, 34, 39, 50, 64), min.cov.per.pool = -1,
                                             max.cov.per.pool = 1e+06, min.maf = -1)

#genobaypass_machado2L <- pooldata2genobaypass(pooldata_machado2L, writing.dir = getwd(), prefix = "",
                     #subsamplesize = -1, subsamplingmethod = "thinning")


pooldata_machado2L_subset_baypass <- pooldata2genobaypass(pooldata_machado2L_subset, writing.dir = getwd(), prefix = "machado2L.subset",
                                                          subsamplesize = -1, subsamplingmethod = "thinning")

pooldata_machado2R_subset_baypass <- pooldata2genobaypass(pooldata_machado2R_subset, writing.dir = getwd(), prefix = "machado2R.subset",
                                                          subsamplesize = -1, subsamplingmethod = "thinning")

pooldata_machado3L_subset_baypass <- pooldata2genobaypass(pooldata_machado3L_subset, writing.dir = getwd(), prefix = "machado3L.subset",
                                                          subsamplesize = -1, subsamplingmethod = "thinning")

pooldata_machado3R_subset_baypass <- pooldata2genobaypass(pooldata_machado3R_subset, writing.dir = getwd(), prefix = "machado3R.subset",
                                                          subsamplesize = -1, subsamplingmethod = "thinning")

pooldata_machadoX_subset_baypass <- pooldata2genobaypass(pooldata_machadoX_subset, writing.dir = getwd(), prefix = "machadoX.subset.pruebaRAL",
                                                          subsamplesize = -1, subsamplingmethod = "thinning")


#### DEF SED COMMAND 20th AUGUST 2020
## 2L
pooldata_machado2L_NA <- vcf2pooldata(vcf.file = "/Users/pogo/Documents/Maria/Baypass/NA/mel_2L_mapping2016.varscanSNP_noindel_dpfilter_biallelic_repeatmask.recode.MODIFIED.vcf", 
                                   poolsizes = machado_info[,3], poolnames = machado_info[,1],
                                   min.cov.per.pool = -1, min.rc = 1, max.cov.per.pool = 1e+06,
                                   min.maf = 0.01, nlines.per.readblock = 1e+06, nthreads = 1)

pooldata_machado2L_subset_NA <- pooldata.subset(pooldata_machado2L_NA, pool.index = c(21, 22, 23, 25, 26, 27, 28, 34, 39, 50, 64), min.cov.per.pool = -1,
                                             max.cov.per.pool = 1e+06, min.maf = -1)
#write.table(pooldata_machado2L_subset@snp.info, file = "/Users/pogo/Documents/Maria/Baypass/NA/baypass/machado2L.subset.names")

pooldata_machado2L_subset_baypass <- pooldata2genobaypass(pooldata_machado2L_subset_NA, writing.dir = getwd(), prefix = "machado2L.subset.pruebaNA",
                                                          subsamplesize = -1, subsamplingmethod = "thinning")
## 2R
pooldata_machado2R_NA <- vcf2pooldata(vcf.file = "/Users/pogo/Documents/Maria/Baypass/NA/mel_2R_mapping2016.varscanSNP_noindel_dpfilter_biallelic_repeatmask.recode.MODIFIED.vcf", 
                                      poolsizes = machado_info[,3], poolnames = machado_info[,1],
                                      min.cov.per.pool = -1, min.rc = 1, max.cov.per.pool = 1e+06,
                                      min.maf = 0.01, nlines.per.readblock = 1e+06, nthreads = 1)

pooldata_machado2R_subset_NA <- pooldata.subset(pooldata_machado2R_NA, pool.index = c(21, 22, 23, 25, 26, 27, 28, 34, 39, 50, 64), min.cov.per.pool = -1,
                                                max.cov.per.pool = 1e+06, min.maf = -1)
#write.table(pooldata_machado2L_subset@snp.info, file = "/Users/pogo/Documents/Maria/Baypass/NA/baypass/machado2L.subset.names")

pooldata_machado2R_subset_baypass <- pooldata2genobaypass(pooldata_machado2R_subset_NA, writing.dir = getwd(), prefix = "machado2R.subset.pruebaNA",
                                                          subsamplesize = -1, subsamplingmethod = "thinning")


## 3L
pooldata_machado3L_NA <- vcf2pooldata(vcf.file = "/Users/pogo/Documents/Maria/Baypass/NA/mel_3L_mapping2016.varscanSNP_noindel_dpfilter_biallelic_repeatmask.recode.MODIFIED.vcf", 
                                      poolsizes = machado_info[,3], poolnames = machado_info[,1],
                                      min.cov.per.pool = -1, min.rc = 1, max.cov.per.pool = 1e+06,
                                      min.maf = 0.01, nlines.per.readblock = 1e+06, nthreads = 1)

pooldata_machado3L_subset_NA <- pooldata.subset(pooldata_machado3L_NA, pool.index = c(21, 22, 23, 25, 26, 27, 28, 34, 39, 50, 64), min.cov.per.pool = -1,
                                                max.cov.per.pool = 1e+06, min.maf = -1)
#write.table(pooldata_machado2L_subset@snp.info, file = "/Users/pogo/Documents/Maria/Baypass/NA/baypass/machado2L.subset.names")

pooldata_machado3L_subset_baypass <- pooldata2genobaypass(pooldata_machado3L_subset_NA, writing.dir = getwd(), prefix = "machado3L.subset.pruebaNA",
                                                          subsamplesize = -1, subsamplingmethod = "thinning")


## 3R
pooldata_machado3R_NA <- vcf2pooldata(vcf.file = "/Users/pogo/Documents/Maria/Baypass/NA/mel_3R_mapping2016.varscanSNP_noindel_dpfilter_biallelic_repeatmask.recode.MODIFIED.vcf", 
                                      poolsizes = machado_info[,3], poolnames = machado_info[,1],
                                      min.cov.per.pool = -1, min.rc = 1, max.cov.per.pool = 1e+06,
                                      min.maf = 0.01, nlines.per.readblock = 1e+06, nthreads = 1)

pooldata_machado3R_subset_NA <- pooldata.subset(pooldata_machado3R_NA, pool.index = c(21, 22, 23, 25, 26, 27, 28, 34, 39, 50, 64), min.cov.per.pool = -1,
                                                max.cov.per.pool = 1e+06, min.maf = -1)
#write.table(pooldata_machado2L_subset@snp.info, file = "/Users/pogo/Documents/Maria/Baypass/NA/baypass/machado2L.subset.names")

pooldata_machado3R_subset_baypass <- pooldata2genobaypass(pooldata_machado3R_subset_NA, writing.dir = getwd(), prefix = "machado3R.subset.pruebaNA",
                                                          subsamplesize = -1, subsamplingmethod = "thinning")


## X
pooldata_machadoX_NA <- vcf2pooldata(vcf.file = "/Users/pogo/Documents/Maria/Baypass/NA/mel_X_mapping2016.varscanSNP_noindel_dpfilter_biallelic_repeatmask.recode.MODIFIED.vcf", 
                                      poolsizes = machado_info[,2], poolnames = machado_info[,1],
                                      min.cov.per.pool = -1, min.rc = 1, max.cov.per.pool = 1e+06,
                                      min.maf = 0.01, nlines.per.readblock = 1e+06, nthreads = 1)

pooldata_machadoX_subset_NA <- pooldata.subset(pooldata_machadoX_NA, pool.index = c(21, 22, 23, 25, 26, 27, 28, 34, 39, 50, 64), min.cov.per.pool = -1,
                                                max.cov.per.pool = 1e+06, min.maf = -1)
#write.table(pooldata_machado2L_subset@snp.info, file = "/Users/pogo/Documents/Maria/Baypass/NA/baypass/machado2L.subset.names")

pooldata_machadoX_subset_baypass <- pooldata2genobaypass(pooldata_machadoX_subset_NA, writing.dir = getwd(), prefix = "machadoX.subset.pruebaNA",
                                                          subsamplesize = -1, subsamplingmethod = "thinning")

