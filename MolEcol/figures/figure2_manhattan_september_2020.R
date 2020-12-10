### Created 31st January 2020

# Manhattan plot for figure 2:

library(qvalue)

## GGPLOT
# https://github.com/drveera/ggman
install.packages("devtools", "gtools", "XML", "digest", "rlang")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BiocGenerics")
BiocManager::install("S4Vectors")
BiocManager::install("Biostrings")
BiocManager::install("SummarizedExperiment")
BiocManager::install("GenomicAlignments")
library(devtools)
library(ggplot2)
library(gtools)
library(XML)
library(digest)
library(rlang)
library(BiocGenerics)
library(S4Vectors)
library(Rsamtools)
library(githubinstall)
install_github("drveera/ggman")
library(ggman)

# Transposable elements for plot:
te_dict <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/for_manhattan_plot/TEcopies_1630_dictionary_manhattan.txt")
#te_dict_v5 <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/NA/corrected_TE_order/NA_all_TEs_coordinates_for_plot_v5_ALL.txt")

####### EUROPE ####### 
#all30.xtx=read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/all30/new/XtX/xtx_star/all30_concatenated_summary_pi_xtx.out", h=F)
all30.xtx=read.table("/Volumes/Maxtor/Maria/Baypass/13-02-20/all30/new/XtX/xtx_star/all30_concatenated_summary_pi_xtx.out", h=F)
names=read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/all30/corrected_TEs/splitted_geno/map/all30_droseu14_auto.geno.names.perfiles.txt")
all30.xtx$V8=names$V1
all30.xtx$V9=names$V2
#xtx.star=read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/all30/new/XtX/xtx_star/xtx.star.txt",h=T)
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
all30.xtx.MAF$qvalue=qvalues.xtx.MAF$qvalues

# Chrom X
#all30.xtx.X = read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/all30/new/XtX/xtx_star/all30Xnoinvsize_summary_pi_xtx.out", h=T)
all30.xtx.X = read.table("/Volumes/Maxtor/Maria/Baypass/13-02-20/all30/new/XtX/xtx_star/all30Xnoinvsize_summary_pi_xtx.out", h=T)
all30.X.names=read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/all30/corrected_TEs/all30_droseu14_X_SNPs_names.txt", h=F)
all30.xtx.X$CHROM = all30.X.names$V1
all30.xtx.X$POS = all30.X.names$V2
#all30.star.xtx.X = read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/all30/new/XtX/xtx_star/xtx.star.X.txt", h=F)
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
all30.xtx.X.MAF$qvalue = qvalues.all30.X.MAF$qvalue

# Manhattan:
xtx_europe_all_significant <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/all30/all30_all_xtx_complete_prueba_short.txt")
# For removing MAF from the original dataset:
xtx_europe_MAF = all30.xtx.MAF[c(8,9,6,11)] # Obtained from the xtx_square.R script calculation
xtx_europe_all_MAF_X = all30.xtx.X.MAF[c(8,9,6,11)]
colnames(xtx_europe_MAF) <- colnames(xtx_europe_all_MAF_X)
xtx_europe_all_MAF <- rbind(xtx_europe_MAF, xtx_europe_all_MAF_X)
xtx_europe_allSNPs_MAF = transform(xtx_europe_all_MAF, names=paste(CHROM, POS, sep=""))
xtx_europe_all_MAF_names = xtx_europe_allSNPs_MAF$names
xtx_europe_all_MAF_chromosome = xtx_europe_allSNPs_MAF$CHROM

# Substitute TE names:
fbti_vector <- grep("FBti", xtx_europe_all_MAF_names)
te_vector <- vector()
for (i in fbti_vector) {
  te_num = grep(xtx_europe_all_MAF_names[i],te_dict$V1)
  te_vector <- c(te_vector, te_num)
}
te_subs <- te_dict$V2[te_vector]
chrom_subs <- te_dict$V3[te_vector]
xtx_europe_allSNPs_MAF_names_changed <- replace(as.character(xtx_europe_all_MAF_names), fbti_vector ,as.character(te_subs))
xtx_europe_allSNPs_MAF$names <- xtx_europe_allSNPs_MAF_names_changed
xtx_europe_allSNPs_MAF_chromosome_changed <- replace(as.character(xtx_europe_all_MAF_chromosome), fbti_vector ,as.character(chrom_subs))
xtx_europe_allSNPs_MAF$chromosome <- xtx_europe_allSNPs_MAF_chromosome_changed
# Continue with the rest of the plot
to_highlight <- as.vector(xtx_europe_all_significant[,1])
# This is the moment to change the MAF background or not
p1 <- ggman(xtx_europe_allSNPs_MAF, snp = "names", bp = "POS", chrom = "chromosome", pvalue = "star", logTransform = FALSE, ymax = 240, sigLine = NA, xlab = "", ylab ="XtX* values", title = "Europe")
highlighted <- ggmanHighlight(p1, highlight = to_highlight, colour = "steelblue1", size = 0.8)
inversions_all <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/all30/all30_ALL_inversions_v6.txt")
inversions_all_def <- as.vector(inversions_all[,1])
all30 <- ggmanHighlight(highlighted, highlight = inversions_all_def, colour = "steelblue4", size = 0.8)
TEs_all <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/for_manhattan_plot/all30_all_TEs_coordinates_for_plot.txt")
TEs_all_def <- as.vector(TEs_all[,1])
all30_all_TE <- ggmanHighlight(all30, highlight = TEs_all_def, colour = "black", size = 0.8)
TEs_all_sig <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_sept_2020_std_corrected_variables/all30/xtx/TEs/all30_significant_TEs_1percent.txt")
TEs_all_sig_def <- as.vector(TEs_all_sig[,1])
all30_TE <- ggmanHighlight(all30_all_TE, highlight = TEs_all_sig_def, colour = "red", size = 1)
label_file <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/all30/all30_labels_for_plot.txt")
all30_TE_labels <- ggmanLabel(all30_TE, labelDfm = label_file, snp = "V1", label = "V2")


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


## Xchrom ##
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

# Manhattan plot
xtx_NA_significant <- read.table("/Users/pogo/Documents/Maria/Baypass/NA/analysis/plots/machadoNA_all_significant_SNPs.txt")
## For removing MAF from the original dataset:
xtx_NA_auto <- xtx.pi.NA.MAF[c(8,9,6,12)]
xtx_NA_X <- xtx.pi.NA.Xchrom.MAF[c(8,9,6,12)]
colnames(xtx_NA_auto) <- colnames(xtx_NA_X)
xtx_NA_MAF <- rbind(xtx_NA_auto,xtx_NA_X)
xtx_NA_all_MAF = transform(xtx_NA_MAF, names=paste(CHROM, POS, sep=""))
xtx_NA_all_MAF_names = xtx_NA_all_MAF$names
xtx_NA_all_MAF_chromosome = xtx_NA_all_MAF$CHR
# # Substitute TE names:
# fbti_vector_NA <- grep("FBti", xtx_NA_all_MAF_names)
# te_vector_NA <- vector()
# for (i in fbti_vector_NA) {
#   te_num_NA = grep(xtx_NA_all_MAF_names[i],te_dict_v5$V1)
#   te_vector_NA <- c(te_vector_NA, te_num_NA)
# }
# te_subs_NA <- te_dict_v5$V2[te_vector_NA]
# chrom_subs_NA <- te_dict_v5$V3[te_vector_NA]
# xtx_NA_all_MAF_names_changed <- replace(as.character(xtx_NA_all_MAF_names), fbti_vector_NA ,as.character(te_subs_NA))
# xtx_NA_all_MAF$names <- xtx_NA_all_MAF_names_changed
# xtx_NA_all_MAF_chromosome_changed <- replace(as.character(xtx_NA_all_MAF_chromosome), fbti_vector_NA ,as.character(chrom_subs_NA))
# xtx_NA_all_MAF$chromosome <- xtx_NA_all_MAF_chromosome_changed

to_highlight_NA <- as.vector(xtx_NA_significant[,1])
p2 <- ggman(xtx_NA_all_MAF, snp = "names", bp = "POS", chrom = "CHROM", pvalue = "xtxstar", logTransform = FALSE, ymax = 100, sigLine = NA, xlab = "", ylab ="XtX* values", title = "North America")
highlighted_NA <- ggmanHighlight(p2, highlight = to_highlight_NA, colour = "steelblue1", size = 0.8)
inversions_NA <- read.table("/Users/pogo/Documents/Maria/Baypass/NA/analysis/core/xtx/inversions/NAmachado_all_inversions_v5.txt")
inversions_NA_def <- as.vector(inversions_NA[,1])
NorthAmerica <- ggmanHighlight(highlighted_NA, highlight = inversions_NA_def, colour = "steelblue4", size = 0.8)
TEs_NA <- read.table("/Users/pogo/Documents/Maria/Baypass/NA/analysis/plots/machadoNA_TE_coordinates.txt")
TEs_NA_def <- as.vector(TEs_NA[,1])
NA_all_TE <- ggmanHighlight(NorthAmerica, highlight = TEs_NA_def, colour = "black", size = 0.8)
NA_all_sig <- read.table("/Users/pogo/Documents/Maria/Baypass/NA/analysis/plots/machadoNA_significant_TEs_qvalue_filter.txt")
NA_all_sig_def <- as.vector(NA_all_sig[,1])
NA_TE <- ggmanHighlight(NA_all_TE, highlight = NA_all_sig_def, colour = "red", size = 1)
label_file_NA <- read.table("/Users/pogo/Documents/Maria/Baypass/NA/analysis/plots/NA_labels_for_plot.txt")
NA_TE_labels <- ggmanLabel(NA_TE, labelDfm = label_file_NA, snp = "V1", label = "V2")






####### SPRING ####### 
spring30.original.xtx=read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/spring30/XtX/star.xtx/spring30_all_concatenated_summary_pi_xtx.txt", h=F)
spring30.names=read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/spring/30samples_batch/splitted_autosomes/spring30_droseu14_auto.map.NAMES.ALL.txt")
spring30.original.xtx$CHROM = spring30.names$V1
spring30.original.xtx$POS = spring30.names$V2
spring30.xtx.star=read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/spring30/XtX/star.xtx/spring30.xtx.star.txt",h=T)
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
spring30.original.xtx.MAF$qvalue = qvalues.spring30.MAF$qvalue
# X CHromosome
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
spring30.original.xtx.X.MAF$qvalues = qvalues.spring.X.MAF$qvalues

# Manhattan plot:
xtx_europe_spring_significant <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/spring30/spring30_all_significant.txt")
to_highlight_spring <- as.vector(xtx_europe_spring_significant[,1])

## For removing MAF from the original dataset:
xtx_spring_auto <- spring30.original.xtx.MAF[c(8,9,6,11)]
xtx_spring_X <- spring30.original.xtx.X.MAF[c(8,9,6,11)]
colnames(xtx_spring_auto) <- colnames(xtx_spring_X)
xtx_spring_MAF <- rbind(xtx_spring_auto,xtx_spring_X)
xtx_spring_all_MAF = transform(xtx_spring_MAF, names=paste(CHR, POS, sep=""))
xtx_spring_all_MAF_names = xtx_spring_all_MAF$names
xtx_spring_all_MAF_chromosome = xtx_spring_all_MAF$CHR
# Substitute TE names:
fbti_vector_spring <- grep("FBti", xtx_spring_all_MAF_names)
te_vector_spring <- vector()
for (i in fbti_vector_spring) {
  te_num_spring = grep(xtx_spring_all_MAF_names[i],te_dict$V1)
  te_vector_spring <- c(te_vector_spring, te_num_spring)
}
te_subs_spring <- te_dict$V2[te_vector_spring]
chrom_subs_spring <- te_dict$V3[te_vector_spring]
xtx_spring_all_MAF_names_changed <- replace(as.character(xtx_spring_all_MAF_names), fbti_vector_spring ,as.character(te_subs_spring))
xtx_spring_all_MAF$names <- xtx_spring_all_MAF_names_changed
xtx_spring_all_MAF_chromosome_changed <- replace(as.character(xtx_spring_all_MAF_chromosome), fbti_vector_spring ,as.character(chrom_subs_spring))
xtx_spring_all_MAF$chromosome <- xtx_spring_all_MAF_chromosome_changed

p3 <- ggman(xtx_spring_all_MAF, snp = "names", bp = "POS", chrom = "chromosome", pvalue = "star", logTransform = FALSE, ymax = 200, sigLine = NA, xlab = "", ylab = "XtX* values", title = "Europe Summer")
highlighted_spring <- ggmanHighlight(p3, highlight = to_highlight_spring, colour = "steelblue1", size = 0.8)
inversions_spring <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/spring30/spring_ALL_inv.txt")
inversions_spring_def <- as.vector(inversions_spring[,1])
spring30 <- ggmanHighlight(highlighted_spring, highlight = inversions_spring_def, colour = "steelblue4", size = 0.8)

TEs_spring <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/for_manhattan_plot/spring30_all_TEs_coordinates_for_plot.txt")
TEs_spring_def <- as.vector(TEs_spring[,1])
spring30_all_TE <- ggmanHighlight(spring30, highlight = TEs_spring_def, colour = "black", size = 0.8)

TEs_spring_sig <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_sept_2020_std_corrected_variables/spring30/xtx/TEs/spring30_significant_TEs_1percent.txt")
TEs_spring_sig_def <- as.vector(TEs_spring_sig[,1])
spring30_TE <- ggmanHighlight(spring30_all_TE, highlight = TEs_spring_sig_def, colour = "red", size = 1)

label_file_spring <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/spring30/spring30_labels_for_plot.txt")
spring30_TE_labels <- ggmanLabel(spring30_TE, labelDfm = label_file_spring, snp = "V1", label = "V2")


####### FALL #######
fall30.original.xtx=read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/fall30/XtX/star.xtx/fall30_all_concatenated_summary_pi_xtx.txt", h=F)
fall30.names=read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/fall/30samples_batch/splitted_autosomes/fall30_droseu14_auto.map.NAMES.ALL.txt")
fall30.original.xtx$CHROM = fall30.names$V1
fall30.original.xtx$POS = fall30.names$V2
fall30.xtx.star=read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/fall30/XtX/star.xtx/fall.star.xtx.txt",h=T)
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
fall30.original.xtx.MAF$qvalues = qvalues.fall30.MAF$qvalues

# X Chromosome
fall30.original.xtx.X=read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/fall30/XtX/star.xtx/fall30Xsize_summary_pi_xtx.out", h=T)
fall30.X.names=read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/only_variants/fall/30samples_batch/fall30_droseu14_X_names_SNPs.txt")
fall30.original.xtx.X$CHR = fall30.X.names$V1
fall30.original.xtx.X$POS = fall30.X.names$V2
fall30.xtx.star.X=read.table("/Users/pogo/Documents/Maria/Baypass/DrosEU/analysis/fall30/XtX/star.xtx/fall30.xtx.star.X.txt",h=T)
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
fall30.original.xtx.X.MAF$qvalues = qvalues.fall.X.MAF$qvalues

# Manhattan plot:
xtx_europe_fall_significant <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/fall30/fall30_all_significant.txt")
to_highlight_fall <- as.vector(xtx_europe_fall_significant[,1])

## For removing MAF from the original dataset:
xtx_fall_auto <- fall30.original.xtx.MAF[c(8,9,6,11)]
xtx_fall_X <- fall30.original.xtx.X.MAF[c(8,9,6,11)]
colnames(xtx_fall_auto) <- colnames(xtx_fall_X)
xtx_fall_MAF <- rbind(xtx_fall_auto,xtx_fall_X)
xtx_fall_all_MAF = transform(xtx_fall_MAF, names=paste(CHR, POS, sep=""))
xtx_fall_all_MAF_names = xtx_fall_all_MAF$names
xtx_fall_all_MAF_chromosome = xtx_fall_all_MAF$CHR
# Substitute TE names:
fbti_vector_fall <- grep("FBti", xtx_fall_all_MAF_names)
te_vector_fall <- vector()
for (i in fbti_vector_fall) {
  te_num_fall = grep(xtx_fall_all_MAF_names[i],te_dict$V1)
  te_vector_fall <- c(te_vector_fall, te_num_fall)
}
te_subs_fall <- te_dict$V2[te_vector_fall]
chrom_subs_fall <- te_dict$V3[te_vector_fall]
xtx_fall_all_MAF_names_changed <- replace(as.character(xtx_fall_all_MAF_names), fbti_vector_fall ,as.character(te_subs_fall))
xtx_fall_all_MAF$names <- xtx_fall_all_MAF_names_changed
xtx_fall_all_MAF_chromosome_changed <- replace(as.character(xtx_fall_all_MAF_chromosome), fbti_vector_fall ,as.character(chrom_subs_fall))
xtx_fall_all_MAF$chromosome <- xtx_fall_all_MAF_chromosome_changed

p4 <- ggman(xtx_fall_all_MAF, snp = "names", bp = "POS", chrom = "chromosome", pvalue = "star", logTransform = FALSE, ymax = 120, sigLine = NA, xlab = "Chromosomes", ylab ="XtX* values", title = "Europe Fall")
highlighted_fall <- ggmanHighlight(p4, highlight = to_highlight_fall, colour = "steelblue1", size = 0.8)
inversions_fall <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/fall30/fall_ALL_inv_v6.txt")
inversions_fall_def <- as.vector(inversions_fall[,1])
fall30 <- ggmanHighlight(highlighted_fall, highlight = inversions_fall_def, colour = "steelblue4", size = 0.8)

TEs_fall <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/for_manhattan_plot/fall30_all_TEs_coordinates_for_plot.txt")
TEs_fall_def <- as.vector(TEs_fall[,1])
fall30_all_TE <- ggmanHighlight(fall30, highlight = TEs_fall_def, colour = "black", size = 0.8)

TEs_fall_sig <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_sept_2020_std_corrected_variables/fall30/results/xtx/TEs/fall30_significant_TEs_1percent.txt")
TEs_fall_sig_def <- as.vector(TEs_fall_sig[,1])
fall30_TE <- ggmanHighlight(fall30_all_TE, highlight = TEs_fall_sig_def, colour = "red", size = 1)
label_file_fall <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/fall30/fall30_labels_for_plot.txt")
fall30_TE_labels <- ggmanLabel(fall30_TE, labelDfm = label_file_fall, snp = "V1", label = "V2")






require(gridExtra)
manhattan_xtx <- grid.arrange(all30_TE_labels, NA_TE_labels, spring30_TE_labels, fall30_TE_labels, nrow=4)
png(manhattan_xtx, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/for_manhattan_plot/manhattan_xtx.png")




require(gridExtra)
manhattan_xtx <- grid.arrange(all30_TE_labels, NA_TE_labels, spring30_TE_labels, fall30_TE_labels, nrow=4)
png(manhattan_xtx, "/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/for_manhattan_plot/manhattan_xtx.png")



