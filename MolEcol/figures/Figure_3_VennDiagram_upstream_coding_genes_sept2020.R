### Created 17th January 2020
### ANALYSIS REDONE ON 19th September 2020

# Ven Diagram

library("VennDiagram")

# NA_genes <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/genes/NA_genes_list_complete.txt")
# Eu_genes <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/genes_top_0.05/all30_all_xtx_top0.05_genes_list.bed")
# Fall_genes <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/genes_top_0.05/fall30_all_xtx_top0.05_genes_list.bed")
# Spring_genes <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/genes_top_0.05/spring30_all_xtx_top0.05_genes_list.bed")

# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)

# #### Two and two (20-01-20)
# myCol <- brewer.pal(3, "Pastel1")
# myCol2 <- brewer.pal(3, "Pastel2")
# 
# venn_NA_Eu <- venn.diagram(
#   x = list(NA_genes[,1], Eu_genes[,1]),
#   category.names = c("North America" , "Europe"),
#   #filename = '/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/genes_top_0.05/venn_diagram_EU_NA.png',
#   #output=TRUE,
#   filename = NULL,
#   output=FALSE,
#   # Numbers
#   cex = 1,
#   fontface = "bold",
#   fontfamily = "sans",
#   # Circles
#   lwd = 2,
#   lty = 'blank',
#   fill = c("#440154ff", '#21908dff'),
# 
#   # Set names
#   cat.cex = 1,
#   cat.fontface = "bold",
#   cat.default.pos = "outer",
#   cat.pos = c(-27, 27),
#   # cat.dist = c(0.055, 0.055),
#   cat.fontfamily = "sans",
#   rotation = 1
# )
# 
# venn_EuF_EuS <- venn.diagram(
#   x = list(Fall_genes[,1], Spring_genes[,1]),
#   category.names = c("Europe Fall" , "Europe Spring"),
#   #filename = '/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/genes_top_0.05/venn_diagram_EuF_EuS.png',
#   #output=TRUE,
#   filename = NULL,
#   output=FALSE,
#   # Numbers
#   cex = 1,
#   fontface = "bold",
#   fontfamily = "sans",
#   # Circles
#   lwd = 2,
#   lty = 'blank',
#   fill = c("darksalmon", 'darkseagreen1'),
#   
#   # Set names
#   cat.cex = 1,
#   cat.fontface = "bold",
#   cat.default.pos = "outer",
#   cat.pos = c(-27, 27),
#   # cat.dist = c(0.055, 0.055),
#   cat.fontfamily = "sans"
#   # rotation = 1
# )
# 
# require(gridExtra)
# venn_diagrams <- grid.arrange(venn_NA_Eu, venn_EuF_EuS)
# library(gridExtra)
# venn_diagrams <-grid.arrange(gTree(children=venn_NA_Eu), gTree(children=venn_EuF_EuS), ncol=2)
# venn_diagrams <-grid.arrange(gTree(children=venn_NA_Eu), gTree(children=venn_EuF_EuS), nrow=2)
# 
#  #### All four in the same diagram 
# myCol <- brewer.pal(4, "Pastel2")
# venn.diagram(
#   x = list(NA_genes[,1], Eu_genes[,1], Fall_genes[,1], Spring_genes[,1]),
#   category.names = c("North America" , "Europe" , "Europe Fall", "Europe Spring"),
#   filename = '/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019/genes/venn_diagram.png',
#   output=TRUE,
#   
#   # # Output features
#   # imagetype="png" ,
#   # height = 480 , 
#   # width = 480 , 
#   # resolution = 300,
#   # compression = "lzw",
#   
#   # Circles
#   lwd = 2,
#   lty = 'blank',
#   fill = myCol,
#   
#   # Numbers
#   cex = .6,
#   fontface = "bold",
#   fontfamily = "sans"
#   
#   # # Set names
#   # cat.cex = 0.6,
#   # cat.fontface = "bold",
#   # cat.default.pos = "outer",
#   # cat.pos = c(-27, 27, 135),
#   # cat.dist = c(0.055, 0.055, 0.085),
#   # cat.fontfamily = "sans",
#   # rotation = 1
# )


##### Venn diagrams for BF only with the 4 variables with more genes:

# Europe:
# temp_eu <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/snpeff/all30/only_genes/temp_all30_upstream_coding_gene_list.txt")
# #rain_eu <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/all30/genes/rain_europe_gene_list.txt")
# evap_eu <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/snpeff/all30/only_genes/evap_all30_upstream_coding_gene_list.txt")
# solar_eu <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/snpeff/all30/only_genes/solar_all30_upstream_coding_gene_list.txt")
# #soil
# wind_eu <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/snpeff/all30/only_genes/wind_all30_upstream_coding_gene_list.txt")
# #daylight_eu <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/all30/genes/daylight_europe_gene_list.txt")

temp_eu<- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_sept_2020_std_corrected_variables/all30/results/genes/all_genes/temperature/all30_temperature_flybase_snpeff_genes_list.txt")
solar_eu <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_sept_2020_std_corrected_variables/all30/results/genes/all_genes/solar_radiation/all30_solarrad_flybase_snpeff_genes_list.txt")
wind_eu <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_sept_2020_std_corrected_variables/all30/results/genes/all_genes/wind/all30_wind_flybase_snpeff_genes_list.txt")

#myCol <- brewer.pal(3, "Pastel2")
myCol <- c("#B3E2CD","#F4CAE4", "#E6F5C9")

venn_Eu_BF <- venn.diagram(
  #x = list(temp_eu[,1], rain_eu[,1], evap_eu[,1], solar_eu[,1], wind_eu[,1], daylight_eu[,1]),
  #category.names = c("Temperature" , "RainFall", "Evaporation", "Solar Raiation", "Wind", "Daylight hours"),
  x = list(temp_eu[,1], solar_eu[,1], wind_eu[,1]),
  category.names = c("Temperature" , "Solar Radiation", "Wind"),
  #filename = '/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/plots/venn_diagram_EU_BF_3var__upstream_coding.png',
  filename = '/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_sept_2020_std_corrected_variables/plots/venn_diagram_EU_BF_3var_upstream_coding.png',
  output=TRUE,
  #height = 1200,
  #width = 1200,
  # Numbers
  cex = 2,
  fontface = "bold",
  fontfamily = "sans",
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Set names
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  #cat.pos = c(11,6,6),
  #cat.dist = 1,
  cat.fontfamily = "sans",
  margin = 0.1,
  # rotation = 1
  
  # # Title
  # main = "A) Europe",
  # main.cex = 2,
  # main.fontfamily = "sans",
  # main.fontface = "bold"
)


# North America:
## Changed on 20th September 2020
#temp_na <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_february_2020_standard/all/temp_ALL_NA_gene_list.txt")
#rain_na <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/snpeff/NA/only_genes/rain_NA_upstream_coding.txt")
#evap_na <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_february_2020_standard/all/evap_ALL_NA_gene_list.txt")
#solar_eu <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/all30/genes/solar_europe_gene_list.txt")
#soil
#wind_na <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_february_2020_standard/all/wind_ALL_NA_gene_list.txt")
#daylight_eu <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/all30/genes/daylight_europe_gene_list.txt")

temp_na <- read.table("/Users/pogo/Documents/Maria/Baypass/NA/analysis/standard/results/genes/all_genes/temperature_machadoNA_flybase_snpeff_genes_list.txt")
rain_na <- read.table("/Users/pogo/Documents/Maria/Baypass/NA/analysis/standard/results/genes/all_genes/rainfall_machadoNA_flybase_snpeff_genes_list.txt")
solar_na <- read.table("/Users/pogo/Documents/Maria/Baypass/NA/analysis/standard/results/genes/all_genes/solarrad_machadoNA_flybase_snpeff_genes_list.txt")


#myCol <- brewer.pal(3, "Pastel2")
myCol <- c("#B3E2CD", "#FDCDAC", "#F4CAE4")

venn_NA_BF <- venn.diagram(
  #x = list(temp_eu[,1], rain_eu[,1], evap_eu[,1], solar_eu[,1], wind_eu[,1], daylight_eu[,1]),
  #category.names = c("Temperature" , "RainFall", "Evaporation", "Solar Raiation", "Wind", "Daylight hours"),
  x = list(temp_na[,1], rain_na[,1], solar_na[,1]),
  category.names = c("Temperature", "RainFall", "Solar Radiation"),
  #filename = '/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/plots/venn_diagram_NA_BF_3var_v5_upstream_coding_corrected_TEs.png',
  filename = '/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_sept_2020_std_corrected_variables/plots/venn_diagram_NA_BF_3var_upstream_coding.png',
  output=TRUE,
  
  # # Output features
  # imagetype="png" ,
  # height = 2000 ,
  # width = 2200 ,
  # resolution = 300,
  # compression = "lzw",
  
  # Numbers
  cex = 2,
  fontface = "bold",
  fontfamily = "sans",
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Set names
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  #cat.pos = c(11,6,6),
  #cat.dist = 1,
  cat.fontfamily = "sans",
  margin = 0.1,
  # rotation = 1
  
  # # Title
  # main = "North America"
)

# Europe Spring:

#temp_sp <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/snpeff/spring/only_genes/temp_ALL_gene_list_spring.txt")
#rain_sp <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/snpeff/spring/only_genes/rain_ALL_gene_list_spring.txt")
#evap_sp <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/snpeff/spring/only_genes/evap_ALL_gene_list_spring.txt")
#solar_sp <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/snpeff/spring/only_genes/solar_ALL_gene_list_spring.txt")
#soil
#wind_na <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/genes/NA/wind_NA_gene_list.txt")
#daylight_eu <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/all30/genes/daylight_europe_gene_list.txt")

temp_sp <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_sept_2020_std_corrected_variables/spring30/results/genes/all_genes/temperature_spring30_flybase_snpeff_genes_list.txt")
rain_sp <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_sept_2020_std_corrected_variables/spring30/results/genes/all_genes/rainfall_spring30_flybase_snpeff_genes_list.txt")
evap_sp <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_sept_2020_std_corrected_variables/spring30/results/genes/all_genes/evaporation_spring30_flybase_snpeff_genes_list.txt")

#myCol <- brewer.pal(3, "Pastel2")
myCol <- c("#B3E2CD","#FDCDAC", "#CBD5E8")

venn_SP_BF <- venn.diagram(
  #x = list(temp_eu[,1], rain_eu[,1], evap_eu[,1], solar_eu[,1], wind_eu[,1], daylight_eu[,1]),
  #category.names = c("Temperature" , "RainFall", "Evaporation", "Solar Raiation", "Wind", "Daylight hours"),
  x = list(temp_sp[,1], rain_sp[,1], evap_sp[,1]),
  category.names = c("Temperature" , "RainFall", "Evaporation"),
  #filename = '/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/plots/venn_diagram_SP_BF_3var_upstream_coding.png',
  filename = '/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_sept_2020_std_corrected_variables/plots/venn_diagram_SP_BF_3var_upstream_coding.png',
  output=TRUE,
  # Numbers
  cex = 2,
  fontface = "bold",
  fontfamily = "sans",
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Set names
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  #cat.pos = c(11,6,6),
  #cat.dist = 1,
  cat.fontfamily = "sans",
  margin = 0.1,
  # rotation = 1
  
  # Title
  #main = "Europe Spring"
)

# Europe Fall:
# Modified on 2nd May 2020
#temp_fl <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/snpeff/fall/only_genes/temp_fin_made_mbogaerts.txt")
#temp_fl <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/snpeff/fall/only_genes/temp_fin_made_mbogaerts.txt")
#rain_fl <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/snpeff/fall/only_genes/rain_ALL_gene_list_fall.txt")
#evap_fl <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/snpeff/fall/only_genes/evap_ALL_gene_list_fall.txt")
#solar_sp <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/spring30/genes/solar_spring_gene_list.txt")
#soil
#wind_fl <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/snpeff/fall/only_genes/wind_fall_gene_list.txt")
#daylight_eu <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/all30/genes/daylight_europe_gene_list.txt")

temp_fl <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_sept_2020_std_corrected_variables/fall30/results/genes/all_genes/temperature_fall30_flybase_snpeff_genes_list.txt")
rain_fl <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_sept_2020_std_corrected_variables/fall30/results/genes/all_genes/rainfall_fall30_flybase_snpeff_genes_list.txt")
wind_fl <- read.table("/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_sept_2020_std_corrected_variables/fall30/results/genes/all_genes/wind_fall30_flybase_snpeff_genes_list.txt")

#myCol <- brewer.pal(3, "Pastel2")
myCol <- c("#B3E2CD", "#FDCDAC","#E6F5C9")

venn_FL_BF <- venn.diagram(
  #x = list(temp_eu[,1], rain_eu[,1], evap_eu[,1], solar_eu[,1], wind_eu[,1], daylight_eu[,1]),
  #category.names = c("Temperature" , "RainFall", "Evaporation", "Solar Raiation", "Wind", "Daylight hours"),
  x = list(temp_fl[,1], rain_fl[,1], wind_fl[,1]),
  category.names = c("Temperature" , "RainFall", "Wind"),
  #filename = '/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_december_2019_standard/plots/venn_diagram_FL_BF_3var_upstream_coding.png',
  filename = '/Users/pogo/Documents/Maria/Baypass/analysis/reanalysis_sept_2020_std_corrected_variables/plots/venn_diagram_FL_BF_3var_upstream_coding.png',
  output=TRUE,
  # Numbers
  cex = 2,
  fontface = "bold",
  fontfamily = "sans",
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Set names
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  #cat.pos = c(11,6,6),
  #cat.dist = 1,
  cat.fontfamily = "sans",
  margin = 0.1,
  # rotation = 1
  
  # Title
  #main = "Europe Fall"
)

#require(gridExtra)
#grid.arrange(gTree(children=g), top="Title", bottom="subtitle")
#venn_diagrams <-grid.arrange(gTree(children=venn_Eu_BF), gTree(children=venn_NA_BF), gTree(children=venn_SP_BF), gTree(children=venn_FL_BF), ncol=2, nrow=2)
 
library(cowplot)
library(grid)
library(ggplotify)

all30_venn <- as.grob(venn_Eu_BF)
NA_venn <- as.grob(venn_NA_BF)
spring30_venn <- as.grob(venn_SP_BF)
fall30_venn <- as.grob(venn_FL_BF)


#650x1200
plot_grid(venn_Eu_BF, venn_NA_BF, venn_SP_BF, venn_FL_BF, nrow = 2, ncol = 2, labels="AUTO")
#plot_grid(all_auto_hm, all_X_hm, NA_auto_hm, NA_X_hm, spring_auto_hm, spring_X_hm, fall_auto_hm, fall_X_hm, nrow = 4, ncol = 2, labels=c("Europe Autosomes", "Europe X chrom", "NA Autosomes", "NA X chrom", "Europe Spring Autosomes", "Europe Spring X chrom", "Europe Fall Autosomes", "Europe Fall X chrom"))



