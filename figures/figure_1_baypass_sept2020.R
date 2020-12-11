# required packages 
library(raster); library(rasterVis); library(rworldxtra); data(countriesHigh)
library(rgdal)


# Read raster files
#period='1986-2010'
#r <- raster(paste('KG_', period, '.grd', sep=''))
r <- raster("/Users/pogo/Documents/Maria/Baypass/Map_KG-Global/KG_1986-2010.grd", sep='')


# Color palette for climate classification
climate.colors=c("#960000", "#FF0000", "#FF6E6E", "#FFCCCC", "#CC8D14", "#CCAA54", "#FFCC00", "#FFFF64", "#007800", "#005000", "#003200", "#96FF00", "#00D700", "#00AA00", "#BEBE00", "#8C8C00", "#5A5A00", "#550055", "#820082", "#C800C8", "#FF6EFF", "#646464", "#8C8C8C", "#BEBEBE", "#E6E6E6", "#6E28B4", "#B464FA", "#C89BFA", "#C8C8FF", "#6496FF", "#64FFFF", "#F5FFFF")

# Legend must correspond to all climate classes, insert placeholders
r0 <- r[1:32]; r[1:32] <- seq(1,32,1)

# Converts raster field to categorical data
r <- ratify(r); rat <- levels(r)[[1]]

# Legend is always drawn in alphabetic order
rat$climate <- c('Af', 'Am', 'As', 'Aw', 'BSh', 'BSk', 'BWh', 'BWk', 'Cfa', 'Cfb','Cfc', 'Csa', 'Csb', 'Csc', 'Cwa','Cwb', 'Cwc', 'Dfa', 'Dfb', 'Dfc','Dfd', 'Dsa', 'Dsb', 'Dsc', 'Dsd','Dwa', 'Dwb', 'Dwc', 'Dwd', 'EF','ET', 'Ocean')

# Remove the placeholders
r[1:32] <- r0; levels(r) <- rat

if(.Platform$OS.type=="windows") {quartz<-function(...) windows(...)}
quartz(width=13, height=10, dpi=100)
print(levelplot(r, col.regions=climate.colors, xlab="", ylab="", 
                scales=list(x=list(limits=c(xmin(r), xmax(r)), at=seq(xmin(r), xmax(r), xat)), 
                            y=list(limits=c(ymin(r), ymax(r)), at=seq(ymin(r), ymax(r), yat))), colorkey=list(space="top", tck=0, maxpixels=ncell(r)))
      +layer(sp.polygons(countriesHigh, lwd=0.25)))
out=paste('/Users/mbogaerts/Documents/Lab/scripts/world_plot.pdf', sep='')
dev.copy2pdf(file=out)

# Select region (Australia)
# x1=80; x2=180; y1=-50; y2=20; xat=5; yat=5	
# Select region (Europe)
# x1=-20; x2=80; y1=30; y2=75; xat=5; yat=5		
# Select region (US)
# x1=-130; x2=-60; y1=20; y2=60; xat=5; yat=5
# Select region (Global)
#x1=-180; x2=180; y1=-90; y2=90; xat=20; yat=10

####### EUROPE ####### 

x1=-20; x2=50; y1=30; y2=75; xat=5; yat=5		
r <- crop(r, extent(x1, x2, y1, y2))

#https://stackoverflow.com/questions/41949694/how-to-overlay-the-cernatin-point-value-in-levelplot-using-r
eu_info <- read.table("/Users/pogo/Documents/Laptop/9-1-20/BayPass/europe_pop.txt", sep="\t", h=T)
eu_points <- data.frame(lat = eu_info$Lat, lon=eu_info$Lon)
coordinates(eu_points) <- ~lon+lat
crs(eu_points) <- projection(r)

# Visualization		
if(.Platform$OS.type=="windows") {quartz<-function(...) windows(...)}
quartz(width=13, height=10, dpi=100)

print(levelplot(r, col.regions=climate.colors, xlab="", ylab="", 
                scales=list(x=list(limits=c(xmin(r), xmax(r)), at=seq(xmin(r), xmax(r), xat)), 
                            y=list(limits=c(ymin(r), ymax(r)), at=seq(ymin(r), ymax(r), yat))), colorkey=list(space="top", tck=0, maxpixels=ncell(r)))
      +layer(sp.polygons(countriesHigh, lwd=0.25)) + layer(sp.points(eu_points, cex=2.2, pch = 23, col = "black", fill="yellow")))

out=paste('/Users/pogo/Documents/Maria/Baypass/paper/figures/EU_plot_dots_yellow_sept2020.pdf', sep='')
dev.copy2pdf(file=out)


####### NORTH AMERICA ####### 
r <- raster("/Users/pogo/Documents/Maria/Baypass/Map_KG-Global/KG_1986-2010.grd", sep='')


# Color palette for climate classification
climate.colors=c("#960000", "#FF0000", "#FF6E6E", "#FFCCCC", "#CC8D14", "#CCAA54", "#FFCC00", "#FFFF64", "#007800", "#005000", "#003200", "#96FF00", "#00D700", "#00AA00", "#BEBE00", "#8C8C00", "#5A5A00", "#550055", "#820082", "#C800C8", "#FF6EFF", "#646464", "#8C8C8C", "#BEBEBE", "#E6E6E6", "#6E28B4", "#B464FA", "#C89BFA", "#C8C8FF", "#6496FF", "#64FFFF", "#F5FFFF")

# Legend must correspond to all climate classes, insert placeholders
r0 <- r[1:32]; r[1:32] <- seq(1,32,1)

# Converts raster field to categorical data
r <- ratify(r); rat <- levels(r)[[1]]

# Legend is always drawn in alphabetic order
rat$climate <- c('Af', 'Am', 'As', 'Aw', 'BSh', 'BSk', 'BWh', 'BWk', 'Cfa', 'Cfb','Cfc', 'Csa', 'Csb', 'Csc', 'Cwa','Cwb', 'Cwc', 'Dfa', 'Dfb', 'Dfc','Dfd', 'Dsa', 'Dsb', 'Dsc', 'Dsd','Dwa', 'Dwb', 'Dwc', 'Dwd', 'EF','ET', 'Ocean')

# Remove the placeholders
r[1:32] <- r0; levels(r) <- rat

x3=-90; x4=-60; y3=20; y4=45; xat=5; yat=5		
r <- crop(r, extent(x3, x4, y3, y4))

na_info <- read.table("/Users/pogo/Documents/Maria/Baypass/paper/figures/sept_2020_data/NA_11pop_machado.txt", sep="\t", h=T)
na_points <- data.frame(lat = na_info$Lat, lon=na_info$Lon)
coordinates(na_points) <- ~lon+lat
crs(na_points) <- projection(r)

# Visualization		
if(.Platform$OS.type=="windows") {quartz<-function(...) windows(...)}
quartz(width=13, height=10, dpi=100)

print(levelplot(r, col.regions=climate.colors, xlab="", ylab="", 
                scales=list(x=list(limits=c(xmin(r), xmax(r)), at=seq(xmin(r), xmax(r), xat)), 
                            y=list(limits=c(ymin(r), ymax(r)), at=seq(ymin(r), ymax(r), yat))), colorkey=list(space="top", tck=0, maxpixels=ncell(r)))
      +layer(sp.polygons(countriesHigh, lwd=0.25)) + layer(sp.points(na_points, cex=2.2, pch = 23, col = "black", fill="yellow")))

out=paste('/Users/pogo/Documents/Maria/Baypass/paper/figures/NA_plot_yellow_11pop.pdf', sep='')
dev.copy2pdf(file=out)



