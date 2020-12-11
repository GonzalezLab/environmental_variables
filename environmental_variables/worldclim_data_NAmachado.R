#  Raster-small tutorial.sh
#  
#
#  Created by Martin Kapun on 21/11/16.
#

## if not installed yet, install raster package:

install.packages("raster")

## load installed package to the R project

require(raster)

## first load WorldClim bio variables at the resolution of 2.5 degrees (see Hijams et al. 2005 for more details)

# The "bio" dataset contains annual averages of 19 interpolated biovaribles
biod <- getData('worldclim', var='bio', res=2.5)
# The "tmin" dataset contains monthly averages of minimum temperatures
tmind <- getData('worldclim', var='tmin', res=2.5)
# The "tmax" dataset contains monthly averages of maximum temperatures
tmaxd <- getData('worldclim', var='tmax', res=2.5)
# The "precd" dataset contains monthly precipitation averages
precd <- getData('worldclim', var='prec', res=2.5)

# now read your input data, which contains the coordinates for the positions you are interested in.

geod<-read.table('/Users/mbogaerts/Lab/BayPass/env_var/NA/worldclim_coord_NA.txt', header=T, stringsAsFactors=F)
droseu<-read.table('/Users/mbogaerts/Lab/BayPass/env_var/DrosEU_coord.txt', header=T, stringsAsFactors=F)

# extact climate data for each pair of Longitude/Latitude coordinates
bio<-extract(biod, geod[,c(2,3)])
tmin<-extract(tmind, geod[,c(2,3)])
tmax<-extract(tmaxd, geod[,c(2,3)])
precd<-extract(precd, geod[,c(2,3)])

# combine coordinates and climate variables
bio.data<-cbind(geod,bio,tmin,tmax,precd)

# save into external file
write.table(bio.data,file='/Users/mbogaerts/Lab/BayPass/env_var/NA/worldclim_coord_NA_clim.txt',sep="\t", row.names=FALSE ,quote=FALSE)

