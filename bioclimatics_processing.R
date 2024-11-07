#install.packages("geodata")
#install.packages("terra")

library(geodata)
library(terra)

# Download global bioclimatic data from worldclim (you may have to set argument 'download=T' for first download, if 'download=F' it will attempt to read from file):
clim <- geodata::worldclim_global(var = 'bio', res = 10, download = F, path = 'data')
# cropping to North America
clim=crop(clim,ext(-166.833,-52.5,40.1666,79.8333))

# terra package offers different functionalities to manipulate the spatial data, 
# for example aggregating the data to coarser resolutions (aggregate), 
# cropping (crop()), and adding spatial layers to a SpatRaster object (c()):
plot(terra::aggregate(clim[[2]], fact=2, fun="mean"))

# Download future climate scenario from 'ACCESS-ESM1-5' climate model.
# Please note that you have to set download=T if you haven't downloaded the data before:
clim_fut <- geodata::cmip6_world(model='ACCESS-ESM1-5', ssp='245', time='2041-2060', var='bioc', download=F, res=10, path='data')
# normalizing names
clim_fut=crop(clim_fut,ext(-166.833,-52.5,40.1666,79.8333))

# renaming biovariables
bionames=scan("~/Dropbox/numerical_ecology/Forester_Wolf_Rcode/bioclimatic_variables_list_abbreviated.txt",what="c")
bionames[14]="Prec_DryM"
bionames[17]="Prec_DryQ"
bionames[3]="Isotherm"
bionames[7]="Ann_TRange"
names(clim)= names(clim_fut)=bionames

# Download fractional tree cover at 30-sec resolution:
# Please note that you have to set download=T if you haven't downloaded the data before:
trees_30sec <- geodata::landcover(var='trees', path='data', download=F)
trees_30sec=crop(trees_30sec,ext(-166.833,-52.5,40.1666,79.8333))

# Aggregate tree cover to 5-min spatial resolution
trees_10min <- terra::aggregate(trees_30sec, fact=20, fun='mean')

# Produce the multi-layer environmental data object with matching extents:
env_cur <- c(clim, trees_10min)
env_fut <- c(clim_fut, trees_10min)

# reducing resolution
env_c <- terra::aggregate(env_cur, fact=3, fun='mean')
env_f <- terra::aggregate(env_fut, fact=3, fun='mean')

# converting rasters to dataframes
envc=as.data.frame(env_c,xy=TRUE)
isna=apply(envc,1,mean)
envc=envc[!is.na(isna),]
envf=as.data.frame(env_f,xy=TRUE)
isna=apply(envf,1,mean)
envf=envf[!is.na(isna),]

# loading wolf data, including latlon - the coordinates of sampled wolves
ll=load("~/Dropbox/RDA-forest/wolf_v4.RData")

# moving one point which is in the sea
latlon[68,]=latlon[68,]+0.5

# Extracting environmental data for a given set of coordinates.
env=terra::extract(x = env_cur, y = data.frame(latlon), cells=F, method="bilinear" )
# checking for misplaced points
which(is.na(env[,2]))
# removing index column ans saving
env=env[,-1]
save(geno,env,latlon,envc,envf,ecotype,file="~/Dropbox/RDA-forest/wolf_v4.RData")
