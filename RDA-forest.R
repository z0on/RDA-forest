rm(list=ls())
# edit this line to point to the path you downloaded the files to:
setwd("~/Dropbox/RDA-forest")

install.packages("~/Dropbox/RDAforest_1.0.1.tar.gz")
library(RDAforest)
library(ggplot2)
packageDescription("RDAforest")
theme_set(theme_bw())


# the next four are only to add outline of coasts to the last plot, at the very last line of this script.
# Skip if having trouble installing these packages.
library(maps)
library(rnaturalearth)
library(rnaturalearthdata)
coasts= ne_coastline(scale="large")

# environmental parameters + spatial coordinates (lat, lon)
env=read.csv("environment.csv")

# covariates (to regress out of the ordination)
# Most commonly these are technical things affecting sequencing quality (i.e. sequencing batch)
# here, covariates are admixture proportions with the other two cryptic lineages existing within the same species
# do NOT put spatial variables here! They should be appended to env (later)
covars=read.csv("covariates.csv")

# spatial coordinates
latlon=read.csv("latlon.csv")

# genetic distances
IBS=as.matrix(read.table("Agaricia_Y.ibsMat"))

# ------ hierarchical clustering tree
# to confirm there are no highly related (possibly clonal) samples (forming forks "hanging down" much lower than others)
# or wrong species collected (well-defined highly dissimilar minor clusters)

hc=hclust(as.dist(IBS),"ave")
plot(hc,cex=0.5)
# if there are clones, remove all but one representative of each clonal group;
# remove wrong species entirely.

# --------------- spatial autocorrelation
# is there correlation between geographic and genetic distances (isolation-by-distance)?
# If yes, spatial predictors must be included in env
# running procrustes test
protest(capscale(dist(latlon)~1),capscale(IBS~1))
# looks like we have very strong IBD!
# lets plot genetic distances against geographic, to convince ourselves
plot(as.dist(IBS)~dist(latlon),col=rgb(0,0,0,alpha=0.2),log="xy")

#  spatial predictors
#  environmental variables must not get "beaten" by those during variable selection stage
# >>>> For random forest, just the geographical coordinates are 100% sufficient 
# to model any spatial grouping of points!!!
# so we will just use lat,lon (later)
# and will not bother with Moran Eigenvector Maps

# # computing uncorrelated spatial coordinates
# xy=data.frame(scores(capscale(dist(latlon)~1),scaling=1)$sites)
# space=names(xy)=c("xx","yy")
# plot(xy,asp=1)

#------- let's see how our ordination looks: computing conditional PCoA

ord.all=capscale(as.dist(IBS)~1+Condition(as.matrix(covars)))
plot(ord.all,scaling=1)
# how many interesting PCs?
plot(ord.all$CA$eig/sum(ord.all$CA$eig))
# looks like beyond PC 10 there is just noise. Will use the first 15 PCs to be on safe side.

#--------------- predictor selection

# the mtrySelJack function will fit two models to each jackknifed replicate,
# with lower and higher mtry settings
# bad predictors ("standing-in" for actually important correlated ones) should
# decrease in raw importance at higher mtry (as per Strobl et al 2008)
# Note -  Y can be a square distance matrix (like we use here)
# or raw measures matrix (rows:samples, columns:features)
# use nreps=15 for exploratory analysis, nreps=25 (default) for real
# adjust top.pcs according to the looks of your eigenvalue plot
# Note #2: we are NOT using spatial variables (lat, lon) here,
# they will have to be included in the final model regardless of their importance.
mm=mtrySelJack(Y=IBS,X=env,covariates=covars, nreps=25, top.pcs=15)
# boxplot of importance differences at different mtry
ggplot(mm$delta,aes(var,values))+
  geom_boxplot()+
  coord_flip()+
  geom_hline(yintercept=0,col="red")

# bar chart of proportion of positive change in response to higher mtry
# good predictors would be the ones above the red line
ggplot(mm$prop.positive,aes(var,prop.positive))+
  geom_bar(stat="identity")+
  coord_flip()+
  geom_hline(yintercept=0.53,col="red")

# which environmental variables pass the selection process?
mm$goodvars

#--------------- jackknifing of mtry-passing variables

# BUT FIRST!
# loading table of environmental values for a grid of spatial points
# where we want to predict how our creatures would adapt.
# see for example Fitzpatrick and Keller 2015, Bay et al 2018 (pdfs included in repo)

ll=load("rasters_XY.RData")
ll
# "rasters" "XY"
# this is the grid (Florida Keys seascape):
plot(XY,pch=".",asp=1)

# if there are no new points to predict, just skip the newX option in the call to spatialBootstrap;
# the predictions will be made for original data points then.

# non-jackknifed model to see how each variable affects each MDS
sb1=makeGF(ord.all,X=cbind(latlon,env[,mm$goodvars]),keep=c(1:15),ntrees=1500)
# how much each MDS is explained overall
sb1$result

# plot turnover curves for the top 4 predictors plus lat, lon
important=names(importance(sb1))[1:4]
important=union(important,c("lat","lon"))
plot(sb1,plot.type="C",show.species=T,common.scale=T,imp.vars=important,
     #     cex.axis=0.6,cex.lab=0.7,
     line.ylab = 0.9, par.args = list(
       mgp = c(1.5,0.5, 0), 
       mar = c(2.5, 1, 0.1, 0.5), 
       omi = c(0,0.3, 0, 0)
     )
)

# running full modeling/prediction with jackknifing, including lat, lon
sb=ordinationJackknife(Y=IBS,X=cbind(latlon,env[,mm$goodvars]),newX=rasters,nreps=25,covariates=covars,top.pcs=15)
# same without lat, lon
sb.env=ordinationJackknife(Y=IBS,X=env[,mm$goodvars],nreps=25,covariates=covars,top.pcs=15)
# just lat, lon
sb.latlon=ordinationJackknife(Y=IBS,X=latlon,nreps=25,covariates=covars,top.pcs=15)
sum(sb$median.importance)
# 0.44
sum(sb.env$median.importance)
# 0.44
sum(sb.latlon$median.importance)
# 0.40


# importance boxplot including space variables
ggplot(sb$all.importances,aes(variable,importance))+geom_boxplot()+coord_flip()+theme_bw()
# importance boxplot without space variables
ggplot(sb$all.importances[!(sb$all.importances$variable %in% space),],aes(variable,importance))+geom_boxplot()+coord_flip()+theme_bw()

save(mm,sb,sb1,env,latlon,file="RDAFresult_aagaY_latlon.RData")

load("RDAFresult_aagaY_latlon.RData")

#------ plotting genetic turnover (adaptation) map
# the following code has been cannibalized from gradientForest vignette. It looks weird but it works!

turnovers=sb$turnovers
raster.vars=colnames(turnovers)

# principal component analysis of the predicted genetic turnover patterns
pc <- prcomp(turnovers)
plot(pc$sdev)
# we will try to visualize the first 3 PCs
pcs2show=c(1,2,3)

# color flippage flags - change between -1 and 1
# to possibly improve color representation in the final map
flip.pc1=(1)
flip.pc2=(-1)
flip.pc3=(1)

flip=""
if(flip.pc1==(-1)) { flip=paste(flip,"1",sep="")}
if(flip.pc2==(-1)) { flip=paste(flip,"2",sep="")}
if(flip.pc3==(-1)) { flip=paste(flip,"3",sep="")}

# magic to color three PCA dimensions on a map
pc1 <- flip.pc1*pc$x[, pcs2show[1]]
pc2 <- flip.pc2*pc$x[, pcs2show[2]]
pc3 <- flip.pc3*pc$x[, pcs2show[3]]
b <- pc1 - pc2
g <- -pc1
r <- pc3 + pc2 - pc1
r <- (r - min(r))/(max(r) - min(r)) * 255
g <- (g - min(g))/(max(g) - min(g)) * 255
b <- (b - min(b))/(max(b) - min(b)) * 255
nvs <- dim(pc$rotation)[pcs2show[1]]
lv <- length(raster.vars)
vind <- rownames(pc$rotation) %in% raster.vars
scal <- 30
xrng <- range(pc$x[, 1], pc$rotation[, pcs2show[1]]/scal) * 1.9
yrng <- range(pc$x[, 2], pc$rotation[, pcs2show[2]]/scal) * 1.9
man.colors=rgb(r, g, b, max = 255)

# -------- plotting the map of predicted adaptive communities

coltxt="hotpink"
important=raster.vars[order(sqrt(pc$rotation[raster.vars, 1]^2+pc$rotation[raster.vars, 2]^2),decreasing=T)][1:min(3,length(raster.vars))]
lv=length(important)
par(mfrow=c(1,2))
plot((pc$x[, pcs2show[1:2]]), xlim = xrng, ylim = yrng, pch = ".", cex = 4, col = rgb(r, g, b, max = 255), asp = 1,xaxt="none",yaxt="none",bty="none",xlab="",ylab="")
jit <- rnorm(lv,0.01,0.002)
arrows(rep(0, lv), rep(0, lv), pc$rotation[important, pcs2show[1]]/scal, pc$rotation[important, pcs2show[2]]/scal, length = 0.0625,col=coltxt)
text(pc$rotation[important, 1]/scal + jit * sign(pc$rotation[important, pcs2show[1]]), pc$rotation[important, pcs2show[2]]/scal + jit * sign(pc$rotation[important, pcs2show[2]]), labels = important,cex=0.7,col=coltxt)
plot(XY, pch=15,cex = 0.5, asp = 1, col = man.colors)
map(coasts,add=T,col="grey80",fill=T,border="grey80",lwd=1)

# contrasting colors = habitats requiring differential adaptation, likely driven by factors in the legend.

# plotting influential variables
library(viridis)
ggplot(XY,aes(x,y,color=rasters$SST_mean.coldest.month))+geom_point()+scale_color_viridis()+coord_equal()
ggplot(XY,aes(x,y,color=rasters$SAL.B_yearly_range))+geom_point()+scale_color_viridis()+coord_equal()
ggplot(XY,aes(x,y,color=rasters$Depth))+geom_point()+scale_color_viridis()+coord_equal()

