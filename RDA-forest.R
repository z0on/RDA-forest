rm(list=ls())
# edit this line to point to the path you downloaded the files to:
setwd("~/Dropbox/RDA-forest")

#install.packages("~/Dropbox/RDAforest_2.0.0.tar.gz")
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
map(coasts)

# loading all data: 
ll=load("Agaricia_yellow.RData")
ll
# latlon: geographical coordinates of samples
# covars: covariates to regress out of the genetic data before computing ordination.
#         Most commonly these are technical things affecting sequencing quality (i.e. sequencing batch)
#         here, covariates are admixture proportions with the other two cryptic lineages existing within the same species
# IBS:    square matrix of genetic distances. In this case, this is identity-by-state matrix produced by ANGSD based on 2bRAD data.
# env:    environmental variables to test for association with coral genetics

plot(lat~lon,latlon,pch=19,cex=1,col="grey70",asp=1)
map(coasts,add=T,col="grey80",fill=T,border="grey80",lwd=1)

# ------ hierarchical clustering tree
# to confirm there are no highly related (possibly clonal) samples (forming forks "hanging down" much lower than others)
# or wrong species collected (well-defined highly dissimilar minor clusters)

hc=hclust(as.dist(IBS),"ave")
plot(hc,cex=0.5)
# if there are clones, remove all but one representative of each clonal group;
# remove wrong species entirely.

# --------------- spatial autocorrelation
# is there correlation between geographic and genetic distances (isolation-by-distance)?
# If yes, spatial predictors must be included in the model
# running procrustes test
protest(capscale(dist(latlon)~1),capscale(IBS~1))
# looks like we have very strong IBD! correlation is almost 0.8
# lets plot genetic distances against geographic, to convince ourselves
plot(as.dist(IBS)~dist(latlon),col=rgb(0,0,0,alpha=0.2),log="xy")

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
# Note #1:  Y can be a square distance matrix (like we use here)
# or raw measures matrix (rows:samples, columns:features)
# use nreps=11 for exploratory analysis, nreps=25 (default) for real
# adjust top.pcs according to the looks of your eigenvalue plot
# Note #2: we are also adding spatial variables (lat, lon) here
# environmental variables must not get "beaten" by those during variable selection stage
# For random forest, just the geographical coordinates are 100% sufficient 
# to model any spatial grouping of points, so we will not bother with Moran Eigenvector Maps
# Note #3: we also supply covariates to be regressed out of IBS before ordination. 

mm=mtrySelJack(Y=IBS,X=cbind(latlon,env),covariates=covars,nreps=15,top.pcs=15)

# which environmental variables pass the selection process?
mm$goodvars

# boxplot of importance differences at higher vs lower mtry
ggplot(mm1$delta,aes(var,values))+
  geom_boxplot()+
  coord_flip()+
  geom_hline(yintercept=0,col="red")

# bar chart of proportion of positive change in response to higher mtry
# good predictors would be the ones above the red line
ggplot(mm$prop.positive,aes(var,prop.positive))+
  geom_bar(stat="identity")+
  coord_flip()+
  geom_hline(yintercept=0.6,col="red")

#save(mm,file="mm.RData")
#load("mm.RData")

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

# if there are no new points to predict, just skip the newX option in the call to ordinationJackknife;
# the predictions will be made for original data points then.

# non-jackknifed model to see how each variable affects each MDS
sb1=makeGF(ord.all,X=cbind(latlon,env)[,mm$goodvars],keep=c(1:15),ntrees=1500)
# how much each MDS is explained overall
sb1$result

# plot turnover curves for the top 4 predictors plus lat, lon
important=names(importance(sb1))[1:4]
plot(sb1,plot.type="C",show.species=T,common.scale=T,imp.vars=important,
     #     cex.axis=0.6,cex.lab=0.7,
     line.ylab = 0.9, par.args = list(
       mgp = c(1.5,0.5, 0), 
       mar = c(2.5, 1, 0.1, 0.5), 
       omi = c(0,0.3, 0, 0)
     )
)

names(XY)=c("lon","lat")

evars=mm$goodvars[!(mm$goodvars %in% c("lat","lon"))]

# running full modeling/prediction with jackknifing, including lat, lon
sb=ordinationJackknife(Y=IBS,X=cbind(latlon,env)[,mm$goodvars],newX=cbind(XY,rasters),oob=0.2,nreps=25,covariates=covars,top.pcs=15)

# same without lat, lon - it woudl be the same in this case since lat, lon did not make it through mtry selection
sb.env=ordinationJackknife(Y=IBS,X=env[,evars],newX=rasters,oob=0.1,nreps=5,covariates=covars,top.pcs=15)

# just lat, lon
sb.latlon=ordinationJackknife(Y=IBS,X=latlon,newX=XY,nreps=15,covariates=covars,top.pcs=15)

# overall R2 of these models
sum(sb$median.importance)
sum(sb.env$median.importance)
sum(sb.latlon$median.importance)

# importance boxplot
ggplot(sb$all.importances,aes(variable,importance))+geom_boxplot()+coord_flip()+theme_bw()

#save(mm,sb,env,latlon,file="RDAFresult_aagaY_latlon.RData")
#load("RDAFresult_aagaY_latlon.RData")

#------ plotting adaptation map

evars=mm$goodvars[!(mm$goodvars %in% c("lat","lon"))]

# predicting genetic PCs straight-up, on full dataset without jackknifing
pcs=scores(ord.all,scaling=1,display="sites",choices=c(1:15))
rfs=predict_rf(Y=pcs,X=env[,evars],newX=rasters)

# alternatively, we can use predictions averaged over ordinationJackknife replicates 
# (not always advantageous, try it either way)
rfs=sb$predictions

# making new ordination and fitting environmental vectors
rd=rda(rfs~1)
bests=names(sb$median.importance)
ee=envfit(rd,cbind(XY,rasters)[,bests])

pp=plot_adaptation(rd,ee,XY,flip1=1,flip2=1,nclust=3,scal=20,cex=0.8,jitscale=0.03,rangeExp=2,coltxt = "coral",cluster.indices = TRUE)
# to make a color-clustered version of the map and PCA, set nclust to some number of clusters 
# to try different color scheme try setting flip1 and/or flip2 to -1
# to make envfit arrows longer, reduce scal
# to make PCA smaller (and leave more space for envfit arrows and labels) increase rangeExp
# cex controls size of arrow text labels
# to make labels sty further away from arrows increase jitscale
# (the distance of label from its arrow in a bit different for each label every time you plot,
# so if labels overlap try plotting the same again until there is a version where they don't)

# # adding map
maps::map(coasts, add = TRUE)
