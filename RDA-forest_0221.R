rm(list=ls())
# edit this line to point to the path you downloaded the files to:
setwd("~/Dropbox/RDA-forest")

#install.packages("~/Dropbox/RDAforest_2.1.4.tar.gz")
library(RDAforest)
library(ggplot2)

theme_set(theme_bw())

# the next four are only to add outline of coasts to the last plot, at the very last line of this script.
# Skip if having trouble installing these packages.
library(maps)
library(mapdata)

# loading all data: 
ll=load("Agaricia_yellow.RData")
ll
# latlon: geographical coordinates of samples
# covars: covariates to regress out of the genetic data before computing ordination.
#         Most commonly these are technical things affecting sequencing quality (i.e. sequencing batch)
#         here, covariates are admixture proportions with the other two cryptic lineages existing within the same species
# IBS:    square matrix of genetic distances. In this case, this is identity-by-state matrix produced by ANGSD based on 2bRAD data.
# env:    environmental variables to test for association with coral genetics

# map of sampling sites
plot(lat~lon,latlon,cex=1,asp=1)
map(database="worldHires", add=TRUE,fill=T,col="grey90")

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

# automatic mtry settings: N/5, 2*N/5
mm=mtrySelJack(Y=IBS,X=cbind(latlon,env),covariates=covars,nreps=15,top.pcs=15)

# which environmental variables pass the selection process?
mm$goodvars

# boxplot of importance differences at higher vs lower mtry
ggplot(mm$delta,aes(var,values))+
  geom_boxplot()+
  coord_flip()+
  ggtitle("mtry 3, 6")+
  geom_hline(yintercept=0,col="red")
ggplot(mm$delta,aes(var,values))+
  geom_boxplot()+
  coord_flip()+
  ggtitle("default")+
  geom_hline(yintercept=0,col="red")

save(mm,file="mms_02212024.RData")
#load("mms_02212024.RData")

#--------------- jackknifing of mtry-passing variables

# BUT FIRST!
# loading table of environmental values for a grid of spatial points
# where we want to predict how our creatures would adapt.
# see for example Fitzpatrick and Keller 2015, Bay et al 2018 (pdfs included in repo)

ll=load("rasters_XY.RData")
ll
# "rasters" "XY"
# this is the grid (Florida Keys seascape):
names(XY)=c("lon","lat")
plot(XY,pch=".",asp=1,col="grey90")
map(database="worldHires", add=TRUE,fill=T,col="grey90")

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

#----- measuring importances while jackknifing the ordination

evars=mm$goodvars[!(mm$goodvars %in% c("lat","lon"))]
# environmental predictors
sb=ordinationJackknife(Y=IBS,X=env[,evars],extra=0.1,newX=rasters,nreps=15,covariates=covars,top.pcs=15)
# just coordinates (XY) as predictors
sbxy=ordinationJackknife(Y=IBS,latlon,extra=0.1,newX=XY,nreps=15,covariates=covars,top.pcs=15)

save(mm,sb,sbxy,file="sbs_02212024.RData")
load("sbs_02212024.RData")

# overall R2 of these models
sum(sb$median.importance) # 0.47
sum(sbxy$median.importance) # 0.41

# importance boxplots
ggplot(sb$all.importances,aes(variable,importance))+geom_boxplot()+coord_flip()+theme_bw()
ggplot(sbxy$all.importances,aes(variable,importance))+geom_boxplot()+coord_flip()+theme_bw()

#------ plotting adaptation map

# predicting genetic PCs straight-up, on full dataset without jackknifing
pcs=scores(ord.all,scaling=1,display="sites",choices=c(1:15))

# straight-up random forest predictions
p0=predict_rf(Y=pcs,X=env[,evars],newX=rasters,extra=0.1)

# subsetting the bigger predictor dataset to match the modeled range
table(p0$goodrows)
xy2=XY[p0$goodrows,]
rfs=p0$preds
ras2=rasters[which(p0$goodrows),]

# gradient forest predictons (turnover curves to be used as clustering guide)
p01=predict_gf(pcs,env[,evars],newX=rasters,extra=0.1)
pp=p01$preds

# making new ordination and fitting environmental vectors
bests=names(sb$median.importance)[1:3]

#---- PLOTTING adaptive neighborhoods

# tips for plot_adaptation function:
# to make a color-clustered version of the map and PCA, set nclust to some number of clusters 
# to try different color scheme try setting flip1 and/or flip2 to -1
# to make envfit arrows longer, reduce scal
# to make PCA smaller (and leave more space for envfit arrows and labels) increase rangeExp
# cex controls size of points AND text labels
# to make labels stand further away from arrows increase jitscale
# (the distance of label from its arrow will be a bit different each time you plot,
# so if labels overlap try plotting the same again and again until there is a version where they don't)

#pdf(file="agaricia_adaptationPlots_rf.pdf",height=4,width=8)

# unclustered adaptation
plot_adaptation(rfs,ras2[,bests],xy2,pcs2show=2,flip1=1,flip2=-1,scal=8,cex=0.8,jitscale=0.05,rangeExp=1.5,coltxt="gray70")
map(database="worldHires", add=TRUE,fill=T,col="grey90")
#plot(coasts,add=TRUE,border="grey60",fill=TRUE,col="grey60",lwd=1)

# clustered based on RF predictions
plot_adaptation(rfs,ras2[,bests],xy2,nclust=7,pcs2show=2,flip1=1,flip2=-1,scal=8,cex=0.8,color.by.medoids = T,jitscale=0.05,rangeExp=1.5,coltxt="gray70",cluster.labels = T)
map(database="worldHires", add=TRUE,fill=T,col="grey90")
#plot(coasts,add=TRUE,border="grey60",fill=TRUE,col="grey80",lwd=1)

# clustered based on gradientForest "turnover curves"
plot_adaptation(rfs,ras2[,bests],xy2,cluster.guide=pp,cluster.merge.threshold=0.001,color.by.medoids = T, nclust=7,pcs2show=2,flip1=1,flip2=-1,scal=8,cex=0.8,jitscale=0.05,rangeExp=1.5,coltxt="gray70",cluster.labels = T)
map(database="worldHires", add=TRUE,fill=T,col="grey90")
#plot(coasts, add=TRUE,col="grey80",fill=T,border="grey80",lwd=1)

# same as above with Wes Anderson-esque color scheme
plot_adaptation(rfs,ras2[,bests],xy2,cluster.guide=pp,cluster.merge.threshold=0,color.by.medoids = F, nclust=7,pcs2show=2,flip1=1,flip2=-1,scal=8,cex=0.8,jitscale=0.05,rangeExp=1.5,coltxt="gray70",cluster.labels = T)
map(database="worldHires", add=TRUE,fill=T,col="grey90")
#plot(coasts, add=TRUE,col="grey80",fill=T,border="grey80",lwd=1)
# dev.off()

# alternatively, we can use predictions averaged over ordinationJackknife replicates 
# (not always advantageous, try it either way)
rfs=sb$predictions

#pdf(file="agaricia_adaptationPlots_sb.pdf",height=4,width=8)
plot_adaptation(rfs,ras2[,bests],xy2,pcs2show=2,flip1=1,flip2=-1,scal=8,cex=0.8,jitscale=0.05,rangeExp=1.5,coltxt="gray70")
map(database="worldHires", add=TRUE,fill=T,col="grey90")
plot_adaptation(rfs,ras2[,bests],xy2,nclust=7,pcs2show=2,flip1=1,flip2=-1,scal=8,cex=0.8,color.by.medoids = T,jitscale=0.05,rangeExp=1.5,coltxt="gray70",cluster.labels = T)
map(database="worldHires", add=TRUE,fill=T,col="grey90")
#plot(coasts,add=TRUE,col="grey30",lwd=1)
plot_adaptation(rfs,ras2[,bests],xy2,cluster.guide=pp,cluster.merge.threshold=0.001,color.by.medoids = T, nclust=7,pcs2show=2,flip1=1,flip2=-1,scal=8,cex=0.8,jitscale=0.05,rangeExp=1.5,coltxt="gray70",cluster.labels = T)
map(database="worldHires", add=TRUE,fill=T,col="grey90")
#plot(coasts, add=TRUE,col="grey80",fill=T,border="grey80",lwd=1)
plot_adaptation(rfs,ras2[,bests],xy2,cluster.guide=pp,cluster.merge.threshold=0.001,color.by.medoids = F, nclust=7,pcs2show=2,flip1=1,flip2=-1,scal=8,cex=0.8,jitscale=0.05,rangeExp=1.5,coltxt="gray70",cluster.labels = T)
map(database="worldHires", add=TRUE,fill=T,col="grey90")
#dev.off()


