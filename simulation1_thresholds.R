# rm(list=ls())
# setwd("~/Dropbox/square_ecology_simulation")
library(RDAforest)
library(ggplot2)

# creating 50x50 environment
a=c(1:50)
b=c(1:50)
env=expand.grid(a,b)
names(env)=c("a","b")
# adding correlated variables
env$c=env$a+rnorm(nrow(env),0,10)
env$d=env$b+rnorm(nrow(env),0,10)
env$e=env$a+rnorm(nrow(env),0,25)
env$f=env$a+rnorm(nrow(env),0,15)
env$g=env$a+rnorm(nrow(env),0,15)
env$h=env$b+rnorm(nrow(env),0,30)
env$i=rnorm(nrow(env),0,40)
env$i[env$i<0]=0
env$j=rnorm(nrow(env),0,20)
env$k=rnorm(nrow(env),0,10)

# generating xy coordinates, rotated by angle theta with respect to envrionment
theta=0.8
rotation=diag(ncol(env))
rotation[1:2,1:2]=t(matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),ncol=2))
xy=t(rotation %*% t(env))[,1:2]

# plotting environment
rre=rda(env~1,scale=F)
nne=length(rre$CA$eig)
pce=scores(rre,scaling=1,display="sites",choices=c(1:nne))
pce0=decostand(pce,method="range")
env.col=decostand(env,method="range")
pdf(file="thresholds_envs.pdf",height=4.5,width=4.5)
plot(xy,pch=18,cex=0.75,col=rgb(pce0[,3],pce0[,2],pce0[,1],1),asp=1,main="environmental PCs",xaxt="na",yaxt="na",ylab="",xlab="",bty="none")
plot(xy,pch=18,cex=0.75,col=rgb(env.col$a,0.1,env.col$b,1),asp=1,main="a and b",xaxt="na",yaxt="na",ylab="",xlab="",bty="none")
plot(xy,pch=18,cex=0.75,col=rgb(env.col$a,0.1,0.1,1),asp=1,main="a",xaxt="na",yaxt="na",ylab="",xlab="",bty="none")
plot(xy,pch=18,cex=0.75,col=rgb(env.col$b,0.1,0.1,1),asp=1,main="b",xaxt="na",yaxt="na",ylab="",xlab="",bty="none")
plot(xy,pch=18,cex=0.75,col=rgb(env.col$c,0.1,0.1,1),asp=1,main="c",xaxt="na",yaxt="na",ylab="",xlab="",bty="none")
plot(xy,pch=18,cex=0.75,col=rgb(env.col$d,0.1,0.1,1),asp=1,main="d",xaxt="na",yaxt="na",ylab="",xlab="",bty="none")
plot(xy,pch=18,cex=0.75,col=rgb(env.col$e,0.1,0.1,1),asp=1,main="e",xaxt="na",yaxt="na",ylab="",xlab="",bty="none")
plot(xy,pch=18,cex=0.75,col=rgb(env.col$f,0.1,0.1,1),asp=1,main="f",xaxt="na",yaxt="na",ylab="",xlab="",bty="none")
plot(xy,pch=18,cex=0.75,col=rgb(env.col$g,0.1,0.1,1),asp=1,main="g",xaxt="na",yaxt="na",ylab="",xlab="",bty="none")
plot(xy,pch=18,cex=0.75,col=rgb(env.col$h,0.1,0.1,1),asp=1,main="h",xaxt="na",yaxt="na",ylab="",xlab="",bty="none")
plot(xy,pch=18,cex=0.75,col=rgb(env.col$i,0.1,0.1,1),asp=1,main="i",xaxt="na",yaxt="na",ylab="",xlab="",bty="none")
plot(xy,pch=18,cex=0.75,col=rgb(env.col$j,0.1,0.1,1),asp=1,main="j",xaxt="na",yaxt="na",ylab="",xlab="",bty="none")
plot(xy,pch=18,cex=0.75,col=rgb(env.col$k,0.1,0.1,1),asp=1,main="k",xaxt="na",yaxt="na",ylab="",xlab="",bty="none")
dev.off()

# clustering tree of predictors
pdf(file="thresholds_envCorrsTree.pdf",height=3.5,width=4)
hc=hclust(as.dist(1-cor(env)),method="ave")
plot(hc)
dev.off()

# setting "species"
sp1=c();sp2=c();sp3=c();sp3a=c();sp4=c();sp5=c();sp0=c()
for(r in 1:nrow(env)){
  i=env$a[r];j=env$b[r]
  sp0=c(sp0,sample(c(0,1),1)) # random
  if(i>=35 & j>=35) { sp1=c(sp1,1) } else { sp1=c(sp1,0) } #interaction 
  if(i>15 & i<35 & j>15 & j<35) { sp2=c(sp2,1) } else { sp2=c(sp2,0) } #interaction
  if(i>15) { sp3=c(sp3,1) } else { sp3=c(sp3,0) } # a-driven, more common
  if(i<25) { sp4=c(sp4,1) } else { sp4=c(sp4,0) } # a-driven, common
  if(j>15 && j<35) { sp5=c(sp5,1) } else { sp5=c(sp5,0) } # b-driven, less common
}

samp=data.frame(cbind(sp1,sp2,sp3,sp4,sp5))
rr=rda(samp~1)
nn=length(rr$CA$eig)
cpcs=scores(rr,scaling=1,display="sites",choices=c(1:nn))
cpcs0=decostand(cpcs,method="range")
pdf(file="thresholds_species.pdf",height=4.5,width=4.5)
plot(xy,pch=18,cex=1,col=rgb(cpcs0[,1],cpcs0[,3],cpcs0[,2],1),asp=1,main="original pcs of community",xaxt="na",yaxt="na",ylab="",xlab="",bty="none")
plot(xy,pch=18,cex=0.75,col=rgb(1-(0.2+samp[,1]/1.5),1-(0.2+samp[,1]/1.5),1-(0.2+samp[,1]/1.5),1),asp=1,main="sp1",xaxt="na",yaxt="na",ylab="",xlab="",bty="none")
plot(xy,pch=18,cex=0.75,col=rgb(1-(0.2+samp[,2]/1.5),1-(0.2+samp[,2]/1.5),1-(0.2+samp[,2]/1.5),1),asp=1,main="sp2",xaxt="na",yaxt="na",ylab="",xlab="",bty="none")
plot(xy,pch=18,cex=0.75,col=rgb(1-(0.2+samp[,3]/1.5),1-(0.2+samp[,3]/1.5),1-(0.2+samp[,3]/1.5),1),asp=1,main="sp3",xaxt="na",yaxt="na",ylab="",xlab="",bty="none")
plot(xy,pch=18,cex=0.75,col=rgb(1-(0.2+samp[,4]/1.5),1-(0.2+samp[,4]/1.5),1-(0.2+samp[,4]/1.5),1),asp=1,main="sp4",xaxt="na",yaxt="na",ylab="",xlab="",bty="none")
plot(xy,pch=18,cex=0.75,col=rgb(1-(0.2+samp[,5]/1.5),1-(0.2+samp[,5]/1.5),1-(0.2+samp[,5]/1.5),1),asp=1,main="sp5",xaxt="na",yaxt="na",ylab="",xlab="",bty="none")
dev.off()

# choosing samples to fit the model (1/10 of total), plotting them
picks=sample(1:nrow(samp),round(nrow(samp)/10))
samp.fit=samp[picks,]
env.fit=env[picks,]
XY=data.frame(xy)
names(XY)=c("x","y")
XY.fit=XY[picks,]
cpcs0.fit=cpcs0[picks,]
pdf(file="thresholds_samples4model.pdf",height=4.5,width=4.5)
plot(XY.fit,pch=18,cex=0.75,col=rgb(cpcs0.fit[,1],cpcs0.fit[,3],cpcs0.fit[,2],1),asp=1,main="samples for fitting",xaxt="na",yaxt="na",ylab="",xlab="",bty="none")
dev.off()

# generating distance matrix for fitting, also scores (for direct RF prediction)
dis=data.frame(as.matrix(vegdist(samp.fit,method="jaccard")))
caps=capscale(dis~1)
nc=length(caps$CA$eig)
sc.fit=data.frame(scores(caps,scaling=1,display="sites",choices=c(1:nc)))

# ------ Gradient Forest: exploring turnover curves

# fitting gradient forest model to ordination
gf=makeGF(caps,env.fit,ntree=1500)
# which PCs are predicted wiht non-zero R2, and how well?
gf$result

# plotting bar chart of importances (scaled to proportion of total variance explained)
imps=data.frame(importance_RDAforest(gf,caps))
# total R2
sum(imps)
names(imps)="R2"
imps$var=row.names(imps)
imps$var=factor(imps$var,levels=imps$var[order(imps$R2)])
pdf(file="thresholds_GFimportances.pdf",width=2,height=2.5)
ggplot(imps,aes(var,R2))+geom_bar(stat="identity")+coord_flip()+theme_bw()
dev.off()

# plot turnover curves (first 2 most important environmental variables)
vars=imps$var[1:2]
pdf(file="thresholds_GF_turnovers.pdf",width=4.4,height=4.4)
plot_gf_turnovers(gf,vars)
dev.off()

#--------- variable selection

mm=mtrySelJack(dis,env.fit,nreps=20,top.pcs = nc)
# selected predictors
mm$goodvars
# median importances at higher mtry (proportion of total variation explained)
mm$importances

library(ggplot2)
# plot mtrySelJack diagnostic plots
pdf(file="thresholds_mtrySelJack.pdf",width=2,height=2.5)
ggplot(mm$prop.positive,aes(var,prop.positive))+
  geom_bar(stat="identity")+
  coord_flip()+
  theme_bw()+
  geom_hline(yintercept=0.5,col="red")
dev.off()

#------- measuring importance of passing variables, forming predictions

sb=ordinationJackknife(Y=dis,X=env.fit[,c("a","b")],newX=env,nreps=10,top.pcs=nc)

# plotting boxplot of importances (across jackknifing replica)
pdf(file="thresholds_sb.pdf",height=1.5,width=3)
important=names(sb$median.importance)
ggplot(sb$all.importances[sb$all.importances$variable %in% important,],aes(variable,importance))+
  geom_boxplot()+
  coord_flip()+
  theme_bw()
dev.off()

#-------------- plotting community maps / adaptive neighborhoods

dirs=sb$predictions.direct
turnovers=sb$predictions.turnover
nbests=2
bests=names(sb$median.importance)[c(1:nbests)]

# subsetting for spots within modeled range
table(sb$goodrows)
xy2=xy[sb$goodrows,]
ras2=env[which(sb$goodrows),]

pdf("thresholds_adaptation.pdf",width=5,heigh=4)
#  RF predictions
pa0=plot_adaptation(dirs,ras2[,bests],xy2,matching.scores=cpcs0[sb$goodrows,],lighten=0,color.scheme="000",ordinate=T,scal=8,cex=0.8,jitscale=0.05,rangeExp=1.5,col.envs="red",main="unclustered RF predictions")
# turnover predictions
pa2=plot_adaptation(turnovers,ras2[,bests],xy2,nclust=12,cluster.merge.level=0.1,lighten=0,color.scheme="000",scal=8,cex=0.8,jitscale=0.05,cluster.labels=F,rangeExp=1.5,col.envs="red",main="clustered by GF, nclust=12")
# RF predictions directly clustered
pa1=plot_adaptation(dirs,ras2[,bests],xy2,matching.scores=cpcs0[sb$goodrows,],lighten=0,nclust=12,color.scheme="000",ordinate=T,cluster.merge.level=0.1,scal=8,cex=0.8,jitscale=0.05,cluster.labels=F,rangeExp=1.5,col.envs="red",main="clustered by RF, nclust=12")
# turnove predictions clustered, clusters merged based on RF
pa3=plot_adaptation(dirs,ras2[,bests],xy2,nclust=12,matching.scores=cpcs0[sb$goodrows,],lighten=0,ordinate=T,color.scheme="000",cluster.merge.level=0.1,cluster.guide = turnovers,cluster.labels=F,col.envs="red",scal=8,cex=0.8,jitscale=0.05,rangeExp=1.5,main="clustered by turnover\nmerged by RF, nclust=12")
dev.off()

#-----predictions based on just spatial coordinates

# choosing 1/5 of total samples to fit 
picks=sample(1:nrow(samp),round(nrow(samp)/5))
samp.fit=samp[picks,]
env.fit=env[picks,]
XY=data.frame(xy)
names(XY)=c("x","y")
XY.fit=XY[picks,]
cpcs0.fit=cpcs0[picks,]

# generating distance matrix for fitting, also scores (for direct RF prediction)
dis=data.frame(as.matrix(vegdist(samp.fit,method="jaccard")))
#caps=capscale(dis~1+Condition(as.matrix(XY.fit)))
caps=capscale(dis~1)
nc=length(caps$CA$eig)
sc.fit=data.frame(scores(caps,scaling=1,display="sites",choices=c(1:nc)))

# predictions based only on spatial coordinates
sbxy=ordinationJackknife(Y=dis,X=XY.fit,newX=XY,nreps=5,top.pcs=nc)
dirs=sbxy$predictions.direct
turnovers=sbxy$predictions.turnover
nbests=2
bests=names(sbxy$median.importance)[c(1:nbests)]

# subsetting for spots within modeled range
table(sbxy$goodrows)
xy2=xy[sbxy$goodrows,]
ras2=XY[which(sbxy$goodrows),]

pdf("thresholds_adaptation_XY_2xSamples_v4.pdf",width=5,heigh=4)
# unclustered
pa0=plot_adaptation(dirs,ras2[,bests],xy2,matching.scores=cpcs0[sbxy$goodrows,],lighten=0,ordinate=F,scal=8,cex=0.8,jitscale=0.05,rangeExp=1.5,col.envs="red",main="unclustered RF predictions")
# clustered directly by RF
pa1=plot_adaptation(dirs,ras2[,bests],xy2,nclust=12,matching.scores=cpcs0[sbxy$goodrows,],cluster.merge.level=0,scal=8,cex=0.8,jitscale=0.05,cluster.labels=F,rangeExp=1.5,col.envs="red",main="clustered by RF, nclust=12")
dev.off()
