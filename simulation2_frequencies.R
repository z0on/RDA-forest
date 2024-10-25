#rm(list=ls())
# setwd("~/Dropbox/square_ecology_simulation")
library(RDAforest)
library(ggplot2)

#setting up environment
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

pdf(file="frequencies_envCorrsTree.pdf",height=3.5,width=3.5)
hc=hclust(as.dist(1-cor(env)),method="ave")
plot(hc)
dev.off()

# rotated x, y coordinates
th=0.8
rot=diag(ncol(env))
rot[1:2,1:2]=t(matrix(c(cos(th),-sin(th),sin(th),cos(th)),ncol=2))
xy=t(rot %*% t(env))[,1:2]
env.col=decostand(env,method="range")

# plotting environments
pdf(file="frequencies_envs.pdf",height=4.5,width=4.5)
plot(xy,pch=18,cex=0.75,col=rgb(env.col$a,0.1,env.col$b,1),asp=1,main="a + b",xaxt="na",yaxt="na",ylab="",xlab="",bty="none")
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

# trends for species distributions
a.trend=(a/max(a))^2
a.trend[a.trend<0.01]=0.01
a.trend[a.trend>0.99]=0.99
plot(a.trend)
b.trend=(b/max(b))^2
b.trend[b.trend<0.01]=0.01
b.trend[b.trend>0.99]=0.99
plot(b.trend)
a.peak=(1-decostand(abs(scale(a)[,1]),method="range")[,1])
b.peak=(1-decostand(abs(scale(b)[,1]),method="range")[,1])
a.peak[a.peak<0.1]=0.1
a.peak[a.peak>0.9]=0.9
b.peak[b.peak<0.1]=0.1
b.peak[b.peak>0.9]=0.9
plot(a.peak)
plot((a.trend*b.trend))

# setting up species
sp1=c();sp2=c();sp3=c();sp2a=c();sp4=c();sp5=c();sp0=sp01=sp02=sp03=sp04=sp5a=sp5b=c()
for(r in 1:nrow(env)){
    i=env[r,"a"];j=env[r,"b"]
    sp0=c(sp0,1) # everywhere (so jaccard distance would work)
    sp01=c(sp01,sample(c(0,1),1)) # random
    sp02=c(sp02,sample(c(0,1),1)) # random
    sp03=c(sp03,sample(c(0,1),1)) # random
    sp04=c(sp04,sample(c(0,1),1)) # random
    sp1=c(sp1,sample(c(0,1),1,prob=c(1-a.trend[i],a.trend[i]))) # a-trend
    sp2=c(sp2,sample(c(0,1),1,prob=c(1-b.trend[j],b.trend[j]))) # b-trend
    sp2a=c(sp2a,sample(c(0,1),1,prob=c(1-b.trend[j]^2,b.trend[j]^2))) # b-trend
    apeakprob=a.peak[i]^2
    sp3=c(sp3,sample(c(0,1),1,prob=c(1-apeakprob,apeakprob))) # a-peak-driven
    abpeakprob=a.peak[i]*b.peak[j]
    if(abpeakprob>0.5) { abprob=0.9}
    if(abpeakprob<0.4) { abprob=0.1}
    sp4=c(sp4,sample(c(0,1),1,prob=c(1-(abpeakprob),(abpeakprob)))) # interaction center
    abprob=(a.trend[i]*b.trend[j])
#    if(i<45 | j<45) { abrpob=0.001 }
    if(abprob>0.6) { abprob=0.9}
    if(abprob<0.4) { abprob=0.1}
    sp5=c(sp5,sample(c(0,1),1,prob=c(1-abprob,abprob))) # interaction top corner
    sp5a=c(sp5a,sample(c(0,1),1,prob=c(1-abprob,abprob))) # interaction top corner
    sp5b=c(sp5b,sample(c(0,1),1,prob=c(1-abprob,abprob))) # interaction top corner
}

samp=data.frame(cbind(sp0,sp2,sp3,sp4,sp5,sp5a,sp5b))
names(samp)
dis0=data.frame(as.matrix((vegdist(samp,method="jaccard"))))
caps0=capscale(dis0~1)
nn=length(caps0$CA$eig)
cpcs=scores(caps0,scaling=1,display="sites",choices=c(1:nn))
cpcs0=decostand(cpcs,method="range")
pdf(file="frequencies_species.pdf",height=4.5,width=4.5)
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
pdf(file="frequencies_samples4model.pdf",height=4.5,width=4.5)
plot(XY.fit,pch=18,cex=1,col=rgb(cpcs0.fit[,1],cpcs0.fit[,3],cpcs0.fit[,2],1),asp=1,main="original pcs of community",xaxt="na",yaxt="na",ylab="",xlab="",bty="none")
dev.off()

# generating distance matrix for fitting, also scores (for direct RF prediction)
dis=data.frame(as.matrix((vegdist(samp.fit,method="jaccard"))))
caps=capscale(dis~1)
nc=length(caps$CA$eig)
sc.fit=data.frame(scores(caps,scaling=1,display="sites",choices=c(1:nc)))

# ------ Gradient Forest: exploring turnover curves

# fitting gradient forest model to ordination
gf=makeGF(caps,env.fit,ntree=1500)
gf$result
nc=length(gf$result)

# plot importances (R2)
imps=data.frame(importance_RDAforest(gf,caps))
# total R2
sum(imps)
names(imps)="R2"
imps$var=row.names(imps)
imps$var=factor(imps$var,levels=imps$var[order(imps$R2)])
pdf(file="frequencies_GFimportances.pdf",width=2,height=2.5)
ggplot(imps,aes(var,R2))+geom_bar(stat="identity")+coord_flip()+theme_bw()
dev.off()

# plot turnover curves (first 4 most important environmental variables)
vars=names(importance(gf))[1:2]
pdf(file="frequencies_GF_turnovers.pdf",width=4.4,height=4.4)
plot_gf_turnovers(gf,vars)
dev.off()

#--------- variable selection

mm=mtrySelJack(dis,env.fit,nreps=20,top.pcs = length(gf$result),importance.type = "GF")
mm$goodvars

# plot mtrySelJack diagnostic plots
pdf(file="frequencies_mtrySelJack.pdf",,width=2,height=2.5)
# bar chart of proportion of positive change in response to higher mtry
# good variables would be the ones above the red line
ggplot(mm$prop.positive,aes(var,prop.positive))+
  geom_bar(stat="identity")+
  coord_flip()+
  theme_bw()+
  geom_hline(yintercept=0.5,col="red")
dev.off()

#------- measuring importance of passing variables, forming predictions

sb=ordinationJackknife(Y=dis,X=env.fit[,mm$goodvars],newX=env,nreps=10,top.pcs=nc)

pdf(file="frequencies_sb.pdf",height=1.5,width=3)
important=names(sb$median.importance)
ggplot(sb$all.importances[sb$all.importances$variable %in% important,],aes(variable,importance))+
  geom_boxplot()+
  coord_flip()+
  #  ylim(0,0.06)+
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
?plot_adaptation
pdf("frequencies_adaptation.pdf",width=5,heigh=4)
# unclustered
pa0=plot_adaptation(dirs,ras2[,bests],xy2,lighten=0,color.scheme="000",ordinate=T,scal=8,cex=0.8,jitscale=0.05,rangeExp=1.5,col.envs="red",main="unclustered RF predictions")
# RF clusters
pa1=plot_adaptation(dirs,ras2[,bests],xy2,lighten=0,nclust=12,color.scheme="000",ordinate=T,cluster.merge.level=0.25,scal=8,cex=0.8,jitscale=0.05,cluster.labels=F,rangeExp=1.5,col.envs="red",main="clustered by RF, nclust=12")
# turnover clusters
pa2=plot_adaptation(turnovers,ras2[,bests],xy2,nclust=12,cluster.merge.level=0.1,lighten=0,color.scheme="000",scal=8,cex=0.8,jitscale=0.05,cluster.labels=F,rangeExp=1.5,col.envs="red",main="clustered by GF, nclust=12")
# turnover clusters, merged by RF
pa3=plot_adaptation(dirs,ras2[,bests],xy2,nclust=12,lighten=0,ordinate=T,color.scheme="000",cluster.merge.level=0.25,cluster.guide = turnovers,cluster.labels=F,col.envs="red",scal=8,cex=0.8,jitscale=0.05,rangeExp=1.5,main="turnover clusters\nmerged by RF, nclust=12")
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
sbxy=ordinationJackknife(Y=dis,X=XY.fit,newX=XY,nreps=7,top.pcs=nc)

dirs=sbxy$predictions.direct
turnovers=sbxy$predictions.turnover
nbests=2
bests=names(sbxy$median.importance)[c(1:nbests)]

# subsetting for spots within modeled range
table(sbxy$goodrows)
xy2=xy[sbxy$goodrows,]
ras2=XY[which(sbxy$goodrows),]

pdf("frequencies_adaptation_XY_2xSamples.pdf",width=5,heigh=4)
# unclustered
pa0=plot_adaptation(dirs,ras2[,bests],xy2,color.scheme="000",scal=8,cex=0.8,jitscale=0.05,rangeExp=1.5,col.envs="red",main="unclustered RF predictions")
dev.off()

