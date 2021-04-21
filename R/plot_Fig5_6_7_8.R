# Script to draw explanatory figures for benchmark building 

# This script can reproduce Fig. 5 to 8 in the research paper:
# 'Using the International Tree-Ring Data Bank (ITRDB) records as century-long benchmarks 
# for land-surface models'

# Author: Jina Jeong (j.jeong@vu.nl)

#- 0. Initialization ####
library(ncdf4)
library(dplR)

source('R/functions_J.R')

load('output/ITRDB_simulation.RData')
simulation=simulation.cny.power
load('data/ITRDB_obs.RData')


#- 1. Growth Trend ####
#-- 1.1 Process data ####
site = 'finl052'
trw.sim = simulation[[site]]
trw.obs = obs[[site]]
trw.age = apply(trw.obs,2,align.trw)
nyears = nrow(trw.obs)

#-- 1.2 Plot ####
png(paste0('figure/trend_',site,'.png'), width=11, height=9, 
    units="in", res=300)
par(mar=c(4,4.5,4.1,2.1))
m <- matrix(c(1,2,3,4,5,5),nrow = 3,ncol = 2,byrow = T)
layout(mat = m,heights = c(0.4,0.4,0.2))

# Trend: Step1
plot(-100,-100,ylim = c(0,1.1*max(trw.obs,trw.sim,na.rm=T)),xlim=c(0,nyears),
     xlab = '',cex.axis = 1.8,ylab='',type='l',col='azure3',xaxt = "n")
apply(trw.age,2, function(x) lines(x,cex=0.7,col='azure3'))

interval = 40
at = seq(1,nyears,by=interval)
axis(1, at=seq(1,nyears,by=interval), labels = at,cex.axis=1.8)
mtext(side=2,'Tree-ring width (mm)',line=2.7,cex=1.3)
mtext(side=1,'Tree age (years)',line=2.7,cex=1.3)

lines(apply(trw.age,1,mean,na.rm=T),lty=2,lwd=2)
apply(trw.sim,2,lines,col='dodgerblue2')
lines(trw.sim[,5],col='dodgerblue2',lwd=2)

# Trend: Step2
plot(-100,-100,ylim = c(0,max(trw.obs,trw.sim,na.rm=T)),xlim=c(0,nyears),
     xlab = '',cex.axis = 1.8,ylab='',type='l',col='azure3',xaxt = "n")

interval = 40
at = seq(1,nyears,by=interval)
axis(1, at=seq(1,nyears,by=interval), labels = at,cex.axis=1.8)
mtext(side=2,'Tree-ring width (mm)',line=2.7,cex=1.3)
mtext(side=1,'Tree age (years)',line=2.7,cex=1.3)

lines(apply(trw.age,1,tbrm),lty=2,lwd=2)
lines(trw.sim[,5],col='dodgerblue2',lwd=2)

# Trend: Step3
sub = apply(trw.age,1,tbrm)-trw.sim[,5]
reg = lm(apply(trw.age,1,tbrm)-trw.sim[,5]~c(1:nyears))

plot(-100,-100,ylim = c(min(sub),max(sub)),xlim=c(0,nyears),
     xlab = '',cex.axis = 1.8,ylab='',type='l',col='azure3',xaxt = "n")

interval = 40
at = seq(1,nyears,by=interval)
axis(1, at=seq(1,nyears,by=interval), labels = at,cex.axis=1.8)
mtext(side=2,'Observation - simulation (mm)',line=2.7,cex=1.3)
mtext(side=1,'Tree age (years)',line=2.7,cex=1.3)

lines(sub,type='l',col='green3',lwd=1.5)

# Trend: Step4
sub = apply(trw.age,1,tbrm)-trw.sim[,5]
reg = lm(apply(trw.age,1,tbrm)-trw.sim[,5]~c(1:nyears))

plot(-100,-100,ylim = c(min(sub),max(sub)),xlim=c(0,nyears),
     xlab = '',cex.axis = 1.8,ylab='',type='l',col='azure3',xaxt = "n")

interval = 40
at = seq(1,nyears,by=interval)
axis(1, at=seq(1,nyears,by=interval), labels = at,cex.axis=1.8)
mtext(side=2,'Observation - simulation (mm)',line=2.7,cex=1.3)
mtext(side=1,'Tree age (years)',line=2.7,cex=1.3)

lines(sub,type='l',col='green3',lwd=1.5)
abline(reg,col='green4',lwd=1.7,lty=2)

# legend 
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend('top',c('Observation','Simulated output','Compiled virtual tree','Residuals',
               'Linear regression of the residuals'),
       col=c('azure3','dodgerblue','black','green3','green4'),
       lty=c(1,1,2,1,2),lwd=c(1,1,2,1.5,1.7),ncol=3,cex=1.3)


dev.off()


#- 2. Diamater growth of matrue trees ####
#- 2.1 Process data ####
site = 'brit021'
trw.sim = simulation[[site]]
trw.obs = obs[[site]]
nyears = nrow(trw.obs)

dia.end.obs = apply(trw.obs,2,cumsum.j)
dia.sim = apply(trw.sim,2,cumsum.j)

#-- 2.2 Plot ####

png(paste0('figure/mature_',site,'.png'), width=11, height=12.5, units="in", res=300)
par(mar=c(5.1,5,3.5,2.1))
m <- matrix(c(1,2,3,4,5,6,7,7),nrow = 4,ncol = 2,byrow = T)
layout(mat = m,heights = c(0.4,0.4,0.4,0.23))


# Mature: Step1
plot(-100,-100,ylim = c(0,max(dia.end.obs,dia.sim,na.rm=T)),xlim=c(0,nyears),
     xlab = '',cex.axis = 1.8,ylab='',type='l',col='azure3',xaxt = "n")
apply(dia.end.obs,2, function(x) lines(x,cex=0.7,col='azure3'))

interval = 40
at = seq(1,nyears,by=interval)
axis(1, at=seq(1,nyears,by=interval), labels = row.names(trw.obs)[at],cex.axis=1.8)
mtext(side=2,'Diameter (mm)',line=3.1,cex=1.5)
mtext(side=1,'Calendar year',line=3.5,cex=1.5)
lines(apply(dia.end.obs,1,mean,na.rm=T),lwd=2.3,lty=2)

apply(dia.sim,2,lines,col='dodgerblue2')
lines(dia.sim[,5],col='dodgerblue2',lwd=2)

# Mature: Step2
plot(-100,-100,ylim = c(0,max(dia.end.obs,dia.sim,na.rm=T)),xlim=c(0,nyears),
     xlab = '',cex.axis = 1.8,ylab='',type='l',col='azure3',xaxt = "n")

interval = 40
at = seq(1,nyears,by=interval)
axis(1, at=seq(1,nyears,by=interval), labels = row.names(trw.obs)[at],cex.axis=1.8)
mtext(side=2,'Diameter (mm)',line=3.1,cex=1.5)
mtext(side=1,'Calendar year',line=3.5,cex=1.5)

lines(apply(dia.end.obs,1,mean,na.rm=T),lwd=2.3,lty=2)
lines(dia.sim[,5],col='dodgerblue2',lwd=2)

# Mature: Step3
plot(-100,-100,ylim = c(0,max(dia.end.obs,dia.sim,na.rm=T)),xlim=c(0,nyears),
     xlab = '',cex.axis = 1.8,ylab='',type='l',col='azure3',xaxt = "n")

interval = 40
at = seq(1,nyears,by=interval)
axis(1, at=seq(1,nyears,by=interval), labels = row.names(trw.obs)[at],cex.axis=1.8)
mtext(side=2,'Diameter (mm)',line=3.1,cex=1.5)
mtext(side=1,'Calendar year',line=3.5,cex=1.5)

lines(apply(dia.end.obs,1,mean,na.rm=T),lwd=2.3,lty=2)
lines(dia.sim[,5],col='dodgerblue2',lwd=2)

x = seq(1,nyears,by=6)
arrows(x,apply(dia.end.obs,1,mean,na.rm=T)[x],x,dia.sim[x,5],length=0.08,angle=30,
       code=3)

# Mature: Step4
plot(-100,-100,ylim = c(0,max(dia.end.obs,dia.sim,na.rm=T)),xlim=c(0,nyears),
     xlab = '',cex.axis = 1.8,ylab='',type='l',col='azure3',xaxt = "n")

interval = 40
at = seq(1,nyears,by=interval)
axis(1, at=seq(1,nyears,by=interval), labels = row.names(trw.obs)[at],cex.axis=1.8)
mtext(side=2,'Diameter (mm)',line=3.1,cex=1.5)
mtext(side=1,'Calendar year',line=3.5,cex=1.5)

lines(apply(dia.end.obs,1,mean,na.rm=T),lwd=2.3,lty=2)
lines(apply(dia.end.obs,1,mean,na.rm=T)[1:50],lwd=2.5,col='white')
lines(dia.sim[,5],col='dodgerblue3',lwd=2)
lines(dia.sim[1:50,5],col='white',lwd=2.5)

x = seq(52,nyears,by=6)
arrows(x,apply(dia.end.obs,1,mean,na.rm=T)[x],x,dia.sim[x,5],length=0.08,angle=30,
       code=3)

# Mature: Step5
sub = apply(dia.end.obs,1,mean,na.rm=T)[51:nyears]-dia.sim[51:nyears,5]
reg = lm(sub~c(51:nyears))

plot(-100,-100,ylim = c(min(sub),max(sub)),xlim=c(51,nyears),
     xlab = '',cex.axis = 1.8,ylab='',type='l',col='azure3',xaxt = "n")

interval = 20
at = seq(53,nyears,by=interval)
axis(1, at=seq(53,nyears,by=interval), labels = row.names(trw.obs)[at],cex.axis=1.8)
mtext(side=2,'Observation - simulation (mm)',line=3.1,cex=1.5)
mtext(side=1,'Calendar year',line=3.5,cex=1.5)

lines(c(51:nyears),sub,type='l',col='green3',lwd=1.5)

# Mature: Step6
plot(-100,-100,ylim = c(min(sub),max(sub)),xlim=c(51,nyears),
     xlab = '',cex.axis = 1.8,ylab='',type='l',col='azure3',xaxt = "n")

interval = 20
at = seq(53,nyears,by=interval)
axis(1, at=seq(53,nyears,by=interval), labels = row.names(trw.obs)[at],cex.axis=1.8)
mtext(side=2,'Observation - simulation (mm)',line=3.1,cex=1.5)
mtext(side=1,'Calendar year',line=3.5,cex=1.5)

lines(c(51:nyears),sub,type='l',col='green3',lwd=1.5)
abline(reg,col='green4',lwd=1.7,lty=2)

# legend 
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend('top',c('Observation','Simulated output','Compiled virtual tree','Residuals      ',
               'Linear regression of residuals',
               'Difference between simulated \n big tree and virtual tree'),
       col=c('azure3','dodgerblue','black','green3','green4',NA),
       lty=c(1,1,2,1,2,NA),lwd=c(1,1,2,1.5,1.7,NA),ncol=3,cex=1.3,
       text.width = c(0.17,0.17,0.17),
       x.intersp= 1,y.intersp=1.2,adj=0)
dev.off()


#- 3 Diamater growth of young trees ####
#-- 3.1 Process data ####
site = 'brit021'
trw.sim = simulation[[site]]
trw.obs = obs[[site]]
nyears = nrow(trw.obs)

trw.age = apply(trw.obs,2,align.trw)
dia.age.obs = apply(trw.age,2,cumsum.j)
dia.sim = apply(trw.sim,2,cumsum.j)

#-- 3.2 Young: Plot #### 
png(paste0('figure/young_',site,'.png'), width=11, height=12.5, units="in", res=300)
par(mar=c(5.1,5,3.5,2.1))
m <- matrix(c(1,2,3,4,5,6,7,7),nrow = 4,ncol = 2,byrow = T)
layout(mat = m,heights = c(0.4,0.4,0.4,0.23))


# Young: Step1
plot(-100,-100,ylim = c(0,max(dia.end.obs,dia.sim,na.rm=T)),xlim=c(0,nyears),
     xlab = '',cex.axis = 1.8,ylab='',type='l',col='azure3',xaxt = "n")

interval = 40
at = seq(1,nyears,by=interval)
axis(1, at=seq(1,nyears,by=interval), labels = at,cex.axis=1.8)
mtext(side=2,'Diameter (mm)',line=3.1,cex=1.5)
mtext(side=1,'Tree age (years)',line=3.5,cex=1.5)

apply(dia.age.obs,2,lines,cex=0.7,col='azure3')
lines(apply(dia.age.obs,1,max,na.rm=T),lwd=2.3,lty=2)
apply(dia.sim,2,lines,col='dodgerblue2')
lines(dia.sim[,5],col='dodgerblue2',lwd=2)

# Young: Step2
plot(-100,-100,ylim = c(0,max(dia.end.obs,dia.sim,na.rm=T)),xlim=c(0,nyears),
     xlab = '',cex.axis = 1.8,ylab='',type='l',col='azure3',xaxt = "n")

interval = 40
at = seq(1,nyears,by=interval)
axis(1, at=seq(1,nyears,by=interval), labels = at,cex.axis=1.8)
mtext(side=2,'Diameter (mm)',line=3.1,cex=1.5)
mtext(side=1,'Tree age (years)',line=3.5,cex=1.5)

lines(apply(dia.age.obs,1,max,na.rm=T),lwd=2.3,lty=2)
lines(dia.sim[,5],col='dodgerblue2',lwd=2)

# Young: Step3
plot(-100,-100,ylim = c(0,max(dia.end.obs,dia.sim,na.rm=T)),xlim=c(0,nyears),
     xlab = '',cex.axis = 1.8,ylab='',type='l',col='azure3',xaxt = "n")

year.index = c(31:nyears)
interval = 40
at = seq(1,nyears,by=interval)
axis(1, at=seq(1,nyears,by=interval), labels = at,cex.axis=1.8)
mtext(side=2,'Diameter (mm)',line=3.1,cex=1.5)
mtext(side=1,'Tree age (years)',line=3.5,cex=1.5)

lines(apply(dia.age.obs,1,max,na.rm=T),lwd=2.3,lty=2)
lines(dia.sim[,5],col='dodgerblue2',lwd=2)

x = seq(1,nyears,by=6)
arrows(x,apply(dia.age.obs,1,max,na.rm=T)[x],x,dia.sim[x,5],length=0.08,angle=30,
       code=3)

# Young: Step4
plot(-100,-100,ylim = c(0,max(dia.end.obs,dia.sim,na.rm=T)),xlim=c(0,nyears),
     xlab = '',cex.lab=1.5,cex.axis = 1.8,ylab='',type='l',col='azure3',xaxt = "n")

year.index = c(31:nyears)
interval = 40
at = seq(1,nyears,by=interval)
axis(1, at=seq(1,nyears,by=interval), labels = at,cex.axis=1.8)
mtext(side=2,'Diameter (mm)',line=3.1,cex=1.5)
mtext(side=1,'Tree age (years)',line=3.5,cex=1.5)

lines(apply(dia.age.obs,1,max,na.rm=T),lwd=2.3,lty=2)
lines(year.index,apply(dia.age.obs,1,max,na.rm=T)[year.index],lwd=2.5,col='white')
lines(dia.sim[,5],col='dodgerblue2',lwd=2)
lines(year.index,dia.sim[year.index,5],col='white',lwd=2.5)

x = seq(1,30,by=6)
arrows(x,apply(dia.age.obs,1,max,na.rm=T)[x],x,dia.sim[x,5],length=0.08,angle=30,
       code=3)

# Young: Step5
sub = apply(dia.age.obs,1,max,na.rm=T)[1:30]-dia.sim[1:30,5]
reg = lm(sub~c(1:30))

plot(-100,-100,ylim = c(min(sub),max(sub)),xlim=c(0,30),
     xlab = '',cex.axis = 1.8,ylab='',type='l',col='azure3',xaxt = "n")

interval = 5
at = seq(5,nyears,by=interval)
axis(1, at=seq(5,nyears,by=interval), labels = at,cex.axis=1.4)
mtext(side=2,'Observation - simulation (mm)',line=3.1,cex=1.5)
mtext(side=1,'Tree age (years)',line=3.5,cex=1.5)

lines(sub,type='l',col='green3',lwd=1.5)

# Young: Step6
plot(-100,-100,ylim = c(min(sub),max(sub)),xlim=c(0,30),
     xlab = '',cex.axis = 1.8,ylab='',type='l',col='azure3',xaxt = "n")

interval = 5
at = seq(5,nyears,by=interval)
axis(1, at=seq(5,nyears,by=interval), labels = at,cex.axis=1.8)
mtext(side=2,'Observation - simulation (mm)',line=3.1,cex=1.5)
mtext(side=1,'Tree age (years)',line=3.5,cex=1.5)

lines(sub,type='l',col='green3',lwd=1.5)
abline(reg,col='green4',lwd=1.7,lty=2)

# legend 
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend('top',c('Observation','Simulated output','Compiled virtual tree','Residuals',
               'Linear regression of residuals',
               'Difference between simulated \n big tree and virtual tree'),
       col=c('azure3','dodgerblue','black','green3','green4',NA),
       lty=c(1,1,2,1,2,NA),lwd=c(1,1,2,1.5,1.7,NA),ncol=3,cex=1.3,
       text.width = c(0.17,0.17,0.17),
       x.intersp= 1,y.intersp=1.2,adj=0)
dev.off()


#- 4. Extreme growths ####
#-- 4.1 Process data ####
site = 'spai006'
trw.obs = obs[[site]]
trw.sim = simulation[[site]]
p=0.25
target_circ=5

idx.age = !is.na(trw.obs[rownames(trw.obs)==1950,]) & apply(trw.obs,2,function(x) sum(!is.na(x))) >50
obs.long = trw.obs[,idx.age]
obs.cut=obs.long[which(as.numeric(rownames(obs.long))>1950),]
obs.cut=rowMeans(obs.cut,na.rm=T)

sim.cut=trw.sim[which(as.numeric(rownames(obs.long))>1950),target_circ]
names(sim.cut)<-names(obs.cut)

ind.obs.low = order(obs.cut)[1:floor(p*length(obs.cut))]
ind.obs.high= rev(order(obs.cut,decreasing=T)[1:floor(p*length(obs.cut))])
ind.sim.low = order(sim.cut)[1:floor(p*length(sim.cut))]
ind.sim.high= rev(order(sim.cut,decreasing=T)[1:floor(p*length(sim.cut))])

obs.extr = scale(obs.cut)[c(ind.obs.low,ind.obs.high),1]
sim.extr = scale(sim.cut)[c(ind.obs.low,ind.obs.high),1]

obs.rank = obs.cut[c(ind.obs.low,ind.obs.high)]
sim.rank = sim.cut[c(ind.sim.low,ind.sim.high)]
obs.rank = obs.rank-mean(obs.rank)
sim.rank = sim.rank-mean(sim.rank)
``

#-- 4.2 Exterem: Plot ####
png(paste0('figure/extreme_',site,'.png'), width=11.5, height=12.5, units="in", res=300)
par(mar=c(5.1,6.6,4.8,1.3))
m <- matrix(c(1,2,3,4,5,6,7,7),nrow = 4,ncol = 2,byrow = T)
layout(mat = m,heights = c(0.4,0.4,0.4,0.28))

# Extreme: Step1
plot(names(obs.cut),obs.cut,type='l',ylim=extendrange(obs.cut,f=0.15),
     ylab='',xlab='',xaxt='n',cex.axis=1.8)

interval = 5
at = seq(1,length(obs.cut),by=interval)
axis(1, at=names(obs.cut)[at], labels = names(obs.cut)[at],cex.axis=1.8)
mtext(side=2,'Annual-mean of long-lived \n observed trees (mm)',line=3,cex=1.3)
mtext(side=1,'Calendar year',line=3.5,cex=1.3)

rect(1900,min(obs.cut[ind.obs.high]),2100,max(obs.cut[ind.obs.high]),
     density=10,border=NA,col='indianred1')
rect(1900,min(obs.cut[ind.obs.low]),2100,max(obs.cut[ind.obs.low]),
     density=10,border=NA,col='lightskyblue')
points(names(obs.cut[ind.obs.high]),obs.cut[ind.obs.high],pch=19,cex=0.7)
points(names(obs.cut[ind.obs.low]),obs.cut[ind.obs.low],pch=19,cex=0.7)

axis(side=1,at=names(obs.cut[ind.obs.high]),tcl=0.4,lwd.ticks=3,mgp=c(0,0.5,0),label=NA,
     col.ticks = 'indianred1')
segments(as.numeric(names(obs.cut[ind.obs.high])),-1,as.numeric(names(obs.cut[ind.obs.high])),
         obs.cut[ind.obs.high],col='indianred1',lwd=0.5,lty=2)

axis(side=1,at=names(obs.cut[ind.obs.low]),tcl=0.4,lwd.ticks=3,mgp=c(0,0.5,0),label=NA,
     col.ticks = 'lightskyblue')
segments(as.numeric(names(obs.cut[ind.obs.low])),-1,as.numeric(names(obs.cut[ind.obs.low])),
         obs.cut[ind.obs.low],col='lightskyblue',lwd=0.5,lty=2)

# Etreme: Step 2
plot(names(sim.cut),sim.cut,type='l',ylim=extendrange(sim.cut,f=0.15),
     ylab='',xlab='',xaxt='n',cex.axis=1.8)

interval = 5
at = seq(1,length(sim.cut),by=interval)
axis(1, at=names(sim.cut)[at], labels = names(sim.cut)[at],cex.axis=1.8)
mtext(side=2,'Simulated tree from the \n biggest diameter class (mm)',line=3,cex=1.3)
mtext(side=1,'Calendar year',line=3.5,cex=1.3)

points(names(sim.cut[ind.obs.high]),sim.cut[ind.obs.high],pch=19,cex=0.7)
points(names(sim.cut[ind.obs.low]),sim.cut[ind.obs.low],pch=19,cex=0.7)

segments(as.numeric(names(sim.cut[ind.obs.high])),-1,as.numeric(names(sim.cut[ind.obs.high])),
         sim.cut[ind.obs.high],col='indianred1',lwd=0.5,lty=2)
segments(as.numeric(names(sim.cut[ind.obs.low])),-1,as.numeric(names(sim.cut[ind.obs.low])),
         sim.cut[ind.obs.low],col='lightskyblue',lwd=0.5,lty=2)

axis(side=1,at=names(obs.cut[ind.obs.high]),tcl=0.4,lwd.ticks=3,mgp=c(0,0.5,0),label=NA,
     col.ticks = 'indianred1')
axis(side=1,at=names(obs.cut[ind.obs.low]),tcl=0.4,lwd.ticks=3,mgp=c(0,0.5,0),label=NA,
     col.ticks = 'lightskyblue')

# Extreme: Step3
plot(obs.extr,sim.extr,pch=19,xlab='',ylab='',cex.axis=1.8,
     xlim=extendrange(range(obs.extr,sim.extr),f=0.1),
     ylim=extendrange(range(obs.extr,sim.extr),f=0.1))
mtext(side=2,'Simulated output at \n extreme years (Normalized)',line=3,cex=1.3)
mtext(side=1,'Observation, 25% and 75% extreme (Normalized)',line=3.5,cex=1.3)

rect(min(scale(obs.cut)[ind.obs.high,1]),-10,max(scale(obs.cut)[ind.obs.high,1]),10,
     density=10,border=NA,col='indianred1')
rect(min(scale(obs.cut)[ind.obs.low,1]),-10,max(scale(obs.cut)[ind.obs.low,1]),10,
     density=10,border=NA,col='lightskyblue')

points(obs.extr,sim.extr,pch=19)

# Extreme: Step4
plot(obs.extr,sim.extr,pch=19,xlab='',ylab='',cex.axis=1.8,
     xlim=extendrange(range(obs.extr,sim.extr),f=0.1),
     ylim=extendrange(range(obs.extr,sim.extr),f=0.1))
mtext(side=2,'Simulated tree from the \n biggest diameter class (mm)',line=3,cex=1.3)
mtext(side=1,'Observation, 25% and 75% extreme (Normalized)',line=3.5,cex=1.3)

abline(a=0,b=1,col='green4')

arrows(obs.extr,sim.extr,obs.extr,obs.extr,length=0.08,angle=30,
       code=3)

# Extreme: Step5
plot(names(sim.cut),sim.cut,type='l',ylim=extendrange(sim.cut,f=0.15),
     ylab='',xlab='',xaxt='n',cex.axis=1.8)

interval = 5
at = seq(1,length(sim.cut),by=interval)
axis(1, at=names(sim.cut)[at], labels = names(sim.cut)[at],cex.axis=1.8)
mtext(side=2,'Simulated tree from the \n biggest diameter class (mm)',line=3,cex=1.3)
mtext(side=1,'Calendar year',line=3.5,cex=1.3)

rect(1900,min(sim.cut[ind.sim.high]),2100,max(sim.cut[ind.sim.high]),
     density=10,border=NA,col='indianred1')
rect(1900,min(sim.cut[ind.sim.low]),2100,max(sim.cut[ind.sim.low]),
     density=10,border=NA,col='lightskyblue')
points(names(sim.cut[ind.sim.high]),sim.cut[ind.sim.high],pch=19,cex=0.7)
points(names(sim.cut[ind.sim.low]),sim.cut[ind.sim.low],pch=19,cex=0.7)

# Extreme: Step6
plot((obs.rank-mean(obs.rank)),(sim.rank-mean(sim.rank)),pch=19,xlab='',ylab='',cex.axis=1.8,
     xlim=extendrange(range(obs.rank-mean(obs.rank),sim.rank-mean(sim.rank)),f=0.1),
     ylim=extendrange(range(obs.rank-mean(obs.rank),sim.rank-mean(sim.rank)),f=0.1))
mtext(side=2,'Simulation, 25% and 75% \n extreme (residual of the mean)',line=3,cex=1.3)
mtext(side=1,'Observation, 25% and 75% extreme',line=3.5,cex=1.3)

abline(a=0,b=1,col='green4')

arrows(obs.rank-mean(obs.rank),sim.rank-mean(sim.rank),
       obs.rank-mean(obs.rank),obs.rank-mean(obs.rank),length=0.08,angle=30,
       code=3)

# legend 
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend('top',c('Compiled virtual tree','Simulated biggest output','1:1 line',NA,
               'Grwoth range \n exceeds the 75 percentile','Growth range \n below the 25 percentile'),
       col=c('black','dodgerblue2','green4','white','white','white'),
       lty=c(1,1,1,NA,NA,NA),fill=c('black','dodgerblue2','green4','white','indianred1','lightskyblue'),
       density = c(0,0,0,0,10,10),border=c(NA,NA,NA,NA,'indianred1','lightskyblue'),
       ncol=3,text.width = c(0.17,0.17,0.17),cex=1.3,
       y.intersp=1.2,adj=0)

dev.off()
