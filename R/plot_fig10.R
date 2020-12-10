# Script to draw benchmarks for selected European biomass network sites 

# This script can reproduce Fig. 10 in the research paper:
# 'Using the International Tree-Ring Data Bank (ITRDB) records as century-long benchmarks 
# for land-surface models'

# Note that 
#   1) build_obs_Rdata.R needs to be run beforehand to create BACI_obs.RData
#   2) You can change sites to be plotted in L28

# Author: Jina Jeong (j.jeong@vu.nl)

#- 0. Initialization ####
source('R/functions_J.R')

# load simulation and obs
load('output/BACI_sim.r5698.RData')
load('data/BACI_obs.RData')

#- 1. Make plot ####
pdf("figure/sites_benchmarked.pdf", width=16,height=19)

layout.matrix <- matrix(c(1:20),nrow=5,byrow = F)
layout(mat=layout.matrix)

# The list can be modified by selecting four sites within : 
#   'DEO','DVN','GIU','HD2','SCH','SOB','TIC','CAN','SOR','TER','ZOF')
site_select = c('DEO','DVN','CAN','SOR')

for(site_name in site_select){
  
  trw.sim = BACI_sim.r5698[[site_name]]
  trw.obs = BACI_obs[[site_name]]
  
  #-- 1.1 Processing OBSERVATED data ####
  dia.end.obs = apply(trw.obs,2,cumsum.j)
  dia.sim = apply(trw.sim,2,cumsum.j)
  trw.age = apply(trw.obs,2,align.trw)
  dia.age.obs = apply(trw.age,2,cumsum.j)
  nyears=nrow(trw.obs)
  
  #-- 1.2 Trend ####
  sub = apply(trw.age,1,mean,na.rm=T)-trw.sim[,5]
  sub = sub[!is.na(sub)]
  reg = lm(sub~c(1:length(sub)))
  
  par(mar=c(5.1,4.1,4.1,2.1))
  plot(-100,-100,ylim = c(min(sub),max(sub)),xlim=c(0,nyears),
       xlab = 'Tree age (years)',cex.lab=1.6,cex.axis = 1.4
       ,ylab='',type='l',col='azure3',xaxt = "n",main=site_name,cex.main=2,
       cex.lab=1.6)
  
  interval = 40
  at = seq(1,nyears,by=interval)
  axis(1, at=seq(1,nyears,by=interval), labels = at,cex.axis=1.4)
  mtext(side=2,'Observation - simulation (mm)',line=2.5,cex=1.1)
  
  lines(sub,type='l',col='green3',lwd=1.5)
  abline(reg,col='green4',lwd=1.7,lty=2)
  
  #-- 1.3 Mature ####
  sub = apply(dia.end.obs,1,mean,na.rm=T)[51:nyears]-dia.sim[51:nyears,5]
  reg = lm(sub~c(51:nyears))
  
  plot(-100,-100,ylim = c(min(sub),max(sub)),xlim=c(51,nyears),
       xlab = 'Calendar year',cex.lab=1.1,cex.axis = 1.4
       ,ylab='',type='l',col='azure3',xaxt = "n",cex.lab=1.6)
  
  interval = 20
  at = seq(53,nyears,by=interval)
  axis(1, at=seq(53,nyears,by=interval), labels = row.names(trw.obs)[at],cex.axis=1.4)
  mtext(side=2,'Observation - simulation (mm)',line=2.5,cex=1.1)
  
  lines(c(51:nyears),sub,type='l',col='green3',lwd=1.5)
  
  abline(reg,col='green4',lwd=1.7,lty=2)
  
  #-- 1.4 Young ####
  sub = apply(dia.age.obs,1,max,na.rm=T)[1:30]-dia.sim[1:30,5]
  reg = lm(sub~c(1:30))
  
  plot(-100,-100,ylim = c(min(sub),max(sub)),xlim=c(0,30),
       xlab = 'Tree age (years)',cex.lab=1.6,cex.axis = 1.4
       ,ylab='',type='l',col='azure3',xaxt = "n",cex.lab=1.6)
  
  interval = 5
  at = seq(5,nyears,by=interval)
  axis(1, at=seq(5,nyears,by=interval), labels = at,cex.axis=1.4)
  mtext(side=2,'Observation - simulation (mm)',line=2.5,cex=1.1)
  
  lines(sub,type='l',col='green3',lwd=1.5)
  abline(reg,col='green4',lwd=1.7,lty=2)
  
  
  #-- 1.5 Extreme ####
  p=0.25
  target_circ=5
  
  # Pick trees aged over 50 and survived after 1950
  idx.age = !is.na(trw.obs[rownames(trw.obs)==1950,]) & 
    apply(trw.obs,2,function(x) sum(!is.na(x))) >50
  obs.long=trw.obs[,idx.age]
  obs.cut.all=obs.long[which(as.numeric(rownames(obs.long))>1950),]
  obs.cut=rowMeans(obs.cut.all,na.rm=T)
  
  sim.cut = trw.sim[which(as.numeric(rownames(obs.long))>1950),target_circ]
  
  # Get indice for extreme values
  ind.obs.low = order(obs.cut)[1:floor(p*length(obs.cut))]
  ind.obs.high= rev(order(obs.cut,decreasing=T)[1:floor(p*length(obs.cut))])
  ind.sim.low = order(sim.cut)[1:floor(p*length(sim.cut))]
  ind.sim.high= rev(order(sim.cut,decreasing=T)[1:floor(p*length(sim.cut))])
  
  # Extreme years based on observation
  obs.extr = scale(obs.cut)[c(ind.obs.low,ind.obs.high),1]
  sim.extr = scale(sim.cut)[c(ind.obs.low,ind.obs.high),1]
  
  # Extreme growth
  obs.rank = obs.cut[c(ind.obs.low,ind.obs.high)]
  sim.rank = sim.cut[c(ind.sim.low,ind.sim.high)]
  obs.rank = obs.rank-mean(obs.rank)
  sim.rank = sim.rank-mean(sim.rank)
  
  rmse.extr = sqrt(mean((sim.extr-obs.extr)^2)/length(sim.rank))
  rmse.rank = sqrt(mean((sim.rank-obs.rank)^2)/length(sim.rank))
  
  #--- 1.5.1 Extreme_ extreme growth ####
  par(mar=c(5.1,5.5,4.1,2.1))
  
  plot(obs.extr,sim.extr,pch=19,
       xlab='',ylab='',cex.lab=1.5,cex.axis=1.3,
       xlim=extendrange(range(obs.extr,sim.extr),f=0.1),
       ylim=extendrange(range(obs.extr,sim.extr),f=0.1))
  mtext(side=2,'Simulated tree from \n the biggest diameter class (mm)',line=2.2,cex=1.1)
  mtext(side=1,'Observation, 25% and 75% extreme \n (Normalized)',line=5,cex=1.1)
  
  abline(a=0,b=1,col='green4')
  
  arrows(obs.extr,sim.extr,obs.extr,obs.extr,length=0.08,angle=30,
         code=3)
  
  #--- 1.5.2 Extreme_ amplitude ####
  plot((obs.rank-mean(obs.rank)),(sim.rank-mean(sim.rank)),pch=19,
       xlab='Observation, 25% and 75% extreme',ylab='',cex.lab=1.5,cex.axis=1.3,
       xlim=extendrange(range(obs.rank-mean(obs.rank),sim.rank-mean(sim.rank)),f=0.1),
       ylim=extendrange(range(obs.rank-mean(obs.rank),sim.rank-mean(sim.rank)),f=0.1),
       cex.lab=1.6)
  mtext(side=2,'Simulation, 25% and 75% extreme \n (residual of the mean)',line=2.2,cex=1.1)
  
  abline(a=0,b=1,col='green4')
  
  arrows(obs.rank-mean(obs.rank),sim.rank-mean(sim.rank),
         obs.rank-mean(obs.rank),obs.rank-mean(obs.rank),length=0.08,angle=30,
         code=3)
  
}

dev.off()
