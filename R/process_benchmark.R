# The code was written to build 4 benchmarks with 2 metrics using ITRDB datasets,
# and visualizing the result.
# This code can reproduce the result (Fig. 9) in the research paper :
# 'Using the International Tree-Ring Data Bank (ITRDB) records as century-long benchmarks 
# for land-surface models'
# 
# Author: Jina Jeong (j.jeong@vu.nl)

rm(list=ls())
library(ncdf4)
library(dplR)
library(plotrix)
library(dplyr)

###############################################################################################
#######################################   Initialize   ########################################

site_list = c('brit019','brit021','finl039','finl052','fran6','germ153','germ214',
              'neth034','spai006','swit188')
col.list=c('black','dodgerblue','indianred1','palegreen2','orange')
pft_list = list(brit019=4,brit021=4,finl039=4,finl052=4,fran6=4,
                germ153=4,germ214=4,neth034=4,spai006=4,swit188=4)
PFT_name = c('SoilBare','BroadLeavedEvergreenTropical','BroadLeavedRaingreenTropical',
             'NeedleleafEvergreenTemperate','BroadLeavedEvergreenTemperate','BroadLeavedSummergreenTemperate',
             'NeedleleafEvergreenBoreal','BroadLeavedSummergreenBoreal','LarixSpBoreal','C3GrassTemperate',
             'C4GrassTemperate','C3AgricultureTemperate','C4Agriculture')

site.info=read.csv2('data/10site_info.csv')
rownames(site.info)=site.info$site

config_list = c('Ndyn','recru','power','basic') 
target_circ=5

setwd('~/surfdrive/work/J.Jeong_GMD_2020/')

## Read functions 
source('R/functions_J.R')

## Load data and simulation output
## each file is the list contains 10 site ITRDB rwl or simulated output
load('data/obs.RData')
for (config in 1:length(config_list)){
  load(paste0('output/simulation.',config_list[config],'.RData'))
  
}

######################################## Compare trend ##################################
var_list=c('slope','err_slope','rmse','err_rmse')

trend_array<-array(0,dim=c(length(site_list),length(var_list),length(config_list)),
                   dimnames = list(site_list,var_list,config_list))

for (config in config_list){
  
  simulation.list = get(paste0('simulation.',config))
  data.trend = data.frame(matrix(NA,ncol=length(var_list),nrow=length(site_list))) %>%
    `colnames<-`(var_list) %>% `rownames<-`(site_list)
  
  for (list in 1:length(site_list)) {
    
    site_name = site_list[list]
    simulation=simulation.list[[site_name]]
    obs_rwl=obs[[site_name]]
    simulation=as.data.frame(simulation)
    
    # Leave-one-out arppoach
    data.trend.sub = data.frame(slope = rep(NA,ncol(obs_rwl)),rmse = rep(NA,ncol(obs_rwl)))
    for (i in c(1:ncol(obs_rwl))){
      
      obs_sub = obs_rwl[,-i]
      data.obs = data.frame(t=c(1:dim(obs_sub)[1]),data= apply(obs_sub,1,tbrm))
      data.obs = data.frame(t=c(1:dim(obs_sub)[1]),
                            data= apply(apply(obs_sub,2,align.trw),1,tbrm))
      data.sim = data.frame(t=c(1:nrow(simulation)),data = simulation[,target_circ])
      sub = (data.obs - data.sim)[,2]
      t = c(1:length(sub))
      test = lm(sub ~ t) # linear regression
      
      rcs.obs = rcs(obs_sub,
                    po = data.frame(ids = colnames(obs_sub), po = rep(1,ncol(obs_sub))),
                    rc.out=T,make.plot=F)
      
      rcs.sim = rcs(simulation,
                    po = data.frame(ids = colnames(simulation), 
                                    po = rep(1,ncol(simulation))),
                    rc.out=T,make.plot=F)
      
      # calculation for Fig.9 a.
      data.trend.sub$slope[i] = summary(test)$coefficients[2,1]
      data.trend.sub$rmse[i] = sqrt(sum((simulation[,target_circ]-data.obs)^2,na.rm=T)/length(t))
      
    }
    
    data.trend$slope[list] <- mean(data.trend.sub$slope)
    data.trend$err_slope[list] <- error_margin(data.trend.sub$slope,0.95)
    data.trend$rmse[list] <- mean(data.trend.sub$rmse)
    data.trend$err_rmse[list] <- error_margin(data.trend.sub$rmse,0.95)
    
  }
  
  trend_array[,,config]=as.matrix(data.trend)
  
}


##########################################################################################
####################### Compare diameter by CALENDAR (MATURE) ############################

var_list=c('rmse','err_rmse','slope','err_slope')
dia_array<-array(0,dim=c(length(site_list),length(var_list),length(config_list)),
                 dimnames = list(site_list,var_list,config_list))

for (config in 1:length(config_list)){
  
  config=config_list[config]
  simulation = get(paste0('simulation.',config))
  data.dia.calendar <- as.data.frame(matrix(NA,nrow=length(site_list),ncol=length(var_list))) %>%
    'colnames<-'(var_list) %>% 'rownames<-'(site_list)
  
  for (list in 1:length(site_list)) {
    
    site_name = site_list[list]
    obs_rwl=obs[[site_name]]
    trw.sim = simulation[[site_name]]
    
    #Leave-one-out
    data.dia.sub = data.frame(rmse = rep(NA,ncol(obs_rwl)),slope = rep(NA,ncol(obs_rwl)))
    for (i in c(1:ncol(obs_rwl))){
      obs_sub = obs_rwl[,-i]
      
      dia.end.obs = apply(obs_sub,2,cumsum.j)
      dia.end.mean.obs = apply(dia.end.obs,1,mean,na.rm=T)[51:nrow(obs_sub)] #mean line
      
      dia.max.sim = (2*cumsum(trw.sim[,target_circ]))[51:nrow(obs_sub)] #maximum sim
      test.dia = lm((dia.end.mean.obs-dia.max.sim)~c(1:length(c(51:nrow(obs_sub)))))
      
      #Fill data for Fig. 9b
      data.dia.sub$rmse[i] = sqrt(sum((dia.end.mean.obs-dia.max.sim)^2)/length(dia.max.sim))
      data.dia.sub$slope[i] = summary(test.dia)$coefficients[2,1]
      
    }
    
    data.dia.calendar$rmse[list] = mean(data.dia.sub$rmse)
    data.dia.calendar$err_rmse[list] = error_margin(data.dia.sub$rmse,0.95)
    data.dia.calendar$slope[list] = mean(data.dia.sub$slope)
    data.dia.calendar$err_slope[list] = error_margin(data.dia.sub$slope,0.95)
    
  } 
  dia_array[,,config] = as.matrix(data.dia.calendar)
}


##########################################################################################
######################## Compare diameter by AGE (YOUNG) #################################

var_list=c('rmse','err_rmse','slope','err_slope')
age_array<-array(0,dim=c(length(site_list),length(var_list),length(config_list)),
                 dimnames = list(site_list,var_list,config_list))

for (config in 1:length(config_list)){
  
  config=config_list[config]
  simulation = get(paste0('simulation.',config))
  data.dia.age = data.frame(matrix(NA,nrow=length(site_list),ncol=length(var_list))) %>%
    `row.names<-`(site_list) %>% `colnames<-`(var_list)
  rownames(data.dia.age) = site_list
  
  for (list in 1:length(site_list)) {
    
    site_name = site_list[list]
    obs_rwl = obs[[site_name]]
    trw.sim=simulation[[site_name]]
    
    # Leave-one-out
    data.dia.sub = data.frame(rmse = rep(NA,ncol(obs_rwl)),slope = rep(NA,ncol(obs_rwl)))
    for ( i in c(1:ncol(obs_rwl))){
      obs_sub = obs_rwl[,-i]
      
      trw.age = apply(obs_sub,2,align.trw)
      dia.age.obs = apply(trw.age,2,cumsum.j)
      dia.age.mean.obs = apply(dia.age.obs,1,max,na.rm=T)[1:30]
      dia.age.max.sim= (2*cumsum(trw.sim[,target_circ]))[1:30]
      
      test.age = lm((dia.age.mean.obs-dia.age.max.sim)~c(1:30))
      # Fill data for Fig. 9c
      data.dia.sub$slope[i] = summary(test.age)$coefficients[2,1]
      data.dia.sub$rmse[i] = sqrt(sum((dia.age.mean.obs-dia.age.max.sim)^2)/30)
    } 
    data.dia.age$rmse[list] = mean(data.dia.sub$rmse)
    data.dia.age$err_rmse[list] = error_margin(data.dia.sub$rmse,0.95)
    data.dia.age$slope[list] = mean(data.dia.sub$slope)
    data.dia.age$err_slope[list] = error_margin(data.dia.sub$rmse,0.95)
    
  }
  age_array[,,config] = as.matrix(data.dia.age)
  
}

##########################################################################################
################################ Extreme growths (EXTREME) ##############################

# p for defining extreme growth
p=0.25

var_list=c('amp','err_amp','extr','err_extr')
extr_array =array(dim=c(length(site_list),length(var_list),length(config_list)),
                  dimnames=list(site_list,var_list,config_list))

for (config in config_list){
  simulation=get(paste0('simulation.',config))
  
  for (site_name in site_list){
    obs_rwl=obs[[site_name]]
    
    # Leave-one-out
    extr.sub = data.frame(amp = rep(NA,ncol(obs_rwl)),extr=rep(NA,ncol(obs_rwl)))
    for (i in c(1:ncol(obs_rwl))){
      obs_sub = obs_rwl[,-i]
      
      # Choose long seriese to avoid age trends
      obs.long=obs_sub[,which(apply(obs_sub,2,function(x) sum(is.na(x))<20)==T)]
      obs.cut=obs.long[which(as.numeric(rownames(obs.long))>1950),]
      obs.cut=rowMeans(obs.cut,na.rm=T)
      
      sim.cut=simulation[[site_name]][which(as.numeric(rownames(obs.long))>1950),target_circ]
      names(sim.cut)<-names(obs.cut)
      
      ind.obs.low = order(obs.cut)[1:floor(p*length(obs.cut))]
      ind.obs.high= rev(order(obs.cut,decreasing=T)[1:floor(p*length(obs.cut))])
      ind.sim.low = order(sim.cut)[1:floor(p*length(sim.cut))]
      ind.sim.high= rev(order(sim.cut,decreasing=T)[1:floor(p*length(sim.cut))])
      
      # Choose by extreme years
      obs.extr = scale(obs.cut)[c(ind.obs.low,ind.obs.high),1]
      sim.extr = scale(sim.cut)[c(ind.obs.low,ind.obs.high),1]
      # Choose by extreme growth
      obs.rank = obs.cut[c(ind.obs.low,ind.obs.high)]
      sim.rank = sim.cut[c(ind.sim.low,ind.sim.high)]
      
      rmse.extr = sqrt(mean((sim.extr-obs.extr)^2)/length(sim.rank))
      rmse.rank = sqrt(mean((sim.rank-obs.rank)^2)/length(sim.rank))
      # Fill data for Fig. 9d
      extr.sub$amp[i] = sqrt(mean((sim.rank-obs.rank)^2)/length(sim.rank))
      extr.sub$extr[i] = sqrt(mean((sim.extr-obs.extr)^2)/length(sim.rank))
    }
    
    extr_array[site_name,'amp',config]=mean(extr.sub$amp)
    extr_array[site_name,'err_amp',config]=error_margin(extr.sub$amp,0.95)
    extr_array[site_name,'extr',config] = mean(extr.sub$extr)
    extr_array[site_name,'err_extr',config] = error_margin(extr.sub$extr,0.95)
    
  }
}    


##################################################################################
################################ FIGURE_COMBIME ##################################

rgb.list=col2rgb(col.list)/255

## To save figure
#png('figure/comb_benchmarks.png',width=11.5, 
#    height=9.8, units="in", res=300)
svg('figure/comb_benchmarks.svg',width=12, height=10)

##! bunch of errors can be made due to too short error bars

par(mfrow=c(2,2))
par(mar=c(5.5,4.9,2.7,2.1))


## Trend
for ( i in c(1:length(config_list))){
  config = config_list[i]
  data=trend_array[,,config]
  
  if (config=='Ndyn'){
    plot(-100,-100,
         xlab='',
         ylab='',
         ylim=c(extendrange(trend_array[,'slope',],f=0.1)[1],1.5*max(trend_array[,'slope',])),
         xlim =c(extendrange(trend_array[,'rmse',],f=0.16)[1],
                 extendrange(trend_array[,'rmse',],f=0.08)[2]),
         xaxt='n')
    
    arrows(data[,'rmse']-data[,'err_rmse'], data[,'slope'], 
           data[,'rmse']+data[,'err_rmse'], data[,'slope'], length=0.02, angle=90, code=3)
    arrows(data[,'rmse'],data[,'slope']-data[,'err_slope'],
           data[,'rmse'],data[,'slope']+data[,'err_slope'], length=0.02, angle=90, code=3)
    
    points(data[,'rmse'],data[,'slope'],pch=19)
    axis.break(1,breakpos=extendrange(trend_array[,'rmse',],f=0.09)[1])
    start = min(trend_array[,'rmse',])-min(trend_array[,'rmse',])%%5
    end = max(trend_array[,'rmse',])-(max(trend_array[,'rmse',]))%%5
    axis(1,at=seq(start,end,length.out=5),
         label=seq(start,end,length.out=5))
    axis(1,at=extendrange(trend_array[,'rmse',],f=0.16)[1],label=0)
    
    text(data[,'rmse'],data[,'slope'],pos=4,cex=1.2,label=c(1:10))
  } else { # if config=='Ndyn'
    arrows(data[,'rmse']-data[,'err_rmse'], data[,'slope'], 
           data[,'rmse']+data[,'err_rmse'], data[,'slope'], length=0.02, angle=90, code=3, col='azure4')
    arrows(data[,'rmse'],data[,'slope']-data[,'err_slope'],
           data[,'rmse'],data[,'slope']+data[,'err_slope'], length=0.02, angle=90, code=3, col='azure4')
    
    points(data[,'rmse'],data[,'slope'],col=col.list[i],pch=19)
    segments(data[,'rmse'],data[,'slope'],trend_array[,'rmse',i-1],trend_array[,'slope',i-1],
             lty=3,lwd=0.8)
    
  }
  
}

mtext('RMSE between \n simulated and observed trends (mm/yr)',1,line=3.5,cex=1.1)
mtext('Temporal slope of the residuals \n from tree-ring width by tree age (mm/yr)',2,line=2.2,cex=1.1)
mtext(expression(bold('Size-related trend')),3,cex=1.2,line=0.4)
legend('topright',config_list,col=col.list[1:length(config_list)],pch=19,bty='n')

## Mature
for ( i in 1:length(config_list)){
  config=config_list[i]
  data=dia_array[,,config]
  cex=rep(1,length(site_list))
  #cex[which(data[,'p_slope']>0.05)]=1.15
  if (config=='Ndyn'){
    plot(data[,'rmse'],data[,'slope'],pch=19,
         xlab='',
         ylab='',
         ylim=extendrange(dia_array[,'slope',],f=0.08),
         xlim = c(extendrange(dia_array[,'rmse',],f=0.12)[1],
                  extendrange(dia_array[,'rmse',],f=0.08)[2]),
         xaxt='n')
    
    arrows(data[,'rmse']-data[,'err_rmse'], data[,'slope'], 
           data[,'rmse']+data[,'err_rmse'], data[,'slope'], length=0.02, angle=90, code=3, col='azure4')
    arrows(data[,'rmse'],data[,'slope']-data[,'err_slope'],
           data[,'rmse'],data[,'slope']+data[,'err_slope'], length=0.02, angle=90, code=3, col='azure4')
    points(data[,'rmse'],data[,'slope'],pch=19)
    axis.break(1,breakpos=extendrange(dia_array[,'rmse',],f=0.08)[1])
    start = min(dia_array[,'rmse',])-min(dia_array[,'rmse',])%%10
    end = max(dia_array[,'rmse',])-(max(dia_array[,'rmse',]))%%10
    axis(1,at=seq(start,end,length.out=5),
         label=seq(start,end,length.out=5))
    axis(1,at=extendrange(range(min(dia_array[,'rmse',]),max(dia_array[,'rmse',])),f=0.12)[1]
         ,label=0)
    
    
  } else { # if config=='Ndyn'
    arrows(data[,'rmse']-data[,'err_rmse'], data[,'slope'], 
           data[,'rmse']+data[,'err_rmse'], data[,'slope'], length=0.02, angle=90, code=3, col='azure4')
    arrows(data[,'rmse'],data[,'slope']-data[,'err_slope'],
           data[,'rmse'],data[,'slope']+data[,'err_slope'], length=0.02, angle=90, code=3, col='azure4')
    points(data[,'rmse'],data[,'slope'],col=col.list[i],pch=19)
    segments(data[,'rmse'],data[,'slope'],dia_array[,'rmse',i-1],dia_array[,'slope',i-1],
             lty=3,lwd=0.8)
    
  }
  
}

text(dia_array[,'rmse','Ndyn'],dia_array[,'slope','Ndyn'],pos=4,cex=1.2,label=c(1:10))

mtext('RMSE between simulated and observed \n diameter for mature trees (mm/yr)',1,line=3.5,cex=1.1)
mtext('Temporal slope of the residuals from
      diameter growth for matured trees (mm/yr)',2,line=2.2,cex=1.1)
mtext(expression(bold('Mature trees')),3,cex=1.2,line=0.4)

## Young
for ( i in 1:length(config_list)){
  config=config_list[i]
  data=age_array[,,config]
  cex=rep(1,length(site_list))
  #cex[which(data[,'p_slope']>0.05)]=1.15
  if (config=='Ndyn'){
    plot(data[,'rmse'],data[,'slope'],pch=19,
         xlab='',
         ylab='',
         ylim=extendrange(age_array[,'slope',],f=0.08),
         xlim = c(extendrange(age_array[,'rmse',],f=0.1)[1],
                  extendrange(age_array[,'rmse',],f=0.1)[2]),
         xaxt='n')
    arrows(data[,'rmse']-data[,'err_rmse'], data[,'slope'], 
           data[,'rmse']+data[,'err_rmse'], data[,'slope'], length=0.02, angle=90, code=3, col='azure4')
    arrows(data[,'rmse'],data[,'slope']-data[,'err_slope'],
           data[,'rmse'],data[,'slope']+data[,'err_slope'], length=0.02, angle=90, code=3, col='azure4')
    points(data[,'rmse'],data[,'slope'],pch=19)
    
    start = min(age_array[,'rmse',])-min(age_array[,'rmse',])%%10
    end = max(age_array[,'rmse',])-(max(age_array[,'rmse',]))%%10
    axis(1,at=seq(start,end,length.out=5),
         label=round(seq(start,end,length.out=5)))
    
  } else {# if config=='Ndyn'
    arrows(data[,'rmse']-data[,'err_rmse'], data[,'slope'], 
           data[,'rmse']+data[,'err_rmse'], data[,'slope'], length=0.02, angle=90, code=3, col='azure4')
    arrows(data[,'rmse'],data[,'slope']-data[,'err_slope'],
           data[,'rmse'],data[,'slope']+data[,'err_slope'], length=0.02, angle=90, code=3, col='azure4')
    
    points(data[,'rmse'],data[,'slope'],col=col.list[i],pch=19)
    segments(data[,'rmse'],data[,'slope'],age_array[,'rmse',i-1],age_array[,'slope',i-1],
             lty=3,lwd=0.8)
    
  }
  
}
text(age_array[,'rmse','Ndyn'],age_array[,'slope','Ndyn'],pos=4,cex=1.2,label=c(1:10))

mtext('RMSE between simulated and \n observed diameter for young trees (mm/yr)',1,line=3.5,cex=1.1)
mtext('Temporal slope of the residuals from
      diameter growth for young trees (mm/yr)',2,line=2.2,cex=1.1)
mtext(expression(bold('Young trees')),3,cex=1.2,line=0.4)


## Extreme
plot(-1,-1,col='white',xlab='',ylab='',
     ylim=extendrange(extr_array[,'amp',],f=0.08),
     xlim = c(0,extendrange(extr_array[,'extr',],f=0.1)[2]),
     xaxt='n')

start = 0
end = max(extr_array[,'extr',])-max(extr_array[,'extr',])%%0.05
axis(1,at=seq(start,end,length.out=5),
     label=seq(start,end,length.out=5))

for(i in c(1:length(config_list))){
  config=config_list[i]
  
  arrows(extr_array[,'extr',config]-extr_array[,'err_extr',config], extr_array[,'amp',config], 
         extr_array[,'extr',config]+extr_array[,'err_extr',config], extr_array[,'amp',config], 
         length=0.02, angle=90, code=3, col='azure4')
  arrows(extr_array[,'extr',config],extr_array[,'amp',config]-extr_array[,'err_amp',config],
         extr_array[,'extr',config],extr_array[,'amp',config]+extr_array[,'err_amp',config], 
         length=0.02, angle=90, code=3, col='azure4')
  
  points(extr_array[,'extr',config],extr_array[,'amp',config],col=col.list[i],pch=19)
  if(i>1){
    segments(extr_array[,'extr',i],extr_array[,'amp',i],
             extr_array[,'extr',i-1],extr_array[,'amp',i-1],
             lty=3,lwd=0.8)
  }
}

text(extr_array[,'extr','Ndyn'],extr_array[,'amp','Ndyn'],
     label=c(1:10),pos=1,cex=1.2)

mtext('RMSE between simulated and observed 25%
      extreme tree-ring width regardless of year (mm/yr)',2,line=2.2,cex=1.1)
mtext('RMSE between simulated and observed 
      tree-ring width at years with extreme growth (scaled)',1,line=3.5,cex=1.1)
mtext(expression(bold('Extreme growth events')),3,cex=1.2,line=0.4)


## To save figure
dev.off()



