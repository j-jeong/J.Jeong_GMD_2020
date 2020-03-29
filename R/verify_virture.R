# The code was written to verify one assumption laid under 4 benchmarks, 
# the overestimation through big-tree selection.
# This code can reproduce Fig. 4 in the research paper :
# 'Using the International Tree-Ring Data Bank (ITRDB) records as century-long benchmarks 
# for land-surface models'
# 
# To use this code, you need to download the BACI dataset.
# You can download it by registration with your email on the site :
# http://www.baci-h2020.eu/

# Author: Jina Jeong (j.jeong@vu.nl)

rm(list=ls())
library(dplyr)
library(dplR)
library(ggplot2)
library(reshape2)
library(plotly)


###############################################################################################
#######################################   Initialize   ########################################

Benchmarks= c('Trend','Mature','Young','Extreme')

## Data
file=list.files(
  '201961114214112_treeringdata_BACI2016_WP3/data/treeringbiomass_network_Europe/raw_data')
site.obs=unlist(strsplit(file,'[.]'))[seq(1,48*2,by=2)]
source('R/functions_J.R')

meta = read.table(
  '201961114214112_treeringdata_BACI2016_WP3/data/treeringbiomass_network_Europe/metadata/ABI_Europe_metadata_sitelevel_BACI2016.txt',
  header=T,stringsAsFactors = F)
meta.tree = read.table(
  '201961114214112_treeringdata_BACI2016_WP3/data/treeringbiomass_network_Europe/metadata/ABI_Europe_metadata_treelevel_BACI2016.txt',
  header=T,stringsAsFactors = F) %>% mutate(Sitecode = as.factor(Sitecode))
# number of trees per site
summary(meta.tree$Sitecode)
# dominant species 
unique(c(meta$DominantSpecies1,meta$DominantSpecies2))
# Select sites dominated by coniferous species
remove.sp = c('LADE','FREX','CASA','FASY')
select = meta$SiteCode[which(meta$DominantSpecies1 %in% remove.sp | meta$DominantSpecies2 %in% remove.sp)]

## Simulation

dist = c(5,7,9,7,5) # distribution for size classes used for simulations
target_circ = 5 # targetted size class
deleuze_p = 0.3 # delueze_p to calculate f_sigma in Eq. 27, Text S1.

site_list = c('brit019','brit021','finl039','finl052','fran6','germ153','germ214',
              'neth034','spai006','swit188')
config_list = c('Ndyn','recru','power','basic')
PFTidx=4


###############################################################################################
##################################   Build benchmarks_BACI   ##################################

# Calculate ratio for selecint trees depends on target circ
ratio_max = sum(dist[1:target_circ])/sum(dist)
ratio_min = sum(dist[1:(target_circ-1)])/sum(dist)

ratio_obs = matrix(nrow=length(select),rep(NA,4*length(select))) %>% 
  `row.names<-`(select) %>% `colnames<-`(c('ratio_trend','ratio_mature','ratio_young','ratio_extr'))
counted=0

for(i in select){
  trw = read.table(paste0(
    '201961114214112_treeringdata_BACI2016_WP3/data/treeringbiomass_network_Europe/raw_data/',
    i,'.txt'
  ),header=T,row.names=1) # read.data
  dia.yr=apply(trw,2,cumsum.j) # diameter aligned by age
  
  # Select sites which have closer stand structure to the model
  if(min(dia.yr[nrow(dia.yr),],na.rm=T)>90){
    counted = counted +1
    message(i,': no trees samller than 9 cm')
    next
  }
  
  trw.age = apply(trw,2,align.trw)
  dia.age = apply(trw.age,2,cumsum)*2
  
  # Make index for big trees in the observed site
  # This choice is based on the size distribution in the simulation
  choose=order(apply(dia.yr[1:nrow(dia.yr),],2,max,na.rm=T),
               decreasing=F)[round(ncol(dia.yr)*ratio_min):round(ncol(dia.yr)*ratio_max)]
  choose=choose[!is.na(choose)]
  
  ############## TREND
  idata=trw.age
  fun=mean #mean for dia.yr, trw.age, max for cumsum trw.age
  virtual = apply(idata,1,fun,na.rm=T)
  
  idata_big=idata[,choose]
  virtual_big = apply(idata_big,1,fun,na.rm=T)
  
  na_ind = !is.na(virtual_big)&!is.na(virtual)
  virtual=virtual[na_ind]
  virtual_big=virtual_big[na_ind]
  
  ratio_obs[i,'ratio_trend'] = mean(virtual_big/virtual,na.rm=T)
  
  ############## YOUNG
  idata=dia.age
  fun=max #mean for dia.yr, trw.age, max for cumsum dia.age
  virtual = apply(idata,1,max,na.rm=T)[1:50]
  
  idata_big=idata[,choose]
  virtual_big = apply(idata_big,1,max,na.rm=T)[1:50]
  
  na_ind = !is.na(virtual_big)&!is.na(virtual)
  virtual=virtual[na_ind]
  virtual_big=virtual_big[na_ind]
  
  ratio_obs[i,'ratio_young'] = mean(virtual_big/virtual)
  
  ############## MATURE
  idata=dia.yr
  fun=mean #mean for dia.yr, trw.age, max for cumsum trw.age
  virtual = apply(idata,1,fun,na.rm=T)[50:nrow(dia.yr)]
  
  idata_big=idata[,choose]
  virtual_big = apply(idata_big,1,fun,na.rm=T)[50:nrow(dia.yr)]
  
  na_ind = !is.na(virtual_big)&!is.na(virtual)
  virtual=virtual[na_ind]
  virtual_big=virtual_big[na_ind]
  
  ratio_obs[i,'ratio_mature'] = mean(virtual_big/virtual)
  
  ############## EXTR
  yr.ind = which(as.numeric(rownames(trw)) %in% c(1951:2000))
  growth.ind = which(apply(trw[yr.ind,],2,function(x) sum(!is.na(x))>(length(yr.ind)-10))==TRUE)
  ind.long=which(apply(trw[,growth.ind],2,function(x) sum(!is.na(x))>80)==TRUE)
  ind = intersect(growth.ind,ind.long)
  
  trw_cut= trw[which(as.numeric(rownames(trw))>1950),ind]
  
  mean.all = apply(trw_cut,1,tbrm)
  ind.extr.low = names(mean.all)[c(order(mean.all)[1:floor(0.25*length(mean.all))])]
  ind.extr.high <- names(mean.all)[rev(order(mean.all,decreasing=T)[1:floor(0.25*length(mean.all))])]
  ind.extr = c(ind.extr.low,ind.extr.high)
  mean.extr.all = mean.all[which(names(mean.all) %in% ind.extr)]
  
  trw_big=trw[,choose[!is.na(choose)]]
  
  growth.ind = which(apply(trw_big[yr.ind,],2,function(x) sum(!is.na(x))>(length(yr.ind)-10))==TRUE)
  ind.long=which(apply(trw_big[,growth.ind],2,function(x) sum(!is.na(x))>80)==TRUE)
  ind = intersect(growth.ind,ind.long)
  
  trw_cut_big= trw_big[which(as.numeric(rownames(trw_big))>1950),ind]
  
  # In case only one tree is picked
  if(length(ind)==1){
    mean.big=trw_cut_big
    names(mean.big) <-rownames(trw_big)[which(as.numeric(rownames(trw_big))>1950)]
  } else {
    mean.big = apply(trw_cut_big,1,tbrm)
  }
  ind.extr.low = names(mean.big)[c(order(mean.big)[1:floor(0.25*length(mean.big))])]
  ind.extr.high <- names(mean.big)[rev(order(mean.big,decreasing=T)[1:floor(0.25*length(mean.big))])]
  ind.extr = c(ind.extr.low,ind.extr.high)
  
  mean.extr.big = mean.big[which(names(mean.big) %in% ind.extr)]
  
  ratio_obs[i,'ratio_extr'] = mean(mean.extr.big/mean.extr.all)
  
} 

message(paste0('Omitted : ',counted,' sites.'))

colnames(ratio_obs)<-c('Trend','Mature','Young','Extreme')
ratio_obs <- ratio_obs[!is.na(rowSums(ratio_obs)),]
ratio_obs_long = melt(ratio_obs,varnames=c('site','Benchmark')) %>% mutate(Benchmark=as.factor(Benchmark))

# stat
mean_obs <- setNames(aggregate(value ~ Benchmark, ratio_obs_long,mean),c('Benchmark','mean'))
sd_obs <- setNames(aggregate(value ~ Benchmark, ratio_obs_long,sd),c('Benchmark','sd'))
stat_obs <- setNames(cbind(mean_obs,sd_obs[,'sd']),c('Benchmark','mean','sd'))


###############################################################################################
##################################   Build benchmarks_ORC  ####################################

ratio_sim_array = array(NA,dim=c(length(site_list),4,length(config_list)),
                   dimnames=list(site_list,Benchmarks,config_list))

for ( config in config_list){
  
  data.eval = data.frame(matrix(NA,ncol=length(Benchmarks),nrow=length(site_list))) %>%
    `colnames<-`(Benchmarks) %>% `rownames<-`(site_list)
  
  load(paste0('output/simulation.',config,'.RData'))
  simulation.all=get(paste0('simulation.',config))
  
  for ( site in site_list){
    simulation = simulation.all[[site]]

    All_trw = rowSums(t(t(simulation)*rep(dist,nrow(simulation))))/sum(dist)
    All_dia = rep(NA,length(All_trw))
    All_dia[!is.na(All_trw)] = 2*cumsum(All_trw[!is.na(All_trw)])
    
    Big_trw = simulation[,target_circ] %>% `names<-`(names(All_trw))
    Big_dia = cumsum(simulation[,target_circ])*2
    
    data.eval[site,'Trend'] = mean(Big_trw/All_trw,na.rm=T)
    data.eval[site,'Mature'] = mean(Big_dia[51:length(Big_dia)]/All_dia[51:length(All_dia)],na.rm=T)
    data.eval[site,'Young'] = mean(Big_dia[1:50]/Big_dia[1:50],na.rm=T)
    
    yr.ind = which(as.numeric(names(All_trw)) %in% c(1951:2000))
    All_cut = All_trw[yr.ind]
    Big_cut = Big_trw[yr.ind]
    
    ind.extr.low = names(All_trw)[c(order(Big_trw)[1:floor(0.25*length(Big_trw))])]
    ind.extr.low2 = names(Big_trw)[c(order(Big_trw)[1:floor(0.25*length(Big_trw))])]
    
    ind.extr.high <- names(All_trw)[rev(order(All_trw,decreasing=T)[1:floor(0.25*length(All_trw))])]
    ind.extr.high2 <- names(Big_trw)[rev(order(Big_trw,decreasing=T)[1:floor(0.25*length(Big_trw))])]
    
    ind.extr = c(ind.extr.low,ind.extr.high)
    ind.extr2 = c(ind.extr.low2,ind.extr.high2)
    
    data.eval[site,'Extreme'] = mean(Big_trw[ind.extr2]/All_trw[ind.extr])
    
  }
  
  assign(paste0(config,'_sim'),data.eval) 
  ratio_sim_array[,,config]<-as.matrix(data.eval)

}


vals_obs_long = melt(vals_obs,variable.name = 'Benchmark')
vals_obs_long$sd = as.vector(t(vals_obs_sd))
ratio_sim_long = melt(ratio_sim_array,varnames=c('Site','Benchmark','config'))


# Stat
mean_sim <- setNames(aggregate(value ~ Benchmark, ratio_sim_long,mean),c('Benchmark','mean'))
sd_sim <- setNames(aggregate(value ~ Benchmark, ratio_sim_long,sd),c('Benchmark','sd'))
stat_sim <- setNames(cbind(mean_sim,sd_sim[,'sd']),c('Benchmark','mean','sd'))

# Merge
stat_sim$type = 'Simulation'
stat_obs$type = 'Data'
stat_all=rbind(stat_sim,stat_obs)


###############################################################################################
#######################################   Make a plot  ########################################

ggplot(stat_all,aes(Benchmark,mean,col=type))+geom_point(size=4,position=position_dodge(w=0.7))+
  theme_bw()+geom_point(data=stat_all,
                        aes(Benchmark,mean,col=type),
                        size=0.7,alpha=0.5,position=position_dodge(w=0.7))+
  ylab('Big_trees/All_tree')+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)),
        axis.title.y=element_text(margin=margin(t=0,r=12,b=0,l=0)))+
  ylim(0.5,3)+
  geom_errorbar(data=stat_all,aes(x=Benchmark,ymin=mean-sd, ymax=mean+sd,col=type), 
                width=.2,position=position_dodge(.7))

#ggsave('figure/verify_benchmarks.png')
