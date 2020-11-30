# Apply modifires to the simulated outputs for the verification of benchmarks.
# This code contains 
# 1) cacluation the modifier for the optimum simulation against big-tree data 
#    using 3 different methods
# 2) apply the same modifier to all-tree data
# 3) calculating the improvements
# , which is performed 8 times (two metrics for each of four benchmarks)
# 
# This code can reproduce the result (Table 2) in the research paper :
# 'Using the International Tree-Ring Data Bank (ITRDB) records as century-long benchmarks 
# for land-surface models'
# 
# Author: Jina Jeong (j.jeong@vu.nl)


#- 0. Initialize ####
site_info = data.frame(site=c('DEO','DVN','GIU','HD2','SCH','SOB','TIC','CAN','SOR','TER','ZOF'),
                       specie = c('PCAB','PCAB','PCAB','PCAB','PCAB','PCAB','PCAB',
                                  'FASY','FASY','FASY','FASY'),
                       PFT=c(4,4,4,4,4,4,4,6,6,6,6))
load('data/BACI_obs.RData')
load('output/BACI_sim.r5698.RData')
load('output/BACI_sim.ccn.r5698.RData')
source('R/functions_J.R')

BACI_sim = BACI_sim.r5698
BACI_sim_ccn = BACI_sim.ccn.r5698

site_list = site_info$site
circ_dist = c(5,7,9,7,5) # distribution of circumference classes used for the simulation
target_circ = 5 # target class for the simulation_big
age_young = 30 # age cut for calculationg young trees benchmarks
p_extr = 0.25 # ratio to calculate extreme growth

# Calculate the ratio to select big trees based on the simulated distribution
ratio_max = sum(circ_dist[1:target_circ])/sum(circ_dist)
ratio_min = sum(circ_dist[1:(target_circ-1)])/sum(circ_dist)

colnames = c(
  apply(expand.grid(apply(expand.grid(c('trend','mature','young'),c('rmse','slope')),1,paste,collapse="_"),
                    c('bf_all','af_all','bf_big','af_big','alpha')),1,paste,collapse='_'),
  apply(expand.grid(c('extr_amp','extr_growth'),c('bf_all','af_all','bf_big','af_big','alpha')),1,paste,
        collapse='_'))
colnames = colnames[order(colnames)]

alpha.data = as.data.frame(matrix(NA,nrow=nrow(site_info),ncol=length(colnames),
                    dimnames = list(site_info$site,colnames)))


#- 1. calculate the optimum modifier ####

for (site in site_info$site){
  obs = BACI_obs[[site]]
  sim =  BACI_sim[[site]]
  ccn = BACI_sim_ccn[[site]]
  
  obs.age = apply(obs,2,align.trw)
  obs.dia.age = apply(obs.age,2,cumsum.j)
  obs.dia.yr = apply(obs,2,cumsum.j)
  sim.dia = apply(sim,2,cumsum)*2
  sim.mean = (rowSums(sim*ccn)/rowSums(ccn))
  
  # Make index for big trees in the observed site
  # This choice is based on the size distribution in the simulation
  choose=order(apply(obs.dia.yr,2,max,na.rm=T),
               decreasing=F)[round(ncol(obs.dia.yr)*ratio_min):round(ncol(obs.dia.yr)*ratio_max)]
  choose=choose[!is.na(choose)]
  
  obs.big = obs[,choose]
  obs.age.big = obs.age[,choose]
  obs.dia.age.big = obs.dia.age[,choose]
  obs.dia.yr.big = obs.dia.yr[,choose]
  
  obs.age.mean = apply(obs.age,1,mean,na.rm=T)
  sim.mean.cut =sim.mean[!is.na(obs.age.mean)]
  obs.age.mean = obs.age.mean[!is.na(obs.age.mean)]
  
  obs.age.big.mean = apply(obs.age.big,1,mean,na.rm=T)
  sim.big.cut  = sim[!is.na(obs.age.big.mean),target_circ]
  obs.age.big.mean = obs.age.big.mean[!is.na(obs.age.big.mean)]
  
  
  #-- 1.1 Trend: RMSE #####
  
  # Find a multiplier minimizes RMSE 
  rmse = sqrt(sum((obs.age.big.mean-sim.big.cut)^2)/length(sim.big.cut))
  for(alpha in seq(0.01,5,by=0.01)){
    rmse2 =  sqrt(sum((obs.age.big.mean-alpha*sim.big.cut)^2)/length(sim.big.cut))
    if(rmse2<rmse){
      alpha_save = alpha
      rmse = rmse2
    }
  }
  
  alpha.data[site,]$trend_rmse_bf_big <- sqrt(sum((obs.age.big.mean-sim.big.cut)^2)/length(sim.big.cut))
  alpha.data[site,]$trend_rmse_af_big <- 
    sqrt(sum((obs.age.big.mean-alpha_save*sim.big.cut)^2)/length(sim.big.cut))
  alpha.data[site,]$trend_rmse_bf_all <- sqrt(sum((obs.age.mean-sim.mean.cut)^2)/length(sim.mean))
  alpha.data[site,]$trend_rmse_af_all <- 
    sqrt(sum((obs.age.mean-alpha_save*sim.mean.cut)^2)/length(sim.mean))
  alpha.data[site,]$trend_rmse_alpha <- alpha_save
  
  
  #-- 1.2 Trend: Slope ####
  nyrs = length(obs.age.big.mean)
  nyrs.all = length(obs.age.mean)
  
  slope = summary(
    lm(obs.age.big.mean - sim.big.cut ~ c(1:nyrs)))$coefficients[2,1]
  
  # Find a modifier minimizes the slope of residuals
  # Since the slope is very small, smaller numbers and intervals were applied.
  for(alpha in seq(-0.1,0.1,by=0.0001)){
    slope2 =  summary(
      lm(obs.age.big.mean - (sim.big.cut-alpha*c(1:nyrs)) ~ c(1:nyrs)))$coefficients[2,1]
    if(abs(slope2)<abs(slope)){
      alpha_save = alpha
      slope = slope2
    }
  }
  
  alpha.data[site,]$trend_slope_bf_big <- summary(
    lm(obs.age.big.mean - sim.big.cut ~ c(1:nyrs)))$coefficients[2,1]
  alpha.data[site,]$trend_slope_af_big <- summary(
    lm(obs.age.big.mean -(sim.big.cut-alpha_save*c(1:nyrs)) ~ c(1:nyrs)))$coefficients[2,1]
  alpha.data[site,]$trend_slope_bf_all <- summary(
    lm(obs.age.mean - sim.mean.cut ~ c(1:nyrs.all)))$coefficients[2,1]
  alpha.data[site,]$trend_slope_af_all <- summary(
    lm(obs.age.mean - (sim.mean.cut-alpha_save*c(1:nyrs.all)) ~
         c(1:nyrs.all)))$coefficients[2,1]
  alpha.data[site,]$trend_slope_alpha <- alpha_save
  
  
  #-- 1.3 Mature: RMSE ####
  rm(rmse,rmse2,slope,slope2,alpha_save)
  
  obs.dia.mean = apply(obs.dia.yr,1,mean,na.rm=T) 
  sim.mean.cut = 2*cumsum(sim.mean)[!is.na(obs.dia.mean)] 
  obs.dia.mean = obs.dia.mean[!is.na(obs.dia.mean)]       
  obs.dia.big.mean = apply(obs.dia.yr.big,1,mean,na.rm=T)
  # adjust lenght of datasets in case that big-tree data is shorter than all-tree data
  sim.big.cut = 2*cumsum(sim)[(is.na(obs.dia.big.mean) | is.infinite(obs.dia.big.mean)) == F ,
                              target_circ]
  obs.dia.big.mean = obs.dia.big.mean[ (is.na(obs.dia.big.mean) | is.infinite(obs.dia.big.mean)) == F ]
  nyrs = length(obs.dia.big.mean)
  
  rmse = sqrt(sum((obs.dia.big.mean[51:nyrs]-sim.big.cut[51:nyrs])^2)/
                length(sim.big.cut[51:nyrs]))
  
  # Find a multiplier minimizes RMSE
  for(alpha in seq(0.01,5,by=0.01)){
    rmse2 =  sqrt(sum((obs.dia.big.mean[51:nyrs]-alpha*sim.big.cut[51:nyrs])^2)/
                    length(sim.big.cut[51:nyrs]))
    if(rmse2<rmse){
      alpha_save = alpha
      rmse = rmse2
    }
  }
  
  alpha.data[site,]$mature_rmse_bf_big <-
    sqrt(sum((obs.dia.big.mean[51:nyrs]-sim.big.cut[51:nyrs])^2)/length(sim.big.cut[51:nyrs]))
  alpha.data[site,]$mature_rmse_af_big <- 
    sqrt(sum((obs.dia.big.mean[51:nyrs]-alpha_save*sim.big.cut[51:nyrs])^2)/length(sim.big.cut[51:nyrs]))
  alpha.data[site,]$mature_rmse_bf_all <- 
    sqrt(sum((obs.dia.mean[51:nyrs]-sim.mean.cut[51:nyrs])^2)/length(sim.mean.cut[51:nyrs]))
  alpha.data[site,]$mature_rmse_af_all <-
    sqrt(sum((obs.dia.mean[51:nyrs]-alpha_save*sim.mean.cut[51:nyrs])^2)/length(sim.mean.cut[51:nyrs]))
  alpha.data[site,]$mature_rmse_alpha <- alpha_save
  
  #-- 1.4 Mature: Slope ####
  nyrs = length(obs.dia.big.mean)
  slope = summary(
    lm(obs.dia.big.mean[51:nyrs] - sim.big.cut[51:nyrs] ~ c(51:nyrs)))$coefficients[2,1]
  
  # Find a modifier minimized slope of the residuals
  for(alpha in seq(-6,4,by=0.005)){
    
    slope2 =  summary(
      lm(obs.dia.big.mean[51:nyrs] - (sim.big.cut[51:nyrs]-alpha*c(1:(nyrs-50))) ~ 
           c(1:(nyrs-50))))$coefficients[2,1]
    if(abs(slope2)<abs(slope)){
      alpha_save = alpha
      slope = slope2
    }
  }
  
  sim.adjust = apply(sim,2,cumsum)*2-alpha_save*c(1:nrow(sim)) # calculate adjusted simulation
  sim.adjust = sim.adjust - sim.adjust[1,1]
  alpha.data[site,]$mature_slope_bf_big <- summary(
    lm(obs.dia.big.mean[51:nyrs] - sim.big.cut[51:nyrs] ~ c(51:nyrs)))$coefficients[2,1]
  alpha.data[site,]$mature_slope_af_big <- summary(
    lm(obs.dia.big.mean[51:nyrs] - (sim.big.cut[51:nyrs]-alpha_save*c(1:(nyrs-50))) ~ 
         c(1:(nyrs-50))))$coefficients[2,1]
  alpha.data[site,]$mature_slope_bf_all <- summary(
    lm(obs.dia.mean[51:nyrs] - sim.mean.cut[51:nyrs] ~ c(51:nyrs)))$coefficients[2,1]
  
  sim.adjust.mean = 2*cumsum(sim.mean)-alpha_save*c(1:nrow(sim)) # calculate adjusted simulation
  sim.adjust.mean = sim.adjust.mean + min(sim.adjust.mean,na.rm=T)
  alpha.data[site,]$mature_slope_af_all <- summary(
    lm(obs.dia.mean[51:nyrs] - (sim.mean.cut[51:nyrs]-alpha_save*c(1:(nyrs-50))) ~ 
         c(1:(nyrs-50))))$coefficients[2,1]
  alpha.data[site,]$mature_slope_alpha <- alpha_save
  
  #-- 1.5 Young: RMSE ####
  rm(rmse,rmse2,slope,slope2,alpha_save)
  
  # Cut the simulation and observation
  obs.dia.mean = apply(obs.dia.age,1,mean,na.rm=T)[1:age_young]
  sim.mean.cut = 2*cumsum(sim.mean)[1:age_young]
  obs.dia.max = apply(obs.dia.age.big,1,max,na.rm=T)[1:age_young]
  sim.big.cut = 2*cumsum(sim[1:age_young,target_circ])
  
  rmse = sqrt(sum((obs.dia.max-sim.big.cut)^2)/length(sim.big.cut))
  
  # Find a multiplier minimized RMSE
  for(alpha in seq(0.01,5,by=0.005)){
    rmse2 =  sqrt(sum((obs.dia.max-alpha*sim.big.cut)^2)/length(sim.big.cut))
    if(rmse2<rmse){
      alpha_save = alpha
      rmse = rmse2
    }
  }
  
  alpha.data[site,]$young_rmse_bf_big <-
    sqrt(sum((obs.dia.max-sim.big.cut)^2)/length(sim.big.cut))
  alpha.data[site,]$young_rmse_af_big <-
    sqrt(sum((obs.dia.max-alpha_save*sim.big.cut)^2)/length(sim.big.cut))
  alpha.data[site,]$young_rmse_bf_all <- 
    sqrt(sum((obs.dia.mean-sim.mean.cut)^2)/length(sim.mean.cut))
  alpha.data[site,]$young_rmse_af_all <- 
    sqrt(sum((obs.dia.mean-alpha_save*sim.mean.cut)^2)/length(sim.mean.cut))
  alpha.data[site,]$young_rmse_alpha <- alpha_save
  
  #-- 1.6 Young: Slope ####
  slope = summary(
    lm(obs.dia.max - sim.big.cut ~ c(1:age_young)))$coefficients[2,1]
  for(alpha in seq(-6,4,by=0.005)){
    slope2 =  summary(
      lm(obs.dia.max - (sim.big.cut-alpha*c(1:age_young)) ~ c(1:age_young)))$coefficients[2,1]
    if(abs(slope2)<abs(slope)){
      alpha_save = alpha
      slope = slope2
    }
  }
  
  sim.adjust = apply(sim,2,cumsum)*2-alpha_save*c(1:nrow(sim))
  sim.adjust = sim.adjust - sim.adjust[1,1]
  sim.adjust.mean = 2*cumsum(sim.mean)-alpha_save*c(1:nrow(sim))
  sim.adjust.mean = sim.adjust.mean + min(sim.adjust.mean,na.rm=T)
  
  alpha.data[site,]$young_slope_bf_big <- summary(
    lm(obs.dia.max - sim.big.cut ~ c(1:age_young)))$coefficients[2,1]
  alpha.data[site,]$young_slope_af_big <- summary(
    lm(obs.dia.max - (sim.big.cut-alpha_save*c(1:age_young)) ~ c(1:age_young)))$coefficients[2,1]
  alpha.data[site,]$young_slope_bf_all <- summary(
    lm(obs.dia.mean - sim.mean.cut ~ c(1:age_young)))$coefficients[2,1]
  alpha.data[site,]$young_slope_af_all <- summary(
    lm(obs.dia.mean - (sim.mean.cut-alpha_save*c(1:age_young)) ~ c(1:age_young)))$coefficients[2,1]
  alpha.data[site,]$young_slope_alpha <- alpha_save  
  
  #-- 1.7 Extreme: Amplitude ####
  # Calculation for big-trees
  # pick long-lived trees and cut the data and simulation after 1950
  idx.age = !is.na(obs.big[rownames(obs.big)==1950,]) & 
    apply(obs.big,2,function(x) sum(!is.na(x))) >50
  obs.long.big = obs.big[,idx.age]
  
  obs.cut.big=obs.long.big[which(as.numeric(rownames(obs.long.big))>1950),]
  obs.cut.big.mean=rowMeans(obs.cut.big,na.rm=T)
  
  sim.cut=sim[which(as.numeric(rownames(obs.long.big))>1950),]
  sim.cut.big = sim.cut[,5]
  
  ind.obs.low = order(obs.cut.big.mean)[1:floor(p_extr*length(obs.cut.big.mean))]
  ind.obs.high= rev(order(obs.cut.big.mean,decreasing=T)[1:floor(p_extr*length(obs.cut.big.mean))])
  ind.sim.low = order(sim.cut.big)[1:floor(p_extr*length(sim.cut.big))]
  ind.sim.high= rev(order(sim.cut.big,decreasing=T)[1:floor(p_extr*length(sim.cut.big))])
  
  obs.rank.big = obs.cut.big.mean[c(ind.obs.low,ind.obs.high)]
  sim.rank.big = sim.cut.big[c(ind.sim.low,ind.sim.high)]
  
  rmse = sqrt(sum(((sim.rank.big-mean(sim.rank.big))-(obs.rank.big-mean(obs.rank.big)))^2)/
                length(sim.rank.big))
  rmse.bf = rmse
  for(alpha in seq(0.01,5,by=0.005)){
    rmse2 =  sqrt(sum(((alpha*sim.rank.big-mean(alpha*sim.rank.big))-
                         (obs.rank.big-mean(obs.rank.big)))^2)/length(sim.rank.big))
    if(rmse2<rmse){
      alpha_save = alpha
      rmse = rmse2
    }
  }
  
  alpha.data[site,]$extr_amp_bf_big = rmse.bf
  alpha.data[site,]$extr_amp_af_big = rmse
  alpha.data[site,]$extr_amp_alpha = alpha_save

  # Calculation for all-trees
  idx.age = !is.na(obs[rownames(obs)==1950,]) & apply(obs,2,function(x) sum(!is.na(x))) >50
  obs.long = obs[,idx.age]

  obs.cut.all=obs.long[which(as.numeric(rownames(obs.long))>1950),]
  obs.cut.mean=rowMeans(obs.cut.all,na.rm=T)
  sim.cut = sim[which(as.numeric(rownames(obs.long))>1950),]
  sim.cut.mean=sim.mean[which(as.numeric(rownames(obs.long))>1950)]
  
  ind.obs.low = order(obs.cut.mean)[1:floor(p_extr*length(obs.cut.mean))]
  ind.obs.high= rev(order(obs.cut.mean,decreasing=T)[1:floor(p_extr*length(obs.cut.mean))])
  ind.sim.low = order(sim.cut.mean)[1:floor(p_extr*length(sim.cut.mean))]
  ind.sim.high= rev(order(sim.cut.mean,decreasing=T)[1:floor(p_extr*length(sim.cut.mean))])
  
  obs.rank = obs.cut.mean[c(ind.obs.low,ind.obs.high)]
  sim.rank = sim.cut.mean[c(ind.sim.low,ind.sim.high)]
  
  alpha.data[site,]$extr_amp_bf_all = sqrt(sum(((sim.rank-mean(sim.rank))-
                                                  (obs.rank-mean(obs.rank)))^2)/length(sim.rank))
  alpha.data[site,]$extr_amp_af_all = sqrt(sum(((alpha_save*sim.rank-mean(alpha_save*sim.rank))-
                                                  (obs.rank-mean(obs.rank)))^2)/length(sim.rank))
  
  #-- 1.8 Extreme: Extreme growth ####
  # Calculation for big-trees
  ind.obs.low.big = order(obs.cut.big.mean)[1:floor(p_extr*length(obs.cut.big.mean))]
  ind.obs.high.big = rev(order(obs.cut.big.mean,decreasing=T)[1:floor(p_extr*length(obs.cut.big.mean))])
  ind.sim.low.big = order(sim.cut.big)[1:floor(p_extr*length(sim.cut.big))]
  ind.sim.high.big = rev(order(sim.cut.big,decreasing=T)[1:floor(p_extr*length(sim.cut.big))])
  
  obs.extr.big = scale(obs.cut.big.mean)[c(ind.obs.low.big,ind.obs.high.big),1]
  sim.extr.big = scale(sim.cut.big)[c(ind.obs.low.big,ind.obs.high.big),1]
  rmse.bf = sqrt(sum((sim.extr.big-obs.rank.big)^2)/length(sim.extr.big))
  
  sim.cut.reorder = sim.cut
  sim.cut.reorder[order(obs.cut.big.mean),] <- sim.cut[order(sim.cut.big),]
  sim.cut.big.reorder = sim.cut.reorder[,5]
  
  sim.extr.big.reorder = scale(sim.cut.big.reorder)[c(ind.obs.low.big,ind.obs.high.big),1]
  rmse.af = sqrt(sum((sim.extr.big.reorder-obs.extr.big)^2)/length(obs.extr.big))
  
  alpha.data[site,]$extr_growth_bf_big = rmse.bf
  alpha.data[site,]$extr_growth_af_big = rmse.af
  
  # Calculation for all-trees
  ind.obs.low = order(obs.cut.mean)[1:floor(p_extr*length(obs.cut.mean))]
  ind.obs.high= rev(order(obs.cut.mean,decreasing=T)[1:floor(p_extr*length(obs.cut.mean))])
  ind.sim.low = order(sim.cut.mean)[1:floor(p_extr*length(sim.cut.mean))]
  ind.sim.high= rev(order(sim.cut.mean,decreasing=T)[1:floor(p_extr*length(sim.cut.mean))])
  
  obs.extr = scale(obs.cut.mean)[c(ind.obs.low,ind.obs.high),1]
  sim.extr = scale(sim.cut.mean)[c(ind.obs.low,ind.obs.high),1]
  
  rmse.bf = sqrt(sum((sim.extr-obs.extr)^2)/length(sim.extr))
  
  sim.cut.mean.reorder = sim.cut.mean
  sim.cut.mean.reorder[order(obs.cut.big.mean)] <- sim.cut.mean[order(sim.cut.big)]
  
  sim.extr.reorder = scale(sim.cut.mean.reorder)[c(ind.obs.low,ind.obs.high),1]
  rmse.af = sqrt(sum((sim.extr.reorder-obs.extr)^2)/length(obs.extr))
  
  alpha.data[site,]$extr_growth_bf_all = rmse.bf
  alpha.data[site,]$extr_growth_af_all = rmse.af

}

improve.data = data.frame(
  site = site_info$site,
  trend_rmse_diff_all = alpha.data$trend_rmse_bf_all - alpha.data$trend_rmse_af_all,
  trend_slope_diff_all = abs(alpha.data$trend_slope_bf_all) - abs(alpha.data$trend_slope_af_all),
  mature_rmse_diff_all = alpha.data$mature_rmse_bf_all - alpha.data$mature_rmse_af_all,
  mature_slope_diff_all = abs(alpha.data$mature_slope_bf_all) - abs(alpha.data$mature_slope_af_all),
  young_rmse_diff_all = alpha.data$young_rmse_bf_all - alpha.data$young_rmse_af_all,
  young_slope_diff_all = abs(alpha.data$young_slope_bf_all) - abs(alpha.data$young_slope_af_all),
  extr_amp_diff_all = alpha.data$extr_amp_bf_all - alpha.data$extr_amp_af_all,
  extr_growth_diff_all = alpha.data$extr_growth_bf_all - alpha.data$extr_growth_af_all
)

write.csv(improve.data,file='report/Table2_data.scv')

