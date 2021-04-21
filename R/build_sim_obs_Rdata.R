# This code is written to
# 1) build observation RData using selected 11 sites
# 2) build simulation RData by calculating annual TRW from daily output
# for further convienience
# The uploaded directory contains finalized simulation.RData, however, 
# to build the observation list (#-2), you need to download the BACI dataset.
# You can download it by registration with your email on the site :
# http://www.baci-h2020.eu/
# All of the simulated output was not uploaded due to its large size, one site,
# DEO was uploaded in ouput directory to give an example.

#- 1. Initialize ####

library(ncdf4)
source('R/fun_cal_trw.R')
BACI_sim.r5698 = list()
BACI_sim.ccn.r5698 = list()

site_info = data.frame(site=c('DEO','DVN','GIU','HD2','SCH','SOB','TIC','CAN','SOR','TER','ZOF'),
                       specie = c('PCAB','PCAB','PCAB','PCAB','PCAB','PCAB','PCAB',
                                  'FASY','FASY','FASY','FASY'),
                       PFT=c(4,4,4,4,4,4,4,6,6,6,6))

file.path='data/201961114214112_treeringdata_BACI2016_WP3/data/treeringbiomass_network_Europe/raw_data/'

#- 2. Build observation list ####

BACI_obs = list()
for (site in site_info$site){
  trw = read.table(paste0(file.path,site,'.txt'),row.names=1,header=T)
  trw = as.matrix(trw)
  
  # The result for the last year of the site 'CAN' was cut since atmospheric forcing
  # utilized for this research (CRU-NCEP v.5.3.2) is ended 2012 
  if(site=='CAN'){
    BACI_obs[[site]] <- trw[1:145,]
    
  }else { 
    BACI_obs[[site]] <- trw}

}

save(BACI_obs, file='data/BACI_obs.RData')

#- 3. Build simulation list #####

load('data/BACI_obs.RData')

#-- 3.1 Calculate tree-ring width ####
for (site in site_info$site){
  ipft = site_info$PFT[which(site_info$site==site)]
  
  # fourd should be F for r5698 output
  trw.sim <- cal_trw(site_name=paste0('BACI.',site),config='r5698',PFT=ipft,fourd=F,init=F,ncirc=5)
  ccn <- extract_ccn(site_name=paste0('BACI.',site),config='r5698',PFT=ipft,fourd=F,ncirc=5)
 
  rownames(trw.sim) <- rownames(BACI_obs[[site]])
  rownames(ccn) <- rownames(BACI_obs[[site]])
  
  BACI_sim.r5698[[site]] <- trw.sim
  BACI_sim.ccn.r5698[[site]] <- ccn
  
  rm(trw.sim,ccn)
}


#-- 3.2 Save calculated output ####
save(BACI_sim.r5698,file='output/BACI_sim.r5698.RData')
save(BACI_sim.ccn.r5698,file='output/BACI_sim.ccn.r5698.RData')


