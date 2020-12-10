# This code is written to
# build observation RData using selected 11 sites for further convenience.
# The uploaded directory contains the finalized simulation.RData, however, 
# to build the observation list (#-2), you need to download the BACI dataset.
# You can download it by registration with your email on the site :
# http://www.baci-h2020.eu/

#- 0. Initialize ####

library(ncdf4)
source('R/fun_cal_trw.R')
BACI_sim.r5698 = list()
BACI_sim.ccn.r5698 = list()

site_info = data.frame(site=c('DEO','DVN','GIU','HD2','SCH','SOB','TIC','CAN','SOR','TER','ZOF'),
                       specie = c('PCAB','PCAB','PCAB','PCAB','PCAB','PCAB','PCAB',
                                  'FASY','FASY','FASY','FASY'),
                       PFT=c(4,4,4,4,4,4,4,6,6,6,6))

file.path='data/201961114214112_treeringdata_BACI2016_WP3/data/treeringbiomass_network_Europe/raw_data/'

#- 1. Build observation list ####

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

#- 2. Save ####
save(BACI_obs, file='data/BACI_obs.RData')

