## Example for the calculation of TRW from the simulated output
# This code is written to give examples to use functions in fun_cal_trw.R 
# with simulated output from ORCHIDEE r5698.
# All outputs were proceeded to output/BACI_sim.r5698.RData and output/BACI_si.ccn.r5698.RData with the
# below code, but all outputs were not uploaded due to the storage of the archive.
# Author: Jina Jeong (j.jeong@vu.nl)

#- 0. Initialize ####
source('R/fun_cal_trw.R')

load('data/BACI_obs.RData')
BACI_sim.r5698 = list()
BACI_sim.ccn.r5698 = list()

#site_info = data.frame(site=c('DEO','DVN','GIU','HD2','SCH','SOB','TIC','CAN','SOR','TER','ZOF'),
#                       specie = c('PCAB','PCAB','PCAB','PCAB','PCAB','PCAB','PCAB',
#                                  'FASY','FASY','FASY','FASY'),
#                       PFT=c(4,4,4,4,4,4,4,6,6,6,6))

# There is only output for the site DEO in the archive.
site_info = data.frame(site='DEO',species='PCAB',PFT=4)

#- 1. Calculate tree-ring width and save ####
for (site in site_info$site){
  ipft = site_info$PFT[which(site_info$site==site)]
  
  # fourd should be F for r5698 output
  trw.sim <- cal_trw(site_name=paste0('BACI.',site),config='r5698',PFT=ipft,fourd=F,init=F,ncirc=5)
  ccn <- extract_ccn(site_name=paste0('BACI.',site),config='r5698',PFT=ipft,fourd=F,ncirc=5)
  
  # The observation from site CAN is one year longer than currently available years from the forcing
  # (CRU-NCEP v.5.3.2) that is used for the simulation, therefore, row names should be cut.
  if (site == 'CAN'){
    rownames(trw.sim) <- rownames(BACI_obs[[site]])[1:145]
    rownames(ccn) <- rownames(BACI_obs[[site]])[1:145]
    
  } else {
    rownames(trw.sim) <- rownames(BACI_obs[[site]])
    rownames(ccn) <- rownames(BACI_obs[[site]])
  }
  
  BACI_sim.r5698[[site]] <- trw.sim
  BACI_sim.ccn.r5698[[site]] <- ccn
  
  rm(trw.sim,ccn)
}

#- 2. Save calculated output ####

#save(BACI_sim.r5698,file='output/BACI_sim.r5698.RData')
#save(BACI_sim.ccn.r5698,file='output/BACI_sim.ccn.r5698.RData')



