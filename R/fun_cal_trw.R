# Jina Jeong
# Documentation added 25.06.2020
# ver 2. 03.07.2020 : extrac_ccn has added 

# cal_trw
#- This function is written to calculate TRW from CCDELTABA (increment of basal area)
#- Please note that this function work only for PFT_def with no age
#- init : Consider initial basal area. If init=F, TRW is re-formated from 0. 
# extract_ccn
#- This function is written to extract CCN (number of individuals per size classes) 
#- which is needed to calculate mean TRW

library(ncdf4)

cal_trw <- function(site_name,config,PFT=4,fourd=T,init=F,ncirc = 3){
  
  if (fourd==T){
    
    infile.4d = paste0('output/',site_name,'_',config,'_4dim.nc')
    PFTidx=PFT
    
    # open a netCDF file & some settings
    ncin.4d = nc_open(infile.4d)
    CCBA.mat = ncvar_get(ncin.4d, 'CCBA')[PFTidx,,]
    DELTA.mat = ncvar_get(ncin.4d, 'CCDELTABA')[PFTidx,,]
    t = ncvar_get(ncin.4d, "time_counter")
    nyears = length(t)/365
    
    treering.from.DELTA <-matrix(NA,nrow=nyears,ncol=ncirc)
    
    for (ii in 1:nrow(CCBA.mat)){
      CCBA.PFT.reshape = matrix(CCBA.mat[ii,], nrow=365) # put every year in a column
      DELTA.PFT.reshape = matrix(DELTA.mat[ii,], nrow=365) 
      
      CCBA.from.DELTA = cumsum(DELTA.mat[ii,])
      if(init){
        CCBA.from.DELTA=CCBA.from.DELTA+CCBA.mat[ii,1]
      }
      
      CCBA.from.DELTA = matrix(CCBA.from.DELTA,nrow = 365)
      R.from.DELTA=sqrt(CCBA.from.DELTA/pi)[365,]
      
      if(init){
        R.from.DELTA=R.from.DELTA-sqrt(CCBA.mat[ii,1]/pi)
      }
      
      treering.from.DELTA[,ii] = 1000*c(R.from.DELTA[1],diff(R.from.DELTA))
      
    }
  } else {
    
    ## for r5698. There is no 4dim.nc 
    infile = paste0('output/',site_name,'_',config,'.nc')
    
    ncin      = nc_open(infile)
    t = ncvar_get(ncin, "time_counter")
    nyears = length(t)/365
    PFTidx = PFT

    #Make matrix to fill
    treering.from.DELTA <-matrix(NA,nrow=nyears,ncol=ncirc)

    # Because 3d file has output for each circ_class,
    # output file needs to be get by each circ_class
    for (ii in 1:ncirc){
      
      if(ii<=9){
        var.name1 = paste0('CCBA_00',ii)          # define variable
        var.name2 = paste0('CCDELTABA_00',ii)
      } else {
        var.name1 = paste0('CCBA_0',ii)          
        var.name2 = paste0('CCDELTABA_0',ii)
      }
      
      CCBA.array = ncvar_get(ncin, var.name1);
      DELTA.array = ncvar_get(ncin, var.name2)
      
      if (length(PFTidx) > 2){
        #Calculation!
        CCBA.PFT = colSums(CCBA.array[PFTidx,]) #pick numbers out (specific pft) 
        DELTA.PFT = colSums(DELTA.array[PFTidx,])
      } else {
        CCBA.PFT = (CCBA.array[PFTidx,]) #pick numbers out (specific pft) 
        DELTA.PFT = (DELTA.array[PFTidx,])  
        
      }
      
      CCBA.PFT.reshape = matrix(CCBA.PFT, nrow=365) 
      DELTA.PFT.reshape = matrix(DELTA.PFT, nrow=365) 
      
      CCBA.from.DELTA = cumsum(DELTA.PFT) 
      CCBA.from.DELTA = matrix(CCBA.from.DELTA,nrow = 365)
      
      CCBA.d0 <-CCBA.PFT.reshape[c(1),1]

      if (init){
        CCBA.from.DELTA = CCBA.from.DELTA + CCBA.d0
        R.from.DELTA = sqrt(CCBA.from.DELTA/pi)
        TRW = c(R.from.DELTA[365,1]-sqrt(CCBA.d0/pi),diff(R.from.DELTA[365,]))*1000
      } else {
        R.from.DELTA = sqrt(CCBA.from.DELTA/pi)
        TRW = c(R.from.DELTA[365,1],diff(R.from.DELTA[365,]))*1000
      }

      treering.from.DELTA[,ii] = TRW
      
    }
    # end loop
    
    nc_close(ncin)
    treering.from.DELTA <- as.data.frame(treering.from.DELTA)
    
  }
  
  return(treering.from.DELTA)
  
}

extract_ccn <- function(site_name,config,PFT=4,fourd=T,ncirc = 3){
  if (fourd==T){
    
    infile.4d = paste0('output/',site_name,'_',config,'_4dim.nc')
    ncin      = nc_open(infile.4d)
    PFTidx=PFT
    CCAN.mat = ncvar_get(ncin,'CCN')[PFT,,]
  } else {
    
    infile = paste0('output/',site_name,'_',config,'.nc')
    
    ncin      = nc_open(infile)
    t = ncvar_get(ncin, "time_counter")
    nyears = length(t)/365
    PFTidx = PFT
    
    #Make matrix to fill
    CCN.mat <-matrix(NA,nrow=nyears,ncol=ncirc)
    
    for (icirc in 1:ncirc){
      if(icirc<=9){
        var.name = paste0('CCN_00',icirc)          # define variable
      } else {
        var.name = paste0('CCN_0',icirc)          
      }
      CCN.PFT = ncvar_get(ncin,var.name)[PFTidx,]
      CCN.mat[,icirc] <- apply(matrix(CCN.PFT,nrow=365),2,mean)
    } # for icirc
    
    
  } # if fourd
  return(CCN.mat)
}
  