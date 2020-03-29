# Function to calculate tree-ring width from ORCHIDEE-output.
# site : target_site
# config : configuration
# fourd : if output is 4-dimentional or not. r5669 doesn't write 4-dim output.
# init : if consider initional diameter or not. 

# Author: Jina Jeong (j.jeong@vu.nl)


cal_trw <- function(site_name,config,PFT=4,fourd=T,init=F){
  
  if (fourd==T){
    
    infile.4d = paste0('output/',site_name,'_',config,'_4dim.nc')
    PFTidx=PFT
    
    # open a netCDF file & some settings
    ncin.4d = nc_open(infile.4d)
    CCBA.mat = ncvar_get(ncin.4d, 'CCBA')[PFTidx,,]
    DELTA.mat = ncvar_get(ncin.4d, 'CCDELTABA')[PFTidx,,]
    TRW.mat = ncvar_get(ncin.4d, 'CCTRW')[PFTidx,,]
    t = ncvar_get(ncin.4d, "time_counter")
    nyears = length(t)/365
    ncirc = 3
    
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
    
    infile = paste0('output/',site_name,'_',config,'.nc')
    
    ncin      = nc_open(infile)
    t = ncvar_get(ncin, "time_counter")
    nyears = length(t)/365
    ncirc = 3
    
    PFT_len=length(grep('PFT',names(ncatt_get(ncin,varid=0))))
    if(PFT_len>14){
      PFTidx=c(10:13)
    } else {
      PFTidx= PFT
    }
    
    #Make matrix to fill
    treering.from.DELTA <-matrix(,nrow=nyears,ncol=ncirc)

    # Loop : circ level
    for (ii in 1:ncirc){
      #ii = 1
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
        CCBA.PFT = colSums(CCBA.array[PFTidx,1:(nyears*365)]) #pick numbers out (specific pft) 
        DELTA.PFT = colSums(DELTA.array[PFTidx,1:(nyears*365)])
      } else {
        CCBA.PFT = (CCBA.array[PFTidx,1:(nyears*365)]) #pick numbers out (specific pft) 
        DELTA.PFT = (DELTA.array[PFTidx,1:(nyears*365)])  
        
      }
      
      CCBA.PFT.reshape = matrix(CCBA.PFT, nrow=365) # put every year in a column
      DELTA.PFT.reshape = matrix(DELTA.PFT, nrow=365) 
      
      CCBA.from.DELTA = cumsum(DELTA.PFT) 
      CCBA.from.DELTA = matrix(CCBA.from.DELTA,nrow = 365)
      
      CCBA.d0 <-matrix(CCBA.PFT.reshape[c(1),]-DELTA.PFT.reshape[c(1),])

      CCBA.from.DELTA.d0 = CCBA.from.DELTA[c(1),]-DELTA.PFT.reshape[c(1),]
      R.d0.from.DELTA = (CCBA.from.DELTA.d0/pi)^(0.5) 
      R.d365.from.DELTA = (CCBA.from.DELTA[365,]/pi)^(0.5)
      
      treering.from.DELTA[,c(ii)] = 1000*(R.d365.from.DELTA-R.d0.from.DELTA)
      
    }
    # end loop
    
    #close.ncdf(ncin)
    treering.from.DELTA <- as.data.frame(treering.from.DELTA)
    
    
  }
  
  return(treering.from.DELTA)
  
}