# Small functions to process 4 tree-ring benchmarks.
# Author: Jina Jeong (j.jeong@vu.nl)


# function to align .rwl dataset by age
align.trw <- function(xx){
  xy = rep(NA,length(xx))
  xy[1:length(na.omit(xx))] = as.numeric(na.omit(xx))
  return(xy)
}

# function to .rwl ring width datset to diameter
cumsum.j <- function(x){
  x[which(is.na(x)==F)]=2*cumsum(x[which(is.na(x)==F)])
  return(x)
}

