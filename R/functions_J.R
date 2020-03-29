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

error_margin <- function(vector, interval) {
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
  return(error)
}
