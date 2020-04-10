### Model equations for population growth
## Gompertz model
## ANthony Davidson
## April2020

## Function
gompertz_model_fun <- function(yr = yr, N = N, r = r, d = d){
  
  for(i in 1:length(yr)){
    
    N[i + 1] <- N[i] * exp((r - d * log(N[i]))) # head(N)
    
  } #i
  
  sim.dat1 <- data.frame(cbind(N[c(1:length(yr))], yr[c(1:length(yr))]))
  names(sim.dat1) <- c("est.N", "year")
  
  # plot
  # par(mfrow=c(1,3))
  return(
    plot(est.N~year, 
         type = "l", 
         lty = "dashed", 
         lwd = 3, 
         col = "red", 
         data = sim.dat1,
         main = "Gops model")
  )
  
}

### Example run
#data storage matrixs for plotting
# yr <-seq(1,50,1)
# N <- matrix(NA,nrow=length(yr),ncol=1)
# 
# length(N)
# 
# #parameters
# r <- 0.5
# d <- 0.15
# N[1] <- 1

# gompertz_model_fun(yr = yr, N = N, r = r, d = d)
# yr <-seq(1,50,1)
# N <- matrix(NA,nrow=length(yr),ncol=1)
# length(N)
#parameters
# r <- 0.5
# d <- 0.15
# N[1] <- 1
# N = matrix(NA,nrow=length(yr),ncol=1)
# yr = seq(1,50,1)
# gompertz_model_fun(yr = yr, N = N, r = r, d = d)


# mapping
#sure purrr to deal with this...
# April2020