## Orginal example dataset
# source(./R)

## Simulate dataset
# list(Y = lY, nGrids = nGrids, nYears = nYears, res = res, lag = 1)
## Y = lY
##nGrids = nGrids
##nYears = nYears,
## res = res
## lag = 1

## Y is a data structure that is ???


##nGrids


##nYears


## res


## lag
#1 would make biological sense
################################################################################
# Gompertz model
# here K = exp(b/d)

# specify the avaliable food parameter lag in the data statement

mod <- "model {

   # Observation model
   for(i in 1:nYears) {
     for(j in 1:nGrids) {
       Y[i, j] ~ dnorm(X[i], tau[1])
     }
   }

   X[1] ~ dnorm(0, 1)
   simX[1] <- X[1]

   # Process model
   for(i in 2:nYears) {

       X[i] ~ dnorm(predX[i], tau[2])
       predX[i] <- X[i-1] + b + beta.res * res[i-lag] - d * X[i-1]

       simX[i] ~ dnorm(predsimX[i], tau[2])    # fitted model
       predsimX[i] <- simX[i-1] + b + beta.res * res[i-lag] - d * simX[i-1]
   }

   b ~ dnorm(0, 0.001)T(0, )
   beta.res ~ dnorm(0, 0.01)
   d ~ dnorm(0, 0.01)T(0, )
   mu.r ~ dnorm(0, 0.001)

   for(i in 1:2) {
     tau[i] <- 1 / (sigma[i] * sigma[i])
     sigma[i] ~ dunif(0, 10)
   }

  }"


# write model
write(mod, "model.txt")

# run the model in JAGS with three chains
mod.jags <- jags(model = "model.txt",
                 data = list(Y = lY, nGrids = nGrids, nYears = nYears, res = res, lag = 1),
                 param = c("b", "d", "beta.res", "sigma", "X", "simX"),
                 n.chains = 3,
                 n.iter = 10000,
                 n.burnin = 5000,
                 parallel = TRUE)

mod.sum2 <- mod.jags$summary
mod.sum2
mod.jags
