## have data already in correct function.
##so the input data is the following
# a list of...
# list(Y = Y, nSites = dim(Y)[2], nYears = dim(Y)[1], nseed = nseed, td = td, valley = valley, se.mice = se.mice)

#y = a matrix

# ################################################################################
# Gompertz model with same parameters for each grid
# Test for a difference between valleys in b

valley <- c(rep(1, 4), rep(2, 4))

mod <- "model {

   # Observation model
   # The Rs are different in each site
   for(i in 1:nYears) {
     for(j in 1:nSites) {
       prec[i, j] <- 1 / (se.mice[i, j] * se.mice[i, j])
       Y[i, j] ~ dnorm(X[i, j], prec[i, j])
       nseed[i, j] ~ dnorm(S[i], tauS)     # model missing seed values as the mean for that year
       se.mice[i, j] ~ dnorm(mu.se, tau.se)T(0, )
     }
     S[i] ~ dnorm(0, 0.01)
   }

   # Process model
    for(j in 1:nSites) {

      X[1, j] ~ dnorm(0, 0.01) # vague normal prior
      simX[1, j] <- X[1, j]       # the fitted model

      for(i in 2:nYears) {

        X[i, j] ~ dnorm(predX[i, j], tauQ)
        predX[i, j] <- X[i-1, j] + (b[valley[j]] + beta.seed * nseed[i, j] - d * X[i-1, j]) * td[i]

        simX[i, j] ~ dnorm(predsimX[i, j], tauQ)    # fitted model
        predsimX[i, j] <- simX[i- 1, j] + (b[valley[j]] + beta.seed * nseed[i, j] - d * simX[i-1, j]) * td[i]
      }
   }


   b[1] ~ dnorm(0, 0.01)T(0, )
   b[2] ~ dnorm(0, 0.01)T(0, )
   d ~ dnorm(0, 0.01)T(0, )
   beta.seed ~ dnorm(0, 0.01)
   tauQ <- 1 / (sigmaQ * sigmaQ)
   sigmaQ ~ dunif(0, 10)
   tauS <- 1 / (sigmaS * sigmaS)
   sigmaS ~ dunif(0, 10)

   mu.se ~ dnorm(0, 0.001)
   tau.se <- 1 / (sigma.se * sigma.se)
   sigma.se ~ dunif(0, 10)

   diff.b <- b[1] - b[2]

  }"


# write model
write(mod, "model.txt")

# run the model in JAGS with three chains
mod.jags <- jags(model = "model.txt",
                 data = list(Y = Y, nSites = dim(Y)[2], nYears = dim(Y)[1], nseed = nseed, td = td, valley = valley, se.mice = se.mice),
                 param = c("b", "diff.b", "d", "beta.seed", "X", "sigmaQ", "simX"),
                 n.chains = 3,
                 n.iter = 10000,
                 n.burnin = 5000,
                 parallel = TRUE)

mod.sum <- mod.jags$summary
mod.sum
mod.jags$DIC

x <- data.frame(mod.sum[substr(rownames(mod.sum), 1, 1) == "X", c(1, 3, 7)])
x <- exp(x)
names(x) <- c("mp", "lcl", "ucl")
x$year.dec <- rep(year.dec, dim(Y)[2])
x$grid <- rep(levels(factor(dat$grid)), each = dim(Y)[1])

simx <- data.frame(mod.sum[substr(rownames(mod.sum), 1, 4) == "simX", c(1, 3, 7)])
simx <- exp(simx)
names(simx) <- c("mp", "lcl", "ucl")
simx$year.dec <- rep(year.dec, dim(Y)[2])
simx$grid <- rep(levels(factor(dat$grid)), each = dim(Y)[1])

ggplot(dat, aes(y = exp(lN), x = year.dec)) +
  geom_point(aes(colour = grid)) +
  geom_line(aes(colour = grid)) +
  geom_line(data = x, aes(y = mp)) +
  geom_line(data = simx, aes(y = mp), colour = "red", lwd = 1.1) +
  facet_wrap(~ grid) +
  theme_bw() +
  theme(legend.position = "none")


# predicted vs observed
act <- select(dat, year.dec, grid, lN)
act <- full_join(simx, act)
act$valley <- substr(act$grid, 1, 3)

ggplot(act, aes(y = log(mp), x = lN)) +
  geom_point(size = 2) +
  facet_wrap(~ valley) +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic()


################################################################################
# Gompertz model with same dynamics on each grid but seed effects
# And difference between valleys in b

mod <- "model {

   # Observation model
   # The Rs are different in each site
   for(i in 1:nYears) {
     for(j in 1:nSites) {
       prec[i, j] <- 1 / (se.mice[i, j] * se.mice[i, j])
       Y[i, j] ~ dnorm(X[i, j], prec[i, j])

       nred[i, j] ~ dnorm(Sr[i], tauR)
       nsilver[i, j] ~ dnorm(Ss[i], tauS)
       se.mice[i, j] ~ dnorm(mu.se, tau.se)T(0, )
     }
     Sr[i] ~ dnorm(0, 0.01)
     Ss[i] ~ dnorm(0, 0.01)
   }

   # Process model
    for(j in 1:nSites) {

      X[1, j] ~ dnorm(0, 0.01) # vague normal prior
      simX[1, j] <- X[1, j]       # the fitted model

      for(i in 2:nYears) {

        X[i, j] ~ dnorm(predX[i, j], tauQ)
        predX[i, j] <- X[i-1, j] + (b[valley[j]] + beta.red * nred[i, j] + beta.silver * nsilver[i, j] - d * X[i-1, j]) * td[i]

        simX[i, j] ~ dnorm(predsimX[i, j], tauQ)    # fitted model
        predsimX[i, j] <- simX[i- 1, j] + (b[valley[j]] + beta.red * nred[i, j] + beta.silver * nsilver[i, j] - d * simX[i-1, j]) * td[i]
      }
   }


   b[1] ~ dnorm(0, 0.01)T(0, )
   b[2] ~ dnorm(0, 0.01)T(0, )
   d ~ dnorm(0, 0.01)T(0, )
   beta.red ~ dnorm(0, 0.01)
   beta.silver ~ dnorm(0, 0.01)
   tauQ <- 1 / (sigmaQ * sigmaQ)
   sigmaQ ~ dunif(0, 10)
   tauS <- 1 / (sigmaS * sigmaS)
   sigmaS ~ dunif(0, 10)
   tauR <- 1 / (sigmaR * sigmaR)
   sigmaR ~ dunif(0, 10)

   mu.se ~ dnorm(0, 0.001)
   tau.se <- 1 / (sigma.se * sigma.se)
   sigma.se ~ dunif(0, 10)

   diff.b <- b[1] - b[2]

  }"


# write model
write(mod, "model.txt")

# run the model in JAGS with three chains
mod.jags <- jags(model = "model.txt",
                 data = list(Y = Y, nSites = dim(Y)[2], nYears = dim(Y)[1], nred = nred, nsilver = nsilver, td = td,
                             valley = valley, se.mice = se.mice),
                 param = c("b", "diff.b", "d", "beta.red", "beta.silver", "X", "sigmaQ", "simX"),
                 n.chains = 3,
                 n.iter = 10000,
                 n.burnin = 5000,
                 parallel = TRUE)

mod.sum <- mod.jags$summary
mod.sum
mod.jags$DIC

x <- data.frame(mod.sum[substr(rownames(mod.sum), 1, 1) == "X", c(1, 3, 7)])
x <- exp(x)
names(x) <- c("mp", "lcl", "ucl")
x$year.dec <- rep(year.dec, dim(Y)[2])
x$grid <- rep(levels(factor(dat$grid)), each = dim(Y)[1])

simx <- data.frame(mod.sum[substr(rownames(mod.sum), 1, 4) == "simX", c(1, 3, 7)])
simx <- exp(simx)
names(simx) <- c("mp", "lcl", "ucl")
simx$year.dec <- rep(year.dec, dim(Y)[2])
simx$grid <- rep(levels(factor(dat$grid)), each = dim(Y)[1])

ggplot(dat, aes(y = exp(lN), x = year.dec)) +
  geom_point(aes(colour = grid)) +
  geom_line(aes(colour = grid)) +
  geom_line(data = x, aes(y = mp)) +
  geom_line(data = simx, aes(y = mp), colour = "red", lwd = 1.1) +
  facet_wrap(~ grid) +
  theme_bw() +
  theme(legend.position = "none")


# predicted vs observed
act <- select(dat, year.dec, grid, lN)
act <- full_join(simx, act)
act$valley <- substr(act$grid, 1, 3)

ggplot(act, aes(y = log(mp), x = lN)) +
  geom_point(size = 2) +
  facet_wrap(~ valley) +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic()


################################################################################
# Gompertz model with same dynamics on each grid but seed effects
# And difference between valleys in b
# Include rats

mod <- "model {

   # Observation model
   # The Rs are different in each site
   for(i in 1:nYears) {
     for(j in 1:nSites) {
       prec[i, j] <- 1 / (se.mice[i, j] * se.mice[i, j])
       Y[i, j] ~ dnorm(X[i, j], prec[i, j])

       nred[i, j] ~ dnorm(Sr[i], tauR)
       nsilver[i, j] ~ dnorm(Ss[i], tauS)
       rat.mna[i, j] ~ dpois(mu.rat)
       se.mice[i, j] ~ dnorm(mu.se, tau.se)T(0, )
     }
     Sr[i] ~ dnorm(0, 0.01)
     Ss[i] ~ dnorm(0, 0.01)
   }

   # Process model
    for(j in 1:nSites) {

      X[1, j] ~ dnorm(0, 0.01) # vague normal prior
      simX[1, j] <- X[1, j]       # the fitted model

      for(i in 2:nYears) {

        X[i, j] ~ dnorm(predX[i, j], tauQ)
        predX[i, j] <- X[i-1, j] + (b[valley[j]] + beta.red * nred[i, j] + beta.silver * nsilver[i, j] + beta.rat * rat.mna[i, j] - d * X[i-1, j]) * td[i]

        simX[i, j] ~ dnorm(predsimX[i, j], tauQ)    # fitted model
        predsimX[i, j] <- simX[i- 1, j] + (b[valley[j]] + beta.red * nred[i, j] + beta.silver * nsilver[i, j] + beta.rat * rat.mna[i, j] - d * simX[i-1, j]) * td[i]
      }
   }


   b[1] ~ dnorm(0, 0.01)T(0, )
   b[2] ~ dnorm(0, 0.01)T(0, )
   d ~ dnorm(0, 0.01)T(0, )
   beta.red ~ dnorm(0, 0.01)
   beta.silver ~ dnorm(0, 0.01)
   beta.rat ~ dnorm(0, 0.01)
   mu.rat ~ dgamma(0.01, 0.01)

   tauQ <- 1 / (sigmaQ * sigmaQ)
   sigmaQ ~ dunif(0, 10)
   tauS <- 1 / (sigmaS * sigmaS)
   sigmaS ~ dunif(0, 10)
   tauR <- 1 / (sigmaR * sigmaR)
   sigmaR ~ dunif(0, 10)

   mu.se ~ dnorm(0, 0.001)
   tau.se <- 1 / (sigma.se * sigma.se)
   sigma.se ~ dunif(0, 10)

   diff.b <- b[1] - b[2]

  }"


# write model
write(mod, "model.txt")

# run the model in JAGS with three chains
mod.jags <- jags(model = "model.txt",
                 data = list(Y = Y, nSites = dim(Y)[2], nYears = dim(Y)[1], nred = nred, nsilver = nsilver, td = td,
                             valley = valley, rat.mna = rat.mna, se.mice = se.mice),
                 param = c("b", "diff.b", "d", "beta.red", "beta.silver", "beta.rat", "X", "sigmaQ", "simX"),
                 n.chains = 3,
                 n.iter = 10000,
                 n.burnin = 5000,
                 parallel = TRUE)

mod.sum <- mod.jags$summary
mod.sum
mod.jags$DIC

x <- data.frame(mod.sum[substr(rownames(mod.sum), 1, 1) == "X", c(1, 3, 7)])
x <- exp(x)
names(x) <- c("mp", "lcl", "ucl")
x$year.dec <- rep(year.dec, dim(Y)[2])
x$grid <- rep(levels(factor(dat$grid)), each = dim(Y)[1])

simx <- data.frame(mod.sum[substr(rownames(mod.sum), 1, 4) == "simX", c(1, 3, 7)])
simx <- exp(simx)
names(simx) <- c("mp", "lcl", "ucl")
simx$year.dec <- rep(year.dec, dim(Y)[2])
simx$grid <- rep(levels(factor(dat$grid)), each = dim(Y)[1])
