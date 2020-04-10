
  library(tidyverse)
  library(jagsUI)
  
  setwd("c:\\users\\s429217\\onedrive\\data\\mpd\\state space\\data")
  
############################################################################
## Read in rat, mouse and possum abundance data on log scale

  rat.pos <- read.csv("c:\\users\\s429217\\onedrive\\data\\mpd\\state space\\data\\Rat possum mouse log N estimates.csv")
  rat.pos <- dplyr:::select(rat.pos, -X)
  glimpse(rat.pos)                                                                                  
  
# choose species
  rat.pos$xx <- rat.pos$rat_abun
  rat.pos$xxse <- rat.pos$rat_sd

# choose treatments
# exclude "Possum + rat removal" treatment when analysing rats
# exclude the possum removal treatments when analysing possums

  rat.pos <- filter(rat.pos, treat %in% c("Possum removal", "Control", "Stoat removal"))
  
# exclude first two trips
  rat.pos <- filter(rat.pos, month > 2)
  
  
  ggplot(rat.pos, aes(y = xx, x = month)) +
    geom_point() +
    geom_line(data = filter(rat.pos, !is.na(xx))) +
    geom_errorbar(aes(ymin = xx - xxse, ymax = xx + xxse)) +
    facet_wrap(~ grid) +
    theme_bw()                                                          
  
############################################################################
## Read in seed trap data
  trap <- read.csv("c:\\users\\s429217\\onedrive\\data\\mpd\\state space\\data\\seed trap data.csv")
  trap <- dplyr:::select(trap, -X)
  glimpse(trap)
  
  trap$month <- ifelse(trap$grid == "Okopeka East" & is.na(trap$month) == T, 31, trap$month)
  trap$month <- ifelse(trap$grid == "Okopeka West" & is.na(trap$month) == T, 31, trap$month)
  
## Summarise by seed trap, summing seed dry weights
  seed.trap <- group_by(trap, seed.trap, grid, month) %>%
               summarise(seed_mass = sum(dry.wt, na.rm=T),
               seed_mass_poss_pref = sum(dry.wt[spp %in% c("Beitaw", "Eladen", "Rubcis", "Melram", "Hedarb", "Pencor")], na.rm=T),
               seed_mass_poss_avoid = sum(dry.wt[spp %in% c("Ripsca", "Prufer")], na.rm=T))

# convert to log scale
  seed.trap <- seed.trap %>%
               mutate(seed_mass = log(seed_mass + 0.001),
                      seed_mass_poss_pref = log(seed_mass_poss_pref + 0.001),
                      seed_mass_poss_avoid = log(seed_mass_poss_avoid + 0.001))
                     
# now summarise by grid and month calculating mean seed trap values
  mean.seed <- group_by(seed.trap, grid, month) %>%
               summarise(seed_mass = mean(seed_mass, na.rm=T),
               seed_mass_poss_pref = mean(seed_mass_poss_pref, na.rm=T),
               seed_mass_poss_avoid = mean(seed_mass_poss_avoid, na.rm=T)) 
               
  ggplot(mean.seed, aes(y = seed_mass, x = month)) +
    geom_point() +
    geom_line(data = filter(mean.seed, !is.na(seed_mass))) +
    facet_wrap(~ grid) +
    theme_bw()
                  
  
############################################################################
  
  dat <- left_join(rat.pos, mean.seed)
  dat <- mutate(dat, grid = factor(grid))
  
## Define years
  dat$year <- 4            
  dat$year <- ifelse(dat$month <= 30, 3, dat$year)
  dat$year <- ifelse(dat$month <= 18, 2, dat$year)
  dat$year <- ifelse(dat$month <= 6, 1, dat$year)
  table(dat$month, dat$year)
  
# convert to decimal years
  dat$year.dec <- dat$month / 12
  

############################################################################
# model rat abundance
################################################################################
# matrices to store data

# N is mean estimated log abundance from mark-recapture model 
  N <- tapply(dat$xx, list(dat$year.dec, dat$grid), mean)

# Nse is standard error of log abundance from mark-recapture model
  Nse <- tapply(dat$xxse, list(dat$year.dec, dat$grid), mean)
# convert Nse to precision
  Ntau <- 1 / Nse^2

  seed <- tapply(dat$seed_mass, list(dat$year.dec, dat$grid), mean)
# centre seed
  nseed <- (seed - mean(seed, na.rm = T))

# measurement time
  year.dec <- as.numeric(rownames(seed))

# time gap between measurements
  td <- year.dec - lag(year.dec)
  
  nYears <- dim(Y)[1]
  nGrids <- dim(Y)[2]

################################################################################
# Gompertz model with same parameters for each grid
# here K = exp(b/d)

  lag.sp <- 0

  mod <- "model {

   # Observation model
   for(i in 1:nYears) {
     for(j in 1:nGrids) {
       Ntau[i, j] ~ dnorm(mu.Ntau, tau[5])       # model missing values for Ntau 
       N[i, j] ~ dnorm(X[i, j], Ntau[i, j])      # observation model for N

       # model seed as an overall mean (alpha1), an offset for each year (alpha2) and an offset for each grid (alpha3)
       seed[i, j] ~ dnorm(mu.r[i, j], tau[1])
       mu.r[i, j] <- alpha1 + alpha2[i] + alpha3[j]
     }
   }
   
  # priors for the observation and seed models
   for(i in 1:nYears) {
     alpha2[i] ~ dnorm(0, tau[2])
   }

   for(j in 1:nGrids) {
     alpha3[j] ~ dnorm(0, tau[3])
     X[1, j] ~ dnorm(0, 0.1) 
     simX[1, j] <- X[1, j]    
   }
   
   # Process model
   for(j in 1:nGrids) {
     for(i in 2:nYears) {

       X[i, j] ~ dnorm(predX[i, j], tau[4])
       predX[i, j] <- X[i-1, j] + (r + beta.seed * (seed[i-lag, j]) - d * X[i-1, j]) * td[i]

       # get independent estimates for the fitted process model
       simX[i, j] ~ dnorm(predsimX[i, j], tau[4])
       predsimX[i, j] <- simX[i-1, j] + (r + beta.seed * (seed[i-lag, j]) - d * simX[i-1, j]) * td[i]
     }
   }

  # priors
   mu.Ntau ~ dnorm(0, 0.001)T(0, )
   r ~ dnorm(0, 0.001)T(0, )
   beta.seed ~ dnorm(0, 0.01)
   
   d ~ dnorm(0, 0.01)T(0, )
   alpha1 ~ dnorm(0, 0.001)

   for(i in 1:5) {
     tau[i] <- 1 / (sigma[i] * sigma[i])
     sigma[i] ~ dunif(0, 100)
   }

  }"


# write model
  write(mod, "model.txt")

# run the model in JAGS with three chains
  mod.jags <- jags(model = "model.txt",
              data = list(N = N, Ntau = Ntau, seed = nseed, nGrids = nGrids, nYears = nYears, lag = lag.sp, td = td),
              param = c("r", "d", "beta.seed", "sigma", "X", "simX", "seed"),
              n.chains = 3,
              n.iter = 20000,
              n.burnin = 5000,
              parallel = TRUE)

  mod.sum <- mod.jags$summary
  mod.sum[1:30, c(1, 3, 7, 8, 9)]
  

# results of modelling missing seed values at each site
  estr <- data.frame(mod.sum[substr(rownames(mod.sum), 1, 4) == "seed", c(1, 4, 6, 8)])
  names(estr) <- c("mp", "lcl", "ucl", "rhat")
  estr$month <- rep(as.numeric(levels(factor(dat$month))), nGrids)
  estr$grid <- rep(levels(factor(dat$grid)), each = nYears)
  
# if rhat is non-missing then seed was estimated
  estr$est <- ifelse(is.na(estr$rhat), 0, 1)
  
  ggplot(estr, aes(y = mp, x = month)) +
    geom_line() +
    geom_point(aes(colour = factor(est))) +
    geom_errorbar(aes(ymin = lcl, ymax = ucl, colour = factor(est)), width = 0, alpha = 0.5) +
    facet_wrap(~ grid) +
    theme_classic() +
    theme(legend.position = "top")

# extract the model fitted to each grid and plot
  x <- data.frame(mod.sum[substr(rownames(mod.sum), 1, 1) == "X", c(1, 3, 7)])
  names(x) <- c("X", "lcl.X", "ucl.X")
  x$month <- rep(as.numeric(levels(factor(dat$month))), nGrids)
  x$grid <- rep(levels(factor(dat$grid)), each = nYears)

  simx <- data.frame(mod.sum[substr(rownames(mod.sum), 1, 4) == "simX", c(1, 3, 7)])
  names(simx) <- c("simX", "lcl.simX", "ucl.simX")
  simx$month <- rep(as.numeric(levels(factor(dat$month))), nGrids)
  simx$grid <- rep(levels(factor(dat$grid)), each = nYears)
  
# merge with the actual site data
  nsite <- left_join(dat, simx)
  nsite <- left_join(nsite, x)


# fitted process model is the red line
# circles are the grid data
# black line is the estimate of the true population having accounted for observation error (e.g. variation between grids at the same site)

  ggplot(nsite, aes(y = xx, x = month)) +
    geom_point(size = 3) +
    geom_line(aes(y = X), lwd = 1) +
    geom_line(aes(y = simX), colour = "red", lwd = 1) +
    facet_wrap(~ grid) +
    theme_classic() +
    theme(legend.position = "none")
    