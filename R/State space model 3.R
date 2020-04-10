
  library(tidyverse)
  library(jagsUI)
  
  setwd("c:\\users\\s429217\\onedrive\\data\\mpd\\state space\\data")
  
############################################################################
## Read in rat, mouse and possum abundance data on log scale

  rat.pos <- read.csv("c:\\users\\s429217\\onedrive\\data\\mpd\\state space\\data\\Rat possum mouse log N estimates.csv")
  rat.pos <- dplyr:::select(rat.pos, -X)
  glimpse(rat.pos)                                                                                  
  
# exclude first few trips
  table(rat.pos$month)
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
  seed.trap <- trap %>%
               group_by(seed.trap, grid, month) %>%
               summarise(seed_mass = sum(dry.wt, na.rm=T),
               seed_mass_poss_pref = sum(dry.wt[spp %in% c("Beitaw", "Eladen", "Rubcis", "Melram", "Hedarb", "Pencor")], na.rm=T),
               seed_mass_poss_avoid = sum(dry.wt[spp %in% c("Ripsca", "Prufer")], na.rm=T))

# convert to log scale
  seed.trap <- seed.trap %>%
               mutate(seed_mass = log(seed_mass + 0.001),
                      seed_mass_poss_pref = log(seed_mass_poss_pref + 0.001),
                      seed_mass_poss_avoid = log(seed_mass_poss_avoid + 0.001))
                     
# now summarise by grid and month calculating mean seed trap values
  mean.seed <- seed.trap %>%
               group_by(grid, month) %>%
               summarise(n = n(),
               seed_mass = mean(seed_mass, na.rm=T),
               seed_mass_poss_pref = mean(seed_mass_poss_pref, na.rm=T),
               seed_mass_poss_avoid = mean(seed_mass_poss_avoid, na.rm=T), 
               seed_mass_se = sd(seed_mass, na.rm = T)) 
               
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
# model abundance
################################################################################
# matrices to store data

# N is mean estimated log abundance from mark-recapture model 
  N.rat <- tapply(dat$rat_abun, list(dat$year.dec, dat$grid), mean)

# Nse is standard error of log abundance from mark-recapture model
  Nse.rat <- tapply(dat$rat_sd, list(dat$year.dec, dat$grid), mean)
# convert Nse to precision
  Ntau.rat <- 1 / Nse.rat^2

  N.poss <- tapply(dat$poss_abun, list(dat$year.dec, dat$grid), mean)
  Nse.poss <- tapply(dat$poss_sd, list(dat$year.dec, dat$grid), mean)
  Ntau.poss <- 1 / Nse.poss^2

  N.mouse <- tapply(dat$mouse_abun, list(dat$year.dec, dat$grid), mean)
  Nse.mouse <- tapply(dat$mouse_sd, list(dat$year.dec, dat$grid), mean)
  Ntau.mouse <- 1 / Nse.mouse^2

################################################################################
  seed <- tapply(dat$seed_mass, list(dat$year.dec, dat$grid), mean)
# centre seed
  nseed <- (seed - mean(seed, na.rm = T))

# measurement time
  year.dec <- as.numeric(rownames(seed))

# time gap between measurements
  td <- year.dec - lag(year.dec)
  
  nYears <- dim(N.rat)[1]
  nGrids <- dim(N.rat)[2]
  
# coding to identify Possum + rat removal grids
  prr <- ifelse(table(dat$treat, dat$grid)[2, ] > 0, 1, 0)

################################################################################
# Gompertz model with same parameters for each grid
# here K = exp(b/d)

  lag.sp <- 0

  mod <- "model {

   # Observation model
   for(i in 1:nYears) {
     for(j in 1:nGrids) {
       Ntau.rat[i, j] ~ dnorm(mu.Ntau.rat, tau[1])        # model missing values for Ntau   
       N.rat[i, j] ~ dnorm(X.rat[i, j], Ntau.rat[i, j])     

       Ntau.poss[i, j] ~ dnorm(mu.Ntau.poss, tau[2])           
       N.poss[i, j] ~ dnorm(X.poss[i, j], Ntau.poss[i, j])     

       Ntau.mouse[i, j] ~ dnorm(mu.Ntau.mouse, tau[2])           
       N.mouse[i, j] ~ dnorm(X.mouse[i, j], Ntau.mouse[i, j])     

       # model seed as an overall mean (alpha1), an offset for each year (alpha2) and an offset for each grid (alpha3)
       seed[i, j] ~ dnorm(mu.r[i, j], tau[3])
       mu.r[i, j] <- alpha1 + alpha2[i] + alpha3[j]
     }
   }
   
  # priors for the observation and seed models
   for(i in 1:nYears) {
     alpha2[i] ~ dnorm(0, tau[4])
   }

   for(j in 1:nGrids) {
     alpha3[j] ~ dnorm(0, tau[5])
     
     X.rat[1, j] ~ dnorm(0, 0.1) 
     simX.rat[1, j] <- X.rat[1, j]    
     
     X.poss[1, j] ~ dnorm(0, 0.1) 
     simX.poss[1, j] <- X.poss[1, j]    

     X.mouse[1, j] ~ dnorm(0, 0.1) 
     simX.mouse[1, j] <- X.mouse[1, j]    
   }
   
   # Process model
   for(j in 1:nGrids) {
     for(i in 2:nYears) {

     # rats
       X.rat[i, j] ~ dnorm(predX.rat[i, j], tau[6])
       predX.rat[i, j] <- X.rat[i-1, j] + (r.rat[j] + beta.seed.rat * (seed[i-lag, j]) - d.rat * X.rat[i-1, j] - m.rat * prr[j] + rat.poss * X.poss[i-1, j]) * td[i]

       simX.rat[i, j] ~ dnorm(predsimX.rat[i, j], tau[6])
       predsimX.rat[i, j] <- simX.rat[i-1, j] + (r.rat[j] + beta.seed.rat * (seed[i-lag, j]) - d.rat * simX.rat[i-1, j] - m.rat * prr[j] + rat.poss * simX.poss[i-1, j]) * td[i]

     # possums
       X.poss[i, j] ~ dnorm(predX.poss[i, j], tau[7])
       predX.poss[i, j] <- X.poss[i-1, j] + (r.poss[j] + beta.seed.poss * (seed[i-lag, j]) - d.poss * X.poss[i-1, j] + poss.rat * X.rat[i-1, j]) * td[i]

       simX.poss[i, j] ~ dnorm(predsimX.poss[i, j], tau[7])
       predsimX.poss[i, j] <- simX.poss[i-1, j] + (r.poss[j] + beta.seed.poss * (seed[i-lag, j]) - d.poss * simX.poss[i-1, j] + poss.rat * simX.rat[i-1, j]) * td[i]

     # mice
       X.mouse[i, j] ~ dnorm(predX.mouse[i, j], tau[8])
       predX.mouse[i, j] <- X.mouse[i-1, j] + (r.mouse[j] + beta.seed.mouse * (seed[i-lag, j]) - d.mouse * X.mouse[i-1, j] + mouse.rat * X.rat[i-1, j]) * td[i]

       simX.mouse[i, j] ~ dnorm(predsimX.mouse[i, j], tau[8])
       predsimX.mouse[i, j] <- simX.mouse[i-1, j] + (r.mouse[j] + beta.seed.mouse * (seed[i-lag, j]) - d.mouse * simX.mouse[i-1, j] + mouse.rat * simX.rat[i-1, j]) * td[i]
     }
   }

  # priors
  # allow rate of increase to vary among grids, modelled hierarchically
  for(i in 1:nGrids) {
    r.rat[i] ~ dnorm(mu.r.rat, tau[9])
    r.poss[i] ~ dnorm(mu.r.poss, tau[10])
    r.mouse[i] ~ dnorm(mu.r.mouse, tau[11])
  }
  
  # overall mean rate of increase across grids
   mu.r.rat ~ dnorm(0, 0.01)
   mu.r.poss ~ dnorm(0, 0.01)
   mu.r.mouse ~ dnorm(0, 0.01)
  
   mu.Ntau.rat ~ dnorm(0, 0.001)T(0, )
   beta.seed.rat ~ dnorm(0, 0.01)         # seed effect on rat rate of increase
   d.rat ~ dnorm(0, 0.01)           # density-dependent mortality of rats
   m.rat ~ dnorm(0, 0.01)           # mortality due to rat removal on grids where rats were killed
   rat.poss ~ dnorm(0, 0.01)              # effect of poss on rat

   mu.Ntau.poss ~ dnorm(0, 0.001)T(0, )
   beta.seed.poss ~ dnorm(0, 0.01)
   d.poss ~ dnorm(0, 0.01)
   poss.rat ~ dnorm(0, 0.01)

   mu.Ntau.mouse ~ dnorm(0, 0.001)T(0, )
   beta.seed.mouse ~ dnorm(0, 0.01)
   d.mouse ~ dnorm(0, 0.01)
   mouse.rat ~ dnorm(0, 0.01)

   alpha1 ~ dnorm(0, 0.001)

   for(i in 1:11) {
     tau[i] <- 1 / (sigma[i] * sigma[i])
     sigma[i] ~ dunif(0, 100)
   }

  }"


# write model
  write(mod, "model.txt")

# run the model in JAGS with three chains
# can modify the observation error associated with the mark-recapture data
  e <- 5

  mod.jags <- jags(model = "model.txt",
              data = list(N.rat = N.rat, Ntau.rat = Ntau.rat/e, 
                          N.poss = N.poss, Ntau.poss = Ntau.poss/e,
                          N.mouse = N.mouse, Ntau.mouse = Ntau.mouse/e,
                          seed = nseed, nGrids = nGrids, nYears = nYears, lag = lag.sp, td = td, prr = prr),
              param = c("mu.r.rat", "d.rat", "m.rat", "beta.seed.rat", "rat.poss", 
                        "mu.r.poss", "d.poss", "beta.seed.poss", "poss.rat", 
                        "mu.r.mouse", "d.mouse", "beta.seed.mouse", "mouse.rat",
                        "sigma", "seed",
                        "X.rat", "simX.rat", 
                        "X.poss", "simX.poss",
                        "X.mouse", "simX.mouse"),
              n.chains = 3,
              n.iter = 20000,
              n.burnin = 5000,
              parallel = TRUE)

  mod.sum <- mod.jags$summary
  mod.sum[1:30, c(1, 3, 7, 8, 9)]
  

################################################################################
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

################################################################################
# extract the model fitted to each grid and plot
  x <- data.frame(mod.sum[substr(rownames(mod.sum), 1, 5) == "X.rat", c(1, 3, 7)])
  names(x) <- c("X.rat", "lcl.X.rat", "ucl.X.rat")
  x$month <- rep(as.numeric(levels(factor(dat$month))), nGrids)
  x$grid <- rep(levels(factor(dat$grid)), each = nYears)

  simx <- data.frame(mod.sum[substr(rownames(mod.sum), 1, 8) == "simX.rat", c(1, 3, 7)])
  names(simx) <- c("simX.rat", "lcl.simX.rat", "ucl.simX.rat")
  simx$month <- rep(as.numeric(levels(factor(dat$month))), nGrids)
  simx$grid <- rep(levels(factor(dat$grid)), each = nYears)
  
# merge with the actual site data
  nsite <- left_join(dat, simx)
  nsite <- left_join(nsite, x)


# fitted process model is the red line
# circles are the grid data
# black line is the estimate of the true population having accounted for observation error (e.g. variation between grids at the same site)

  ggplot(nsite, aes(y = exp(rat_abun), x = month)) +
    geom_point(size = 3) +
    geom_line(aes(y = exp(X.rat)), lwd = 1) +
    geom_line(aes(y = exp(simX.rat)), colour = "red", lwd = 1) +
    facet_wrap(~ grid) +
    theme_classic() +
    theme(legend.position = "none")

################################################################################
# extract the model fitted to each grid and plot
  x <- data.frame(mod.sum[substr(rownames(mod.sum), 1, 6) == "X.poss", c(1, 3, 7)])
  names(x) <- c("X.poss", "lcl.X.poss", "ucl.X.poss")
  x$month <- rep(as.numeric(levels(factor(dat$month))), nGrids)
  x$grid <- rep(levels(factor(dat$grid)), each = nYears)

  simx <- data.frame(mod.sum[substr(rownames(mod.sum), 1, 9) == "simX.poss", c(1, 3, 7)])
  names(simx) <- c("simX.poss", "lcl.simX.poss", "ucl.simX.poss")
  simx$month <- rep(as.numeric(levels(factor(dat$month))), nGrids)
  simx$grid <- rep(levels(factor(dat$grid)), each = nYears)
  
# merge with the actual site data
  nsite <- left_join(dat, simx)
  nsite <- left_join(nsite, x)


# fitted process model is the red line
# circles are the grid data
# black line is the estimate of the true population having accounted for observation error (e.g. variation between grids at the same site)

  ggplot(nsite, aes(y = exp(poss_abun), x = month)) +
    geom_point(size = 3) +
    geom_line(aes(y = exp(X.poss)), lwd = 1) +
    geom_line(aes(y = exp(simX.poss)), colour = "red", lwd = 1) +
    facet_wrap(~ grid) +
    theme_classic() +
    theme(legend.position = "none")
    
################################################################################
# extract the model fitted to each grid and plot
  x <- data.frame(mod.sum[substr(rownames(mod.sum), 1, 7) == "X.mouse", c(1, 3, 7)])
  names(x) <- c("X.mouse", "lcl.X.mouse", "ucl.X.mouse")
  x$month <- rep(as.numeric(levels(factor(dat$month))), nGrids)
  x$grid <- rep(levels(factor(dat$grid)), each = nYears)

  simx <- data.frame(mod.sum[substr(rownames(mod.sum), 1, 10) == "simX.mouse", c(1, 3, 7)])
  names(simx) <- c("simX.mouse", "lcl.simX.mouse", "ucl.simX.mouse")
  simx$month <- rep(as.numeric(levels(factor(dat$month))), nGrids)
  simx$grid <- rep(levels(factor(dat$grid)), each = nYears)
  
# merge with the actual site data
  nsite <- left_join(dat, simx)
  nsite <- left_join(nsite, x)


# fitted process model is the red line
# circles are the grid data
# black line is the estimate of the true population having accounted for observation error (e.g. variation between grids at the same site)

  ggplot(nsite, aes(y = exp(mouse_abun), x = month)) +
    geom_point(size = 3) +
    geom_line(aes(y = exp(X.mouse)), lwd = 1) +
    geom_line(aes(y = exp(simX.mouse)), colour = "red", lwd = 1) +
    facet_wrap(~ grid) +
    theme_classic() +
    theme(legend.position = "none")        