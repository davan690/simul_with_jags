
## Grassland Earless Dragon mark-recapture data from Cookanalla, Jerra East, and Majura Training Area

#load necessary libraries

  library(tidyverse)
  library(lubridate)
  library(jagsUI)

  unlogit <- function(x) exp(x) / (1 + exp(x))

# setwd("C:/Users/s438134/OneDrive/Casual Contract UC/Organised Directory/data")
#  setwd("c:\\users\\s429217\\onedrive\\pgrad\\emily stringer\\geds")
# setwd("C:/Users/s438134/Downloads/OneDrive-2019-11-08")
## read dataframe and rename
  ged <- read.csv("./data/ged_aug2019_Airport.csv", stringsAsFactors = F) %>%
         mutate(svl = ifelse(svl %in% c("Not handled", "", " "), NA, svl),
                svl = ifelse(is.na(svl), NA, as.numeric(svl)))

# remove blanks in site name
  ged$site <- str_remove_all(ged$site, " ")

  glimpse(ged)

  names(ged)[names(ged) == 'grid'] <- 'grid.id'

  ged <- mutate(ged,
                grid = paste(site, grid.id),
                date = ymd(date),
                month = as.character(month(date, label = T)))

  glimpse(ged)

  table(ged$age)
#  ged <- filter(ged, age %in% c("SA", "A"))

################################################################################
# calculate the number of surveys on each grid in each trip

  surveys <- read.csv("./data/ged_survey_numbers2019_2.csv", stringsAsFactors = F)
  glimpse(surveys)

# 2019 surveys
  surveys19 <- data.frame(site = rep(c("Majura", "JerraEast", "JerraWest", "Cookanalla"), each = 4),
                          grid.id = c("A", "B", "C", "D", "O", "P", "U", "Y", "E", "F", "G", "H", "I", "J", "K", "L"),
                          year = 2019,
                          checks = rep(c(18, 18, 18, 11), each = 4)) %>%
               mutate(site = as.character(site),
                      grid.id = as.character(grid.id))
  glimpse(surveys19)

  surveys <- bind_rows(surveys, surveys19) %>%
             mutate(grid = paste(site, grid.id),
                    gridyear = paste(grid, year)) %>%
             select(gridyear, checks) %>%
             arrange(gridyear)

  head(surveys)
################################################################################
## for each individual, calculate the number of times it was captured on each grid in each year
  ind <- ged %>%
         mutate(gridyear = paste(grid, year)) %>%
         group_by(gridyear, age, em.animal) %>%
         summarise(n = n(),
                   svl = min(svl, na.rm = T)) %>%
         mutate(recap = ifelse(n > 1, 1, 0),
                svl = ifelse(svl == "Inf", NA, svl))

# merge with number of checks for each survey
  ind <- full_join(ind, surveys)

# drop Bonshaw Z in 2011 because this is not in the survey list
  ind <- filter(ind, gridyear != "Bonshaw Z 2011")

  table(ind$n)
  table(ind$recap)
  table(ind$age, ind$recap)

  glimpse(ind)

#############################################################################
# the minimum number of animals caught on each grid on each trip
# and number of recaptures
  min.ged <- group_by(ind, gridyear) %>%
             summarise(recap = sum(recap),
                       n = n())

################################################################################
# use the data frame ind to create an array with the number of captures for each individual on each grid
# maximum number caught on a grid
  max.n <- max(min.ged$n)
  max.n

# set up array Y with rows i = grid x year, columns = individuals
# set up Y columns to be 5 times the length of max.n to allow sufficient number of augmented individuals
  Y <- matrix(0, nrow = nrow(min.ged), ncol = max.n * 5 + 1)

# SVL for each animal
  SVL <- matrix(NA, nrow = nrow(min.ged), ncol = max.n * 5 + 1)

  # fill in Y with the individuals captured on each grid in each year
  # first get a list of grids x years
  gy <- levels(factor(ind$gridyear))


# check these are all equal
  length(gy)
  nrow(min.ged)
  nrow(surveys)


# loop through these to get values
  for(i in 1:length(gy)) {
    sub.ind <- filter(ind, gridyear == gy[i])
    cap <- sub.ind$n
    ncap <- length(cap)
    Y[i, 1:ncap] <- cap
    svl <- sub.ind$svl
    SVL[i, 1:ncap] <- svl
  }

# set the missing values to 0
  Y <- ifelse(is.na(Y), 0, Y)

# centre and scale SVL
  SVL <- (SVL - mean(SVL, na.rm = T)) / sd(SVL, na.rm = T)

  Y[1:20, 1:10]
  SVL[1:20, 1:10]

# some checks
  table(Y)
  table(ind$n)

# vector with the number of checks for each year
  nc <- ind %>%
        select(gridyear, checks) %>%
        unique() %>%
        arrange(gridyear)

  J <- nc$checks
  J
  length(J)

################################################################################
# row i is grid x year
# column j is individual

# latent indicator variable
  z <- ifelse(Y > 0, 1, NA)

################################################################################
# set a prior on the number of individuals that could be present linked to the number
# caught by allowing the number of augmented individuals to vary
# = 4 times the number caught or a minimum of 5
# to do this we create a variable N_animals that is the number of augmented individuals
# for each grid x year to sample up to

  N_animals <- min.ged$n * 5
  N_animals <- ifelse(N_animals < 5, 5, N_animals)

# replace z values that are >N_animals with 0 as these are assumed absent animals
  for(i in 1:nrow(z)) {
    z[i, (N_animals[i] + 1):ncol(z)] <- 0
  }

  N_gridyear <- nrow(Y)


# loop to get values
  for(i in 1:length(gy)) {
    sub.ind <- filter(ind, gridyear == gy[i])
    cap <- sub.ind$n
    ncap <- length(cap)
    if(ncap > 0) Y[i, 1:ncap] <- cap
  }
################################################################################
# JAGS model

  CR_ged <- "model {

    for(i in 1:N_gridyear){
      for(j in 1:N_animals[i]){

      # model missing values of SVL
        SVL[i, j] ~ dnorm(mu.svl, tau.svl)

        Y[i, j] ~ dbin(mu[i, j], J[i])
        mu[i, j] <- z[i, j] * p[i, j]
        z[i, j] ~ dbern(psi[i])
        logit(p[i, j]) <- lp[i, j]

        # allow for heterogeneity among individuals in capture probability
        lp[i, j] ~ dnorm(mu.ind[i, j], tau.ind)

        # model heterogeneity as a function of SVL
        mu.ind[i, j] <- mu.p + b.svl * SVL[i, j] + mu.gridyear[i]
      } #j

      # allow for differences among grids x years in capture probability
      mu.gridyear[i] ~ dnorm(0, tau.gridyear)

      # prior for presence probability differs betwen trips and grids
      psi[i] ~ dunif(0, 1)

      # Derived parameters: the estimated number of individuals on each grid at each trip
      N[i] <- sum(z[i, ])

    } #i

    # priors

    mu.p ~ dnorm(0, 0.01)
    mu.svl ~ dnorm(0, 0.01)
    b.svl ~ dnorm(0, 0.01)

    tau.ind <- 1 / (sigma.ind * sigma.ind)
    sigma.ind ~ dunif(0, 10)
    tau.svl <- 1 / (sigma.svl * sigma.svl)
    sigma.svl ~ dunif(0, 100)
    tau.gridyear <- 1 / (sigma.gridyear * sigma.gridyear)
    sigma.gridyear ~ dunif(0, 10)

      }"  #model finish

################################################################################
#RUN MODEL------------------------------------------------

  write(CR_ged, "CR_ged.txt")

  mod <- jags(model.file = "CR_ged.txt",
              data = list(Y = Y, J = J, z = z, N_gridyear = N_gridyear, N_animals = N_animals,
                          SVL = SVL),
              parameters.to.save = c("N", "mu.p", "b.svl", "sigma.ind", "sigma.gridyear", "mu.gridyear"),
              n.chains = 3,
              n.iter = 5000,
              n.burnin = 1000,
              n.thin = 1,
              parallel = T)

  mod
  mod.sum <- mod$summary

################################################################################
# overall capture probability
  cp <- mod.sum[substr(rownames(mod.sum), 1, 4) == "mu.p", 1]
  unlogit(cp)

# capture probability for each grid x year
  gyp <- mod.sum[substr(rownames(mod.sum), 1, 4) == "mu.g", 1]
  gyp <- unlogit(cp + gyp)

  hist(gyp)

  outp <- data.frame(p = gyp)
  outp$gridyear <- gy
  a <- strsplit(gy, split = " ")
  outp$site <- unlist(lapply(a, function(x) x[1]))
  outp$grid.id <- unlist(lapply(a, function(x) x[2]))
  outp$year <- as.numeric(unlist(lapply(a, function(x) x[3])))

  head(outp)
  plot(outp$p ~ factor(outp$site))

################################################################################
# extract the N values
  out <-  data.frame(mod.sum[substr(rownames(mod.sum), 1, 1) == "N", c(1, 2, 3, 7)])
  names(out) <- c("N", "se", "lcl", "ucl")
  out$gridyear <- gy
  a <- strsplit(gy, split = " ")
  out$site <- unlist(lapply(a, function(x) x[1]))
  out$grid.id <- unlist(lapply(a, function(x) x[2]))
  out$year <- as.numeric(unlist(lapply(a, function(x) x[3])))
 write.csv(out, file = "C:/phd-repository/ged-project/data/RD_final_Output.csv")
# renumber grids
  gn <- out %>%
        select(site, grid.id) %>%
        unique() %>%
        group_by(site) %>%
        mutate(grid.n = row_number())

  out <- full_join(out, gn)

  out$n.cap <- min.ged$n

  head(out)

  ggplot(out, aes(y = N, x = year, colour = site)) +
    geom_line(aes(y = n.cap), col = "red") +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = lcl, ymax = ucl), colour = "grey", width = 0) +
    facet_grid(site ~ grid.n) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none")


# for eacg site
  ss <- levels(factor(out$site))

  #pdf(file = "c:\\users\\s429217\\onedrive\\pgrad\\emily stringer\\geds\\Abudnance estimates by site.pdf")
  pdf(file = "C:\\Users\\s438134\\OneDrive\\Pictures\\Abudnance estimates by site.pdf")
  for(i in 1:length(ss)) {

  p1 <- ggplot(filter(out, site == ss[i]), aes(y = N, x = year)) +
    geom_point(aes(y = n.cap), col = "red") +
    geom_line(aes(y = n.cap), col = "red") +
    geom_errorbar(aes(ymin = lcl, ymax = ucl), colour = "grey", width = 0) +
    geom_point() +
    geom_line() +
    facet_wrap(~ grid.id) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none") +
    ggtitle(ss[i])

    print(p1)
  }
  dev.off()
###############################################################################################################
####################################                            ###############################################
#--------------------------------------- state-space model -----------------------------------------------------
  out2 <- read.csv("C:/phd-repository/ged-project/data/GED_abundance_model_output.csv")
  head(out2)
  summary(out2$site)
  dat.rain <- read.csv("CanberraAairport_monthlyweather.csv")

  # use rainfall in particular months of each year
  use.rain <- dat.rain %>%
    group_by(year) %>%
    summarise(rain = sum(rain[month %in% 8:12]))

  # abundance per gird
  # remove low population grids first
  n <- out2 %>%
    ungroup() %>%
    filter((site %in% c("Cookanalla", "JerraEast", "Majura"))) %>%
    filter(!(site %in% c("Cookanalla", "JerraEast") & grid.no. %in% c(3, 4))) %>%
    mutate(sitegrid = paste(site, grid.id)) %>%
    select(year, sitegrid, N) %>%
    spread(sitegrid, N)


  # get rainfall values for each year
  rain <- use.rain$rain[match(n$year, use.rain$year)]
  year <- n$year

  Y <- as.matrix(n[, -1])
  lY <- log(Y)

  nGrids <- dim(Y)[2]
  nYears <- dim(Y)[1]
  myplotdata <- out2[paste(out2$site, out2$grid.id) %in% colnames(n),]
  myplotdata2 <- merge(myplotdata, use.rain2, by = "year", all.x = TRUE)

  ggplot(myplotdata2, aes(x = year, y = N, fill = site)) +
    geom_bar(stat="identity") +
    facet_wrap(~site)
    theme_minimal()
  ggplot(use.rain, aes(x = year, y = rain)) +
    geom_point()+
    geom_line()
  use.rain2$rainlag <- 0
   use.rain2 <- data.frame(use.rain)
  for (i in 1:(nrow(use.rain2)-1)) {
    lagrain  <- use.rain2[i,2]
    use.rain2$rainlag[i+1] <- lagrain
  }
  use.rain2
   myplotdata2

  mdata <-  myplotdata2 %>%
   group_by(year) %>%
     summarise(mean_N = mean(N), mean_rain = mean(rainlag))

   plot(mean_N ~ mean_rain, data = mdata)
   abline(lmodel)
   lmodel <-glm(mean_N ~ mean_rain, data = data.frame(mdata)[-14])
   plot(lmodel, which = 1)

  ################################################################################
  # Gompertz model
  # here K = exp(b/d)

  # specify the rainfall lag in the data statement

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
       predX[i] <- X[i-1] + b + beta.rain * rain[i-lag] - d * X[i-1]

       simX[i] ~ dnorm(predsimX[i], tau[2])    # fitted model
       predsimX[i] <- simX[i-1] + b + beta.rain * rain[i-lag] - d * simX[i-1]
   }

   b ~ dnorm(0, 0.001)T(0, )
   beta.rain ~ dnorm(0, 0.01)
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
                   data = list(Y = lY, nGrids = nGrids, nYears = nYears, rain = rain, lag = 1),
                   param = c("b", "d", "beta.rain", "sigma", "X", "simX"),
                   n.chains = 3,
                   n.iter = 10000,
                   n.burnin = 5000,
                   parallel = TRUE)

  mod.sum2 <- mod.jags$summary
  mod.sum2
  mod.jags


  install.packages("dotwhisker")
  library(dotwhisker)

  dwplot(tidy(mod.sum2))
  dwplot(mod.sum2, vline = geom_vline(xintercept = 0, colour = "grey50", linetype = 2)) +
    xlab("Coefficient")
  summary(model)
  para <- mod.jags$parameters
  whiskerplot(mod.jags, parameters = para[1:length(para)-1], quantiles=c(0.025,0.975), zeroline=TRUE)
  # mean over all grids
  av <- data.frame(n = rowMeans(Y, na.rm = T), year = n$year)


  x <- data.frame(mod.sum[substr(rownames(mod.sum), 1, 1) == "X", c(1, 3, 7)])
  x <- exp(x)
  names(x) <- c("mp", "lcl", "ucl")
  x$year <- n$year

  simx <- data.frame(mod.sum[substr(rownames(mod.sum), 1, 4) == "simX", c(1, 3, 7)])
  simx <- exp(simx)
  names(simx) <- c("mp", "lcl", "ucl")
  simx$year <- n$year

  ggplot(av, aes(y = n, x = year)) +
    geom_point() +
    geom_line(data = x, aes(y = mp)) +
    geom_line(data = simx, aes(y = mp), colour = "red", lwd = 1.1) +
    theme_bw() +
    theme(legend.position = "none")

  plot(log(simx$mp) ~ av$n)
  model <- lm(log(simx$mp) ~ av$n)
  plot(model, which = c(1,2,4))
  summary(model)

