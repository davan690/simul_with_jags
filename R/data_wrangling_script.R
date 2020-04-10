# data modification
library(tidyverse)
library(jagsUI)
library(kableExtra)
############################################################################
## Read in rat, mouse and possum abundance data on log scale
rat.pos <- read.csv("./data/raw/Rat possum mouse log N estimates.csv")

rat.pos <- dplyr:::select(rat.pos, -X)
# glimpse(rat.pos)

# choose species
rat.pos$xx <- rat.pos$rat_abun
rat.pos$xxse <- rat.pos$rat_sd

# choose treatments
# exclude "Possum + rat removal" treatment when analysing rats
# exclude the possum removal treatments when analysing possums

#remove
rat.pos <- filter(rat.pos, treat %in% c("Possum removal", "Control", "Stoat removal"))

# exclude first two trips
rat.pos <- filter(rat.pos, month > 2)

############################################################################
## Read in seed trap data
trap <- read.csv("./data/raw/seed trap data.csv")
trap <- dplyr:::select(trap, -X)
# glimpse(trap)

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

