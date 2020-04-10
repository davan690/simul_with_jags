############################################################################
# model rat abundance
# import data

################################################################################
# matrices to store data
Y <- tapply(dat$xx, list(dat$year.dec, dat$grid), mean)
Yse <- tapply(dat$xxse, list(dat$year.dec, dat$grid), mean)
seed <- tapply(dat$seed_mass, list(dat$year.dec, dat$grid), mean)
year.dec <- as.numeric(rownames(seed))
td <- year.dec - lag(year.dec)


Ytau <- 1 / Yse^2
td <- year.dec - lag(year.dec)
nseed <- (seed - mean(seed, na.rm = T))

nYears <- dim(Y)[1]
nGrids <- dim(Y)[2]

## select only grids with killr = 
levels(dat$treat)

dat.rat.treat <- filter(dat, treat == "Possum + rat removal")

table(dat.rat.treat$grid, dat.rat.treat$month)

dat.table <- mutate(dat, treat.p = ifelse(treat == "Possum + rat removal", 1, 0))


killr <- tapply(dat.table$treat.p, 
                list(dat.table$year.dec,
                     dat.table$grid), 
                mean)

# killp <- tapply(dat$seed_mass, list(dat$year.dec, dat$grid), mean)
# killm <- tapply(dat$seed_mass, list(dat$year.dec, dat$grid), mean)