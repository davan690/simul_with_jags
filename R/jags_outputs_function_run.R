##plotting results from model fit
#April2020
#Anthony davidson
# mod_out <- read_rds("./data/model_output_funtion1.RDS")
# mod.sum <- mod_out$summary
# mod.sum
# mod_out$DIC
# mod.sum[1:30, c(1, 3, 7, 8, 9)]

jags_outputs_plots_func <- function(mod.sum = mod.sum){
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

}

mod_out <- read_rds("./data/mod_jags_april3-TO-1v3.rds")
mod.sum <- mod_out$summary
mod.sum
mod_out$DIC
mod.sum[1:30, c(1, 3, 7, 8, 9)]

jags_outputs_plots_func(mod.sum = mod.sum)
