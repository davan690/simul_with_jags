#jags model run function
#Anthony
#March2020

# Y = Y,
# Ytau = Ytau,
# seed = nseed,
# nGrids = nGrids,
# nYears = nYears,
# lag = lag.sp,
# td = td,
# killr = killr

creating_jags_output_rds <- function(model = model_jags, list_data) {
  
  mod.jags01 <- jags(
    model = model_jags,
    data = list_data,
    param = c("r", "d", "beta.seed", "sigma", "X", "simX", "seed", "rt", "g"),
    n.chains = 3,
    n.iter = 20000,
    n.burnin = 5000,
    parallel = TRUE
  )
  
  write_rds(mod.jags01, "./data/model_output_funtion1.RDS")
  
}

#example
##Source code
source("./R/data_wrangling_script.R")
source("./R/model_data_setup.R")

model_jags <- "./R/model_test01"

list_data <- list(
  Y = Y,
  Ytau = Ytau,
  seed = nseed,
  nGrids = nGrids,
  nYears = nYears,
  lag = lag.sp,
  td = td,
  killr = killr
)

creating_jags_output_rds(model = model_jags, list_data)




