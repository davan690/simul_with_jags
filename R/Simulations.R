# Author: Kevin See
# Purpose: Simulate data for PVAs and estimate parameters with MARSS
# Created: 9/21/10
# Last Modified: 4/9/14
# sets paths for computers with other software components
# this code block produces an enviroment warning
# but not sure how to sort depenancies still
# Feb2020
# myPaths <- .libPaths("C:/Program Files/R/R-3.6.2/library")
# myPaths <- c(myPaths)
# .libPaths(myPaths)  # add new path
# .libPaths()
library(MARSS)

#---------------------------------------------
# Fix known parameters
#---------------------------------------------
# set restrictions / constraints
sigPlusTau = 0.1
sig2Tau = c(0.05, 0.2, 1, 5)

# calculate values
sig = (sig2Tau * sigPlusTau) / (1 + sig2Tau)
tau = sigPlusTau - sig

# list parameters
Qs = sig	# Qs are variances
Rs = tau	# Rs are variances
Us = c(-0.02, -0.04)
tSteps = c(5, 15, 30, 45, 60)	# number of years for parameterization
maxTsteps = max(tSteps)
totobs = c(15,30,60,75,90)	# number of total observations

nobs = matrix(NA, length(totobs), length(tSteps), dimnames=list(totobs, tSteps))
for(i in 1:length(totobs)) {
	for(j in 1:length(tSteps)) {
		nobs[i,j] = totobs[i] / tSteps[j]
	}
}

maxN = max(nobs)
burn = 100 		# length of burn in period
B = matrix(1, 1)
x0 = matrix(10, 1, dimnames=list('x0'))



#---------------------------------------------
# Run simulations and estimate parameters
#---------------------------------------------

set.seed(3)
nsim = 500
nSim = nsim * length(totobs) * length(tSteps)


for(i in 1:length(Us)) {
	for(j in 1:length(Qs)) {
	# simulate some data
		U = matrix(Us[i], 1, dimnames=list('U'))
		Q = matrix(Qs[j], 1, dimnames=list('Q'))
		R = matrix(Rs[j], 1, 1, dimnames=list('diag'))
		# set up array to hold estimated paramters
		ParEst = array(NA, dim=c(nSim, 8), dimnames=list(1:nSim, c('U', 'Q2', 'R2', 'TotObs', 'TSteps', 'EstU', 'EstQ', 'EstR')))

		ParEst[1:nSim,1:3] = matrix(c(U,Q,Rs[j]), nSim, 3, byrow=T)	
		for(t in 1:length(tSteps)) {
			time.steps = tSteps[t]
			for(l in 1:length(totobs)) {
				N = totobs[l] / time.steps
				miss.locs = array(1, dim=c(ceiling(N), burn+time.steps, nsim))
				for(sim in 1:nsim) {
					locs = sort(sample(1:(ceiling(N) * time.steps), totobs[l]))
					test = array(1, dim=c(1, ceiling(N)*time.steps))
					test[-locs] = NA
					miss.locs[,(burn+1):(burn+time.steps),sim] = matrix(test, ceiling(N), byrow=time.steps)
				}
				# simulate data
				constraints = list(Z=rep(as.factor(1), ceiling(N)), R='diagonal and equal', Q='diagonal and equal', A='zero', B='identity')
				# this sets up the correct structure of the marssMLE object
				my.marssMLE = MARSS(matrix(rnorm(totobs[l], x0, sqrt(Q+R[1,1])), ceiling(N)), model=constraints, fit=F)						# set the parameters
				my.marssMLE$par$Q = Q
				my.marssMLE$par$R = R
				my.marssMLE$par$U = U
				my.marssMLE$par$x0 = x0
				my.marssMLE$par$Z = my.marssMLE$par$A = my.marssMLE$par$B = my.marssMLE$par$V0 = matrix(NA, 0, 1)
				# simulate states and observations
				sims = MARSSsimulate(my.marssMLE, nsim=nsim, tSteps=burn+time.steps, miss.loc=miss.locs)
				# extract the simulated observations
				obs = array(sims$sim.data[,(burn+1):(burn+time.steps),], dim=c(ceiling(N), time.steps, nsim))	
				# save simulated observations
				write.table(obs, file=paste('SimulatedData_U', Us[i], '_', j, '_TS_', time.steps, '_TotObs_', totobs[l], '.txt', sep=''))
				# estimate parameters
				for(h in 1:nsim) {
					# what row are we recording in the parameter estimation array?
					startRow = nsim * length(totobs) * (which(tSteps==tSteps[t]) - 1) + (nsim * (which(totobs==totobs[l]) - 1))
					ParEst[startRow + h,4:5] = c(totobs[l], time.steps)

					# show progress
					print(c(Us[i], j, time.steps, totobs[l], h))

					# constraints = list(Z=factor(c(rep(1,ceiling(N)))), R='diagonal and equal')
					# fixedlist = list(A=array(0, dim=c(ceiling(N),1)))
					my.MLEobj=MARSS(obs[1:ceiling(N),1:time.steps,h], model=constraints, silent=T, method='kem')
					ParEst[startRow + h,6:8] = c(my.MLEobj$par$U, my.MLEobj$par$Q, my.MLEobj$par$R[1,1])	
				}
			}
		}
	# save results
	write.table(ParEst, paste('Estimates//U',Us[i], '_', j,'.txt', sep=''))               
	}
}


#---------------------------------------------
# Estimate probability of quasi-extinction
#---------------------------------------------

# Function to estimate the probability of quasi-extinction
# some of the parameters must be set in your code beforehand (e.g. Th,) and "x" contains mu and sig2Q
ProbExt3 = function(x, perDec){
	a= -log(1-perDec)
	u = as.numeric(-x[1] * sqrt(Th) / sqrt(x[2]))
	v = as.numeric(a / (sqrt(x[2]) * sqrt(Th)))
	pe = pnorm(u-v) + exp(2*u*v)*pnorm(-(u+v))
	Pe = ifelse(exp(2*u*v) == Inf, 0, pe)
	return(Pe)
}

#---------------------------------------------
# set percent declines
PerDec = c(0.3, 0.5, 0.8)
# set time horizons
TimeHz = c(5,10,20,30,50,100)

# calculate risk of decline for 3 levels of decline
PE=array(NA, dim=c(dim(ParEst)[1], length(PerDec), length(TimeHz)), dimnames=list(NULL, PerDec, TimeHz))
for(i in 1:length(PerDec)){
	for(j in 1:length(TimeHz)){
		Th = TimeHz[j]
		PE[,i,j] = apply(ParEst[,6:7],1,ProbExt3, perDec=PerDec[i])
		}
	}
write.table(PE, '../data/PE_Estimates.txt')

