# Author: Kevin See
# Purpose: Compare estimates of mu, sigma and tau using MARSS, REML 1 and REML 2
# Created: 10/23/10
# Last Modified: 4/9/14
# Notes: REML 1 based on 1st differences, REML 2 based on 2nd differences

library(MASS)
library(MARSS)

#--------------------------------------------------------------
# likelihood functions for REML 1 and REML 2
#--------------------------------------------------------------
#  Log-likelihood for REML estimation from 1st differences. 
negloglike.reml1=function(theta,wt,pt)  
{ 
   muu=theta[1];
   sigmasq=exp(theta[2]);         #  Constrains ssq > 0.  
   tausq=exp(theta[3]);           #  Constrains tsq > 0. 
   q=length(wt)/pt; 
   qp1=q+1; 
   vx=matrix(0,qp1, qp1)
   for(ti in 1:q){
   	vx[(ti+1):qp1, (ti+1):qp1]=matrix(1,1, (qp1-ti))*tt[ti+1]
   	}
   Sigma.mat=sigmasq*vx			# this is the covariance matrix of X
   J=matrix(1/pt,pt,pt);           #  pXp matrix of 1/pt's. 
   I=diag(qp1);            #  (q+1)X(q+1) identity matrix. 
   D=kronecker(I,J);      #  kronecker product (block diagonal). 
   I.2=cbind(matrix(0,q*pt,pt),diag(q*pt)); 
   I.2=rbind(I.2,matrix(0,pt,pt*qp1)); 
   D=-D+I.2;     
   D=D[1:(q*pt),];   #  This is the D transformation matrix for REML. 
   j=matrix(1,pt,1);           #  pX1 column vector of ones. 
   j.2=matrix(1,pt*qp1,1);   #  p*(q+1) X 1 column vector of ones. 
   C=kronecker(I, j);     #  C matrix in the multivariate normal distribution of YP.reml . 
   V=C%*%Sigma.mat%*%t(C)+tausq*diag(qp1*pt);   #  Var-cov matrix of YP.reml . 
   Phi.mat=D%*%V%*%t(D);     #  Var-cov matrix of WP.t. 
   Phiinv.mat=ginv(Phi.mat); 
   llikew=-(p*q/2)*log(2*pi)-(1/2)*log(det(Phi.mat))- 
      (1/2)*t(wt)%*%Phiinv.mat%*%wt;   #  REML log-likelihood. 
   ofn=-llikew; 
   return(ofn); 
} 

#  Log-likelihood for REML estimation from 2nd differences. 
negloglike.reml2=function(theta,ut,pt)  
{ 
#   muu=theta[1];
   sigmasq=exp(theta[1]);         #  Constrains ssq > 0.  
   tausq=exp(theta[2]);           #  Constrains tsq > 0. 
   qm1=(length(ut)/pt);
   q=qm1+1 
   qp1=q+1;
   vx=matrix(0,qp1, qp1)
   for(ti in 1:q){
   	vx[(ti+1):qp1, (ti+1):qp1]=matrix(1,1, (qp1-ti))*tt[ti+1]
   	}
   Sigma.mat=sigmasq*vx			# this is the covariance matrix of X
   J=matrix(1/pt,pt,pt);           #  pXp matrix of 1/pt's. 
   I=diag(qp1);            #  (q+1)X(q+1) identity matrix. 
   D=kronecker(I,J);      #  kronecker product (block diagonal). 
   I.2=cbind(matrix(0,q*pt,pt),diag(q*pt)); 
   I.2=rbind(I.2,matrix(0,pt,pt*qp1)); 
   D=-D+I.2;     
   D=D[1:(q*pt),];   #  This is the D transformation matrix for 1st differences REML. 
   D2=D[1:(pt*qm1), 1:(pt*q)]; # this is the D2 transformation matrix for 2nd differences REML
   j=matrix(1,pt,1);           #  pX1 column vector of ones. 
   j.2=matrix(1,pt*qp1,1);   #  p*(q+1) X 1 column vector of ones. 
   C=kronecker(I, j);     #  C matrix in the multivariate normal distribution of YP.reml . 
   V=C%*%Sigma.mat%*%t(C)+tausq*diag(qp1*pt);   #  Var-cov matrix of YP.reml . 
   Phi.mat=D2%*%D%*%V%*%t(D)%*%t(D2);     #  Var-cov matrix of UP.t. 
   Phiinv.mat=ginv(Phi.mat); 
   llikew=-(p*qm1/2)*log(2*pi)-(1/2)*log(det(Phi.mat))- 
      (1/2)*t(ut)%*%Phiinv.mat%*%ut;   #  REML log-likelihood. 
   ofn=-llikew; 
   return(ofn); 
} 


#--------------------------------------------------------------
# Fix known parameters
#--------------------------------------------------------------
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

#--------------------------------------------------------------
# now run many simulations to compare REML and MARSS estimates 
#--------------------------------------------------------------
# how many simulations to run
nsim = 1000

# fix some of the parameters used for simulations
U = matrix(Us[2], 1, dimnames=list('U'))
time.steps = tSteps[3]

# loop over different process:non-process variance ratios and obs / time step
set.seed(4)
for(j in 1:4) {
	for(k in totobs[c(2,3,5)]) {
		N = totobs[k] / time.steps
		# set some initial parameter values for REML optimization
		mu0=Us[2];            #  Value of "a" to be used. 
		ssq0=Qs[j];           #  Value of "ssq" to be used. 
		tsq0=Rs[j];           #  Value of "tsq" to be used. 
		
		# simulate the data
		Q = matrix(Qs[j], 1, dimnames=list('Q'))
		R = matrix(Rs[j], 1, 1, dimnames=list('diag'))
		# set model constraints
		constraints = list(Z=rep(as.factor(1), ceiling(N)), R='diagonal and equal', Q='diagonal and equal', A='zero', B='identity')
		# this sets up the correct structure of the marssMLE object
		my.marssMLE = MARSS(matrix(rnorm(totobs[2], x0, sqrt(Q+R[1,1])), ceiling(N)), model=constraints, fit=F)						
		# set the parameters
		my.marssMLE$par$Q = Q
		my.marssMLE$par$R = R
		my.marssMLE$par$U = U
		my.marssMLE$par$x0 = x0
		my.marssMLE$par$Z = my.marssMLE$par$A = my.marssMLE$par$B = my.marssMLE$par$V0 = matrix(NA, 0, 1)
		# simulate the states and observations
		sims = MARSSsimulate(my.marssMLE, nsim=nsim, tSteps=burn+time.steps)
		obs = array(sims$sim.data[,(burn+1):(burn+time.steps),], dim=c(ceiling(N), time.steps, nsim))
		
		# set up array to capture parameter estimates
		ParEst = array(NA, dim=c(nsim, 12), dimnames=list(NULL, c('U.reml1', 'U.reml2', 'U.em', 'Q.reml1', 'Q.reml2', 'Q.em', 'R.reml1', 'R.reml2', 'R.em', 'X0.reml1', 'X0.reml2', 'X0.em')))
		
		# set which simulation to test
		for(h in 1:nsim) {
			# set up the data correctly
			qplus1 = dim(obs)[2]
			q=qplus1-1
			YP.t = obs[,,h]
			p = dim(obs)[1]
			Time.t=c(1:dim(obs)[2])
			tt=Time.t-Time.t[1]
			
			#  Matrices for calculating multivariate normal  likelihood for REML estimates. 
			YP.reml=matrix(YP.t,p*qplus1,1);   #  Replicated samples "stacked" in a vector. 
			J.p=matrix(1/p,p,p);           #  pXp matrix of 1/p's. 
			I.p=diag(qplus1);            #  (q+1)X(q+1) identity matrix. 
			D.reml=kronecker(I.p,J.p);   #  kronecker product (block diagonal). 
			I.temp=cbind(matrix(0,q*p,p),diag(q*p)); 
			I.temp=rbind(I.temp,matrix(0,p,p*qplus1)); 
			D.reml=-D.reml+I.temp; 
			D.reml=D.reml[1:(q*p),];     #  This is the D transformation matrix 
			                            #   for REML. 
			D2.reml = D.reml[1:(p*(q-1)), 1:(p*q)]; # this is the 2nd difference transformation matrix
			WP.t=D.reml%*%YP.reml ;       #  The REML-transformed observations. 
			UP.t=D2.reml%*%WP.t                            
			j.p=matrix(1,p,1);           #  pX1 column vector of 1's. 
			#j.pXqp1=matrix(1,p*qplus1,1);   #  p*(q+1) X 1 column vector of ones. 
			C.ml=kronecker(I.p, j.p);     #  C matrix in the multivariate normal 
			                            #    distribution of YP.reml . 
			
			
			# Calculate REML estimates                            
			GSSRSreml=optim(par=c(mu0, log(ssq0), log(tsq0)), negloglike.reml1 ,NULL,method="Nelder-Mead",wt=WP.t,pt=p)
			GSSRSreml=c(GSSRSreml$par[1],exp(GSSRSreml$par[2]), exp(GSSRSreml$par[3]),-GSSRSreml$value)
			mu.reml=GSSRSreml[1];             #  These are the REML estimates. 
			ssq.reml=GSSRSreml[2];           #          -- 
			tsq.reml=GSSRSreml[3];           #          -- 
			lnlike.reml=GSSRSreml[4];        #  This is the maximized log-likelihood 
			
			#  Calculate REML estimate of the parameter "beta0". 
			vx=matrix(0,qplus1, qplus1)
			for(ti in 1:q){
			vx[(ti+1):qplus1, (ti+1):qplus1]=matrix(1,1, (qplus1-ti))*tt[ti+1]
			}
			Sigma.mat=ssq.reml*vx	
			V.mat=C.ml%*%Sigma.mat%*%t(C.ml)+tsq.reml*diag(qplus1*p); 
			Vinv.mat=ginv(V.mat);
			j1 = matrix(1,p*qplus1,1) 
			beta0.reml = (t(j1)%*%Vinv.mat%*%(YP.reml-rep(tt*mu.reml,each=p))) / (t(j1)%*%Vinv.mat%*%j1)
			
			print(c(h, 'REML 1 done'))
			
			# Calculate REML estimates with 2nd differences                          
			GSSRSreml2=optim(par=c(log(ssq0), log(tsq0)), negloglike.reml2 ,NULL,method="Nelder-Mead",ut=UP.t,pt=p)
			GSSRSreml2=c(exp(GSSRSreml2$par[1]), exp(GSSRSreml2$par[2]),-GSSRSreml2$value)
			ssq.reml2=GSSRSreml2[1];           # This is the 2nd diff REML estimate 
			tsq.reml2=GSSRSreml2[2];			 # This is the 2nd diff REML estimate 
			lnlike.reml2=GSSRSreml2[3];        #  This is the maximized log-likelihood
			
			# Calculate REML 2nd diff estimate of mu
			Sigma.mat=ssq.reml2*vx	
			V.mat=C.ml%*%Sigma.mat%*%t(C.ml)+tsq.reml2*diag(qplus1*p); 
			Vinv.mat=ginv(V.mat);
			j2 = matrix(1, p*q,1)
			mu.reml2 = t(j2)%*%D.reml%*%Vinv.mat%*%t(D.reml)%*%WP.t / (t(j2)%*%D.reml%*%Vinv.mat%*%t(D.reml)%*%j2)
			
			# Calculate REML 2nd diff estimate of 'beta0'
			beta0.reml2 = (t(j1)%*%Vinv.mat%*%(YP.reml-rep(tt*mu.reml2,each=p))) / (t(j1)%*%Vinv.mat%*%j1)
			
			print(c(h, 'REML 2 done'))
			
			# get ML estimates with MARSS
			my.MLEobj=MARSS(YP.t, model=constraints, silent=T)
			MLpar = my.MLEobj$par
			
			print(c(h, 'MARSS done'))
			
			ParEst[h,] = c(mu.reml, mu.reml2, MLpar$U, ssq.reml, ssq.reml2, MLpar$Q, tsq.reml, tsq.reml2, MLpar$R[1,1], beta0.reml, beta0.reml2, MLpar$x0)
		
		}
		
		write.table(ParEst, paste('ParEst_U_', mu0, '_Q_', round(ssq0,3), '_', N, 'obs.txt', sep=''))
		rm(ParEst, N)
	}
}

