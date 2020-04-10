# Author: Kevin See
# Purpose: Create figures for PVA manuscript: "Reducing bias and improving precision in species extinction forecasts"
# Created: 10/9/14


library(ggplot2)
library(gridExtra)
library(reshape2)
library(plyr)

setwd('/Users/kevin/Dropbox/Personal/Ph_D/Chapters/Ch_1/TimeVsTotObs')

# this computes the coefficient of variation, based on a known "mean" value
cv = function(x, y){
	sd(x, na.rm=T) / abs(y)
}

# this calculates the estimated prob of quasi-extict based on estimates of mu and process variance
# returns the probability of reaching the specified percent decline at some point within the specified time horizon, given the trend and process error variance provided
# x: a vector of parameter estimates: first trend (mu), second process variance (sigma^2)
# perDec: the quasi-extinction threhold (e.g. 0.3 = risk of a 30% decline)
# Th: time horizon over which we are projecting this risk (e.g. 30 years)
ProbExt = function(x, perDec, Th){
	a= -log(1-perDec)
	u = as.numeric(-x[1] * sqrt(Th) / sqrt(x[2]))
	v = as.numeric(a / (sqrt(x[2]) * sqrt(Th)))
	pe = pnorm(u-v) + exp(2*u*v)*pnorm(-(u+v))
	Pe = ifelse(exp(2*u*v) == Inf, 0, pe)
	return(Pe)
}

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
Qs = sig		# Qs are variances
Rs = tau		# Rs are variances
Us = c(-0.02, -0.04) # Us are mean rates of growth
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
# Read in parameter estimates
#---------------------------------------------
ParEst = NULL
for(i in 1:length(Us)){
	for(j in 1:4){
		# read estimates back in
		temp = read.table(paste('Estimates//U',Us[i], '_', j,'.txt', sep=''), header=T)
		ParEst = rbind(ParEst, temp)
	}
}
ParEst$Q.R.ratio = as.factor(with(ParEst, Q2/R2))
levels(ParEst$Q.R.ratio) = paste('sigma^2 / tau^2 ==', sig2Tau)

param.est.df = melt(ParEst, id.vars=c('TotObs', 'TSteps', 'U', 'Q2', 'R2', 'Q.R.ratio'), measure.vars=c('EstU', 'EstQ.biased', 'EstR.biased'), variable.name='Parameter', value.name='Estimate')
levels(param.est.df$Parameter) = c('U', 'Q', 'R')
param.est.df$True.Value = NA
param.est.df$True.Value[param.est.df$Parameter=='U'] = param.est.df$U[param.est.df$Parameter=='U']
param.est.df$True.Value[param.est.df$Parameter=='Q'] = param.est.df$Q2[param.est.df$Parameter=='Q']
param.est.df$True.Value[param.est.df$Parameter=='R'] = param.est.df$R2[param.est.df$Parameter=='R']
param.est.df$Abs.Error = with(param.est.df, Estimate - True.Value)
param.est.df$Rel.Error = with(param.est.df, (Estimate - True.Value) / True.Value)
param.est.df$plot.group = with(param.est.df, as.factor(paste(TotObs, TSteps, sep='.')))
for(t.step in tSteps) {
	for(tot.obs in totobs) param.est.df$plot.group = relevel(param.est.df$plot.group, ref=paste(tot.obs, t.step, sep='.'))
}
levels(param.est.df$Parameter) = c('mu', 'sigma^2', 'tau^2')
param.est.df$U = as.factor(param.est.df$U)
levels(param.est.df$U) = paste('mu ==', rev(Us))

#---------------------------------------------
# Assess the accuracy of parameter estimates
#---------------------------------------------
# focus on trend parameter
mu.err.p = ggplot(subset(param.est.df, Parameter=='mu'), aes(y=Abs.Error, x=as.factor(TSteps), group=plot.group)) +
  geom_boxplot(aes(fill=log(TotObs/TSteps)), outlier.size=0.5) +
  geom_hline(yintercept=0) +
  theme_bw() + 
  theme(legend.position='bottom') +
  scale_fill_gradient2(low='blue', mid='white', high='red', midpoint=0, space='Lab', name='Obs / Yr', labels=c(0.5, 1, 2, 6, 18), breaks=log(c(0.5, 1, 2, 6, 18))) +
  facet_grid(Q.R.ratio~U, scales='free', labeller=label_parsed) +
  labs(y='Absolute Error', x='Length of Observation Period', title=expression(paste('Esimtates of ', mu)))

print(mu.err.p)

# focus on variance estimates (Q and R)
var.est.p1 = ggplot(subset(param.est.df, Parameter!='mu' & U==levels(param.est.df$U)[2]), aes(y=Abs.Error, x=as.factor(TSteps), group=plot.group)) +
  geom_boxplot(aes(fill=log(TotObs/TSteps)), outlier.size=0.5) +
  geom_hline(yintercept=0) +
  theme_bw() + 
  theme(legend.position='bottom') +
  scale_fill_gradient2(low='blue', mid='white', high='red', midpoint=0, space='Lab', name='Obs / Yr', labels=c(0.5, 1, 2, 6, 18), breaks=log(c(0.5, 1, 2, 6, 18))) +
  facet_grid(Q.R.ratio~Parameter, scales='free', labeller=label_parsed) +
  labs(y='Absolute Error', x='Length of Observation Period', title=expression(atop('Variance Estimates', atop(paste(mu, '= -0.02'), ''))))

var.est.p2 = ggplot(subset(param.est.df, Parameter!='mu' & U==levels(param.est.df$U)[1]), aes(y=Abs.Error, x=as.factor(TSteps), group=plot.group)) +
  geom_boxplot(aes(fill=log(TotObs/TSteps)), outlier.size=0.5) +
  geom_hline(yintercept=0) +
  theme_bw() + 
  theme(legend.position='bottom') +
  scale_fill_gradient2(low='blue', mid='white', high='red', midpoint=0, space='Lab', name='Obs / Yr', labels=c(0.5, 1, 2, 6, 18), breaks=log(c(0.5, 1, 2, 6, 18))) +
  facet_grid(Q.R.ratio~Parameter, scales='free', labeller=label_parsed) +
  labs(y='Absolute Error', x='Length of Observation Period', title=expression(atop('Variance Estimates', atop(paste(mu, '= -0.04'), ''))))

print(var.est.p1)
print(var.est.p2)

#---------------------------------------------
# Assess the precision of parameter estimates
#---------------------------------------------
# estimate precision of parameter estimates at each time step / total obs combination
param.summ = ddply(param.est.df, .(Parameter, U, Q.R.ratio, TSteps, TotObs), summarise,
  CV.Estimate = sd(Estimate) / abs(mean(True.Value)))

# smooth the contours for plotting purposes
my.grid = with(param.summ, expand.grid(list(TSteps=seq(min(TSteps), max(TSteps),1), TotObs=seq(min(TotObs), max(TotObs), 1))))
plot.param.prec = NULL
for(my.u in unique(param.summ$U)) {
	for(my.q.r in unique(param.summ$Q.R.ratio)) {
		for(par in levels(param.summ$Parameter)) {
			my.subset = subset(param.summ, Parameter==par & U==my.u & Q.R.ratio==my.q.r)
			my.loess = loess(CV.Estimate ~ TSteps + TotObs, data=my.subset, degree=1, span=0.9)
			my.predict = predict(my.loess, newdata=my.grid)
			plot.data = melt(my.predict)
			plot.data$Parameter = par
			plot.data$U = unique(my.subset$U)
			plot.data$Q.R.ratio = as.factor(unique(my.subset$Q.R.ratio))
			plot.data$TSteps = gsub('TSteps= ', '', plot.data$TSteps)
			plot.data$TSteps = gsub('TSteps=', '', plot.data$TSteps)
			plot.data$TSteps = as.numeric(plot.data$TSteps)
			plot.data$TotObs = gsub('TotObs= ', '', plot.data$TotObs)
			plot.data$TotObs = gsub('TotObs=', '', plot.data$TotObs)
			plot.data$TotObs = as.numeric(plot.data$TotObs)			
			plot.param.prec = rbind(plot.param.prec, plot.data)
		}
	}
}
plot.param.prec$Parameter = factor(plot.param.prec$Parameter, levels=c('mu', 'sigma^2', 'tau^2'))
levels(plot.param.prec$Q.R.ratio) = paste('sigma^2 / tau^2 ==', sig2Tau)

param.prec.p1 = ggplot(subset(plot.param.prec, U==levels(plot.param.prec$U)[2]), aes(x=TSteps, y=TotObs, z=value)) + 
  geom_tile(aes(fill=value)) +
  stat_contour(aes(color=..level..), bins=10) +
  scale_color_gradient(low='gray80', high='gray20', guide=F) +
  scale_fill_gradient(low='white', high='black', name='CV of\nEstimate') +
  facet_grid(Parameter~Q.R.ratio, scales='free', labeller=label_parsed) +
  geom_abline(intercept=0, slope=1, lty=2, col='black') +
  theme_bw() +
  theme(legend.position='bottom', legend.key.width=unit(1.5, units='cm')) +
  labs(y='Total Number of Observations', x='Length of Observation Period', title=expression(atop('Precision of Parameter Estimates', atop(paste(mu, '= -0.02'), ''))))

param.prec.p2 = ggplot(subset(plot.param.prec, U==levels(plot.param.prec$U)[1]), aes(x=TSteps, y=TotObs, z=value)) + 
  geom_tile(aes(fill=value)) +
  stat_contour(aes(color=..level..), bins=10) +
  scale_color_gradient(low='gray80', high='gray20', guide=F) +
  scale_fill_gradient(low='white', high='black', name='CV of\nEstimate') +
  facet_grid(Parameter~Q.R.ratio, scales='free', labeller=label_parsed) +
  geom_abline(intercept=0, slope=1, lty=2, col='black') +
  theme_bw() +
  theme(legend.position='bottom', legend.key.width=unit(1.5, units='cm')) +
  labs(y='Total Number of Observations', x='Length of Observation Period', title=expression(atop('Precision of Parameter Estimates', atop(paste(mu, '= -0.04'), ''))))

print(param.prec.p1)
print(param.prec.p2)

#---------------------------------------------
# Examine estimates of probability of quasi-extinction
#---------------------------------------------
# create data frame with Pe estimates
# set time horizon
time.horiz = 30
PE.df.30 = NULL
for(perc.decl in c(0.3, 0.5, 0.8)) PE.df.30 = rbind(PE.df.30, data.frame(ParEst, Perc.Decline = perc.decl*100, True.Prob = apply(ParEst[,1:2],1,ProbExt, perDec=perc.decl, Th=time.horiz), Est.Prob = apply(ParEst[,6:7],1,ProbExt, perDec=perc.decl, Th=time.horiz)))
PE.df.30$Abs.Error = with(PE.df.30, Est.Prob - True.Prob)
PE.df.30$Rel.Error = with(PE.df.30, (Est.Prob - True.Prob) / True.Prob)
PE.df.30$plot.group = with(PE.df.30, as.factor(paste(TotObs, TSteps, sep='.')))
for(t.step in tSteps) {
  for(tot.obs in totobs) PE.df.30$plot.group = relevel(PE.df.30$plot.group, ref=paste(tot.obs, t.step, sep='.'))
}
levels(PE.df.30$Q.R.ratio) = paste('sigma^2 / tau^2 ==', sig2Tau)
PE.df.30$Perc.Decl.label = as.factor(PE.df.30$Perc.Decline)

# change time horizon
time.horiz = 50
PE.df.50 = NULL
for(perc.decl in c(0.3, 0.5, 0.8)) PE.df.50 = rbind(PE.df.50, data.frame(ParEst, Perc.Decline = perc.decl*100, True.Prob = apply(ParEst[,1:2],1,ProbExt, perDec=perc.decl, Th=time.horiz), Est.Prob = apply(ParEst[,6:7],1,ProbExt, perDec=perc.decl, Th=time.horiz)))
PE.df.50$Abs.Error = with(PE.df.50, Est.Prob - True.Prob)
PE.df.50$Rel.Error = with(PE.df.50, (Est.Prob - True.Prob) / True.Prob)
PE.df.50$plot.group = with(PE.df.50, as.factor(paste(TotObs, TSteps, sep='.')))
for(t.step in tSteps) {
  for(tot.obs in totobs) PE.df.50$plot.group = relevel(PE.df.50$plot.group, ref=paste(tot.obs, t.step, sep='.'))
}
levels(PE.df.50$Q.R.ratio) = paste('sigma^2 / tau^2 ==', sig2Tau)
PE.df.50$Perc.Decl.label = as.factor(PE.df.50$Perc.Decline)

# combine both time horizons
PE.df = rbind(data.frame(PE.df.30, Time.Horiz=30), data.frame(PE.df.50, Time.Horiz=50)) 
PE.df$Perc.Decl.label = as.factor(PE.df$Perc.Decline)
levels(PE.df$Perc.Decl.label) = paste('Perc.Decl. ==', unique(PE.df$Perc.Decline))

# Accuracy, measured by absolute error
pe.acc.box.p = vector('list', 4)
i = 1
for(time.horiz in unique(PE.df$Time.Horiz)) {
  for(my.u in unique(PE.df$U)) {
    pe.acc.box.p[[i]] = ggplot(subset(PE.df, Time.Horiz==time.horiz & U==my.u), aes(y=Abs.Error, x=as.factor(TSteps), group=plot.group)) +
      geom_boxplot(aes(fill=log(TotObs/TSteps)), outlier.size=0.5) +
      geom_hline(yintercept=0) +
      scale_fill_gradient2(low='blue', mid='white', high='red', midpoint=0, space='Lab', name='Obs / Yr', labels=c(0.5, 1, 2, 6, 18), breaks=log(c(0.5, 1, 2, 6, 18))) +
      facet_grid(Perc.Decl.label~Q.R.ratio, scales='free', labeller=label_parsed) +
      theme_bw() +
      labs(y='Absolute Error', x='Length of Observation Period', title=bquote(atop(paste('Estimating a Decline over ', .(time.horiz), ' Years'), atop(paste(mu, '=', .(my.u)), ''))))
    
    i = i+1
  }
}

pe.acc.box.30 = arrangeGrob(pe.acc.box.p[[1]], pe.acc.box.p[[2]], ncol=1)
pe.acc.box.50 = arrangeGrob(pe.acc.box.p[[3]], pe.acc.box.p[[4]], ncol=1)

print(pe.acc.box.30)
print(pe.acc.box.50)

# summarize PE with median to interpolate contour plot
pe.summ = ddply(PE.df, .(U, Q.R.ratio, Time.Horiz, Perc.Decline, TotObs, TSteps), summarise,
  CV.pe = sd(Est.Prob) / mean(True.Prob),
  Abs.Error = median(Abs.Error),
  Rel.Error = median(Rel.Error))
  
  
# smooth the contours for plotting purposes
my.grid = with(PE.df, expand.grid(list(TSteps=seq(min(TSteps), max(TSteps),1), TotObs=seq(min(TotObs), max(TotObs), 1))))
my.form = as.formula('Abs.Error ~ TSteps + TotObs')
plot.pe.acc = NULL
for(time.horiz in unique(PE.df$Time.Horiz)) {
  for(my.u in unique(PE.df$U)) {
  	for(my.q.r in unique(PE.df$Q.R.ratio)) {
  	  for(perc.decl in unique(PE.df$Perc.Decline)) {
        my.subset = subset(pe.summ, U==my.u & Q.R.ratio==my.q.r & Perc.Decline==perc.decl & Time.Horiz==time.horiz)  
  			my.loess = loess(my.form, data=my.subset, degree=1, span=0.5)
  			my.predict = predict(my.loess, newdata=my.grid)
  			plot.data = melt(my.predict)
        plot.data$Time.Horiz = time.horiz
  			plot.data$U = my.u
  			plot.data$Q.R.ratio = as.factor(unique(my.subset$Q.R.ratio))
  			plot.data$Perc.Decline = perc.decl
  			plot.data$TSteps = gsub('TSteps= ', '', plot.data$TSteps)
  			plot.data$TSteps = gsub('TSteps=', '', plot.data$TSteps)
  			plot.data$TSteps = as.numeric(plot.data$TSteps)
  			plot.data$TotObs = gsub('TotObs= ', '', plot.data$TotObs)
  			plot.data$TotObs = gsub('TotObs=', '', plot.data$TotObs)
  			plot.data$TotObs = as.numeric(plot.data$TotObs)			
  			plot.pe.acc = rbind(plot.pe.acc, plot.data)
  		}
  	}
  }
}
plot.pe.acc$Perc.Decl.label = as.factor(plot.pe.acc$Perc.Decline)
levels(plot.pe.acc$Perc.Decl.label) = paste('Perc.Decl. ==', unique(plot.pe.acc$Perc.Decline))

pe.acc.p = vector('list', 4)
i = 1
for(time.horiz in unique(plot.pe.acc$Time.Horiz)) {
  for(my.u in unique(plot.pe.acc$U)) {
    pe.acc.p[[i]] = ggplot(subset(plot.pe.acc, U==my.u & Time.Horiz==time.horiz), aes(x=TSteps, y=TotObs, z=value)) + 
      geom_tile(aes(fill=value)) +
      # stat_contour(aes(color=..level..), bins=10) +
      # scale_color_gradient2(low='darkblue', mid='gray', high='darkred', guide=F) +
      scale_fill_gradient2(low='blue', mid='white', high='red', name='Absolute\nError') +
      facet_grid(Perc.Decl.label~Q.R.ratio, scales='free', labeller=label_parsed) +
      # geom_abline(intercept=0, slope=1, lty=2, col='black') +
      # geom_abline(intercept=0, slope=2, lty=2, col='darkred') +
      theme_bw() +
      theme(legend.position='bottom') +
      labs(y='Total Number of Observations', x='Length of Observation Period', title=bquote(atop(paste('Probability of Quasi-Extinction in ', .(time.horiz), ' Years'), atop(paste(mu, ' = ', .(my.u)), ''))))
    
    i = i+1
  }
}

pe.acc.30 = arrangeGrob(pe.acc.p[[1]], pe.acc.p[[2]], ncol=1)
pe.acc.50 = arrangeGrob(pe.acc.p[[3]], pe.acc.p[[4]], ncol=1)

print(pe.acc.30)
print(pe.acc.50)

#----------------------------------------
# Look at percent improvement 
#----------------------------------------
# use reference case of 30 years, 30 total observations
ref.case = subset(pe.summ, TotObs==30 & TSteps==30)
pe.ref.all = NULL
for(my.u in unique(pe.summ$U)) {
	for(my.q.r in unique(pe.summ$Q.R.ratio)) {
		pe.ref = subset(pe.summ, U==my.u & Q.R.ratio==my.q.r, c('Time.Horiz', 'Perc.Decline', 'CV.pe'))
		my.ref = subset(ref.case, U==my.u & Q.R.ratio==my.q.r, c('Time.Horiz', 'Perc.Decline', 'CV.pe'))
		for(perc.decl in unique(pe.summ$Perc.Decline)) {
          for(time.horiz in unique(pe.summ$Time.Horiz)) {
            pe.ref$CV.pe[with(pe.ref, Time.Horiz==time.horiz & Perc.Decline==perc.decl)] = 100 * (pe.ref$CV.pe[with(pe.ref, Time.Horiz==time.horiz & Perc.Decline==perc.decl)] - my.ref$CV.pe[with(my.ref, Time.Horiz==time.horiz & Perc.Decline==perc.decl)]) / my.ref$CV.pe[with(my.ref, Time.Horiz==time.horiz & Perc.Decline==perc.decl)]
          }
        }
		pe.ref.all = rbind(pe.ref.all, pe.ref)
	}
}
pe.ref.all = cbind(pe.summ[,c(1,2,5,6)], pe.ref.all)


my.grid = with(pe.ref.all, expand.grid(list(TSteps=seq(min(TSteps), max(TSteps),1), TotObs=seq(min(TotObs), max(TotObs), 1))))
my.form = as.formula(CV.pe ~ TSteps + TotObs)
plot.data.pe.imp = NULL
for(my.u in unique(pe.ref.all$U)) {
	for(my.q.r in unique(pe.ref.all$Q.R.ratio)) {
	  for(perc.decl in unique(pe.ref.all$Perc.Decline)) {
	    for(time.horiz in unique(pe.ref.all$Time.Horiz)) {
        my.subset = subset(pe.ref.all, U==my.u & Q.R.ratio==my.q.r & Perc.Decline==perc.decl & Time.Horiz==time.horiz)
			
  			my.loess = loess(my.form, data=my.subset, degree=2, span=0.7)
  
  			my.predict = predict(my.loess, newdata=my.grid)
  			plot.data = melt(my.predict)
  			plot.data$U = unique(my.subset$U)
  			plot.data$Q.R.ratio = as.factor(unique(my.subset$Q.R.ratio))
  			plot.data$Perc.Decline = perc.decl
        plot.data$Time.Horiz = time.horiz
  			plot.data$TSteps = gsub('TSteps= ', '', plot.data$TSteps)
  			plot.data$TSteps = gsub('TSteps=', '', plot.data$TSteps)
  			plot.data$TSteps = as.numeric(plot.data$TSteps)
  			plot.data$TotObs = gsub('TotObs= ', '', plot.data$TotObs)
  			plot.data$TotObs = gsub('TotObs=', '', plot.data$TotObs)
  			plot.data$TotObs = as.numeric(plot.data$TotObs)			
  			plot.data.pe.imp = rbind(plot.data.pe.imp, plot.data)
	    }
		}
	}
}
plot.data.pe.imp$Perc.Decl.label = as.factor(plot.data.pe.imp$Perc.Decline)
levels(plot.data.pe.imp$Perc.Decl.label) = paste('Perc.Decl. ==', unique(plot.data.pe.imp$Perc.Decline))
plot.data.pe.imp$Q.R.ratio.simp = plot.data.pe.imp$Q.R.ratio
levels(plot.data.pe.imp$Q.R.ratio.simp) = sig2Tau

# make 2 versions of plots
pe.imp.p.col = pe.imp.p.bw.1 = vector('list', 4)
# pick one percent decline to focus on
perc.decl = unique(plot.data.pe.imp$Perc.Decline)[2]
i = 1
for(time.horiz in unique(plot.pe.acc$Time.Horiz)) {
  for(my.u in unique(plot.pe.acc$U)) {
    pe.imp.p.col[[i]] = ggplot(subset(plot.data.pe.imp, U==my.u & Time.Horiz==time.horiz), aes(x=TSteps, y=TotObs, z=value)) + 
      geom_tile(aes(fill=value)) +
      stat_contour(aes(color=..level..), bins=10) +
      geom_point(aes(x=30, y=30), pch=17) +
      scale_color_gradient2(low='darkblue', mid='gray', high='darkred', guide=F) +
      scale_fill_gradient2(low='blue', mid='white', high='red', name='Percent\nChange\n in CV') +
      facet_grid(Perc.Decl.label~Q.R.ratio, scales='free', labeller=label_parsed) +
      geom_abline(intercept=0, slope=1, lty=2, col='black') +
      geom_abline(intercept=0, slope=2, lty=3, col='black') +
      theme_bw() +
      theme(legend.position='right', legend.key.height=unit(1.5, units='cm')) +
      labs(y='Total Number of Observations', x=bquote(atop('Length of Observation Period', paste(mu, ' = ', .(my.u)))), title=bquote(atop('Difference in Precision of the Probability',  paste('of Quasi-Extinction in ', .(time.horiz), ' Years'))))
 
    pe.imp.p.bw.1[[i]] = ggplot(subset(plot.data.pe.imp, U==my.u & Time.Horiz==time.horiz & Perc.Decline==perc.decl), aes(x=TSteps, y=TotObs, z=abs(value))) + 
      geom_tile(aes(fill=abs(value))) +
      geom_point(aes(x=30, y=30), pch=17) +
      scale_fill_gradient(low='white', high='black', name='Percent\nDifference in CV') +
      facet_grid(Q.R.ratio~., labeller=label_parsed) +
      geom_abline(intercept=0, slope=1, lty=2, col='black') +
      geom_abline(intercept=0, slope=2, lty=3, col='black') +
      theme_bw() +
      theme(legend.position='right') +
      labs(y='Total Number of Observations', x=bquote(atop('Length of Observation Period', paste(mu, ' = ', .(my.u)))), title=bquote(atop('Difference in Precision of the Probability',  paste('of a ', .(perc.decl), '% Decline over ', .(time.horiz), ' Years'))))
    
    i = i+1
  }
}

print(arrangeGrob(pe.imp.p.col[[1]], pe.imp.p.col[[2]]))
print(arrangeGrob(pe.imp.p.col[[3]], pe.imp.p.col[[4]]))
print(pe.imp.p.bw.1[[2]])


#---------------------------------------------
# Estimation comparison (EM, REML 1 & REML 2)
#---------------------------------------------
# read in ALL estimates
mu0=Us[2]
variables = c('U', 'Q', 'R', 'X0')
par.est.methods = NULL
for(n in 1:3){
  N=n
  for(j in 1:4){
    ssq0=Qs[j]
    tsq0=Rs[j]
    # read estimates back in
    temp = read.table(paste('REML_Test/ParEst_U_', mu0, '_Q_', round(ssq0, 3), '_', N, 'obs.txt', sep=''), header=T)
    temp = cbind(temp, N, mu0, ssq0, tsq0)
    par.est.methods = rbind(par.est.methods, temp)
  }
}

plot.data.methods = melt(par.est.methods, id.vars=c('N', 'mu0', 'ssq0', 'tsq0'))
names(plot.data.methods) = c('N', 'U', 'Q', 'R', 'variable', 'Estimate')
plot.data.methods$Parameter = as.factor(sapply(strsplit(as.character(plot.data.methods$variable), '[.]'), function(x) x[1]))
plot.data.methods$Method = sapply(strsplit(as.character(plot.data.methods$variable), '[.]'), function(x) x[2])
plot.data.methods$Method = as.factor(plot.data.methods$Method)
levels(plot.data.methods$Method) = c('EM', 'REML.1', 'REML.2')
plot.data.methods$Q.R = as.factor(with(plot.data.methods, Q/R))
plot.data.methods$True = NA
plot.data.methods$True[which(plot.data.methods$Parameter=='X0')] = x0
for(variable in c('U', 'Q', 'R')) plot.data.methods$True[which(plot.data.methods$Parameter==variable)] = plot.data.methods[which(plot.data.methods$Parameter== variable), variable]
levels(plot.data.methods$Parameter) = c('sigma^2', 'tau^2', 'mu', 'x_0')
plot.data.methods$Parameter = factor(plot.data.methods$Parameter, levels=c('mu', 'sigma^2', 'tau^2', 'x_0'))
levels(plot.data.methods$Q.R) = paste('sigma^2 / tau^2 ==', sig2Tau)
plot.data.methods$Rel.Error = with(plot.data.methods, (Estimate-True) / True)
plot.data.methods$Abs.Error = with(plot.data.methods, Estimate-True)

methods.p = ggplot(subset(plot.data.methods, Parameter!='x_0' & abs(Rel.Error)<25), aes(x=as.factor(N), y=Rel.Error)) +
  geom_boxplot(aes(fill=Method)) +
  geom_hline(yintercept=0) +
  facet_grid(Parameter ~ Q.R, labeller=label_parsed, scales='free') +
  theme_bw() +
  scale_fill_grey(start=0.2, end=0.9) + 
  labs(x='Number of Observations / Year', y='Relative Error', title='Comparison of Methods') +
  theme(legend.position='bottom')

print(methods.p)
