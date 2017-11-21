# Simulate data with JAGS

Recently, I have been struggling with simulating data from complex hierarchical models. After several unsuccessful attempts in `R`, I remembered the good old times when I was using `WinBUGS` (more than 10 years already!) and the possibility to simulate data with it. I'm using `Jags` now, and a quick search in Google with 'simulating data with jags' led me to [a complex example](https://www.georg-hosoya.de/wordpress/?p=799) and [a simple example](https://stackoverflow.com/questions/38295839/simulate-data-in-jags-r2jags).

Here, I illustrate the possibility to use Jags to simulate data with two examples that might be of interest to population ecologists: first a linear regression, second a Cormack-Jolly-Seber capture-recapture model to estimate animal survival (formulated as a state-space model). 

Simulating data with `Jags` is convenient because you can use (almost) the same code for simulation and inference, and you can carry out simulation studies (bias, precision, interval coverage) in the same environment (namely `Jags`). 
