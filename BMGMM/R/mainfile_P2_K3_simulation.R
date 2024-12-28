#################################################################
#################################################################
## Code to simulate and fit multivariate growth mixture models ##
## Wenyi Lin (wel316@ucsd.edu)						                     ##
#################################################################
#################################################################

#rm(list = ls())
#setwd("~/Dropbox/UCSD/Projects/Integral Curve/Integral-Curve/R/Scratch")
args = commandArgs(F)
dir = "/home/wel316/simulation/20201208_chain1"
setwd(dir)
Sys.getenv("SGE_TASK_ID")

rng_seed <- 2020
set.seed(rng_seed)

###############
###############
## functions ##
###############
###############

source("mgmm_simulation.R")
spaghetti <- function(y,times,ids,ylab="y",xlab="times",
                      xlim=c(min(times,na.rm=T),max(times,na.rm=T)),
                      ylim=c(min(y,na.rm=T),max(y,na.rm=T)),
                      cols=rep(1,length(unique(ids))),lwd=1){
  plot(times,y,type="n",xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab)
  for(i in 1:length(unique(ids))){
    if(length(y[ids==unique(ids)[i]])==1) lines(times[ids==unique(ids)[i]],
                                                y[ids==unique(ids)[i]],type="p",col=cols[i],pch=15,cex=.25)
    if(length(y[ids==unique(ids)[i]])>1) lines(times[ids==unique(ids)[i]],
                                               y[ids==unique(ids)[i]],type="l",col=cols[i],lwd=lwd)
  }
}

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

####################
####################
## load libraries ##
####################
####################

library(MASS)
library(mvtnorm)

###########################
###########################
## simulation parameters ##
###########################
###########################

## simulate data
N <- 200  			# subjects
P <- 2  			# outcomes
TT <- 6     			# maximum time points
K <- 3     			# number of latent classes
age0 = rnorm(N, 0, 5) 	# baseline ages

## initiate simulated data
data.sim = data.frame(id = rep(1:N,each=TT*P),outcome = rep(rep(1:P,each=TT),N), visit = rep(1:TT,N*P), age0 = rep(age0,each=TT*P))
data.sim$time = NA
for(i in 1:N){
  for(p in 1:P)
  {
    data.sim$time[data.sim$id == i&data.sim$outcome==p] = c(0,sort(runif(TT-1, 0, 5)))  # TT number of time points
  }
}
data.sim$age = data.sim$age0 + data.sim$time

Z = list() # covariates latent class membership probabilities
X = list() # fixed effects covariates
for(i in 1:N){
  Z[[i]] = cbind(1,rnorm(1,0,1)) 
  ## only save the matrix once for fixed effects
  X[[i]] = as.matrix(cbind(rnorm(1,0,1),data.sim$age0[data.sim$id == i&data.sim$outcome==1]) )
}
Q_X = dim(X[[1]])[2]	# number of fixed-effect covariates for linear model
Q_Z = dim(Z[[1]])[2]	# number of fixed-effect covariates for multinomial logistic (class-membership) model 
Z_matrix = do.call(rbind, Z)

### set true parameter values & simulate data
s_sq_eps_true = .1 

# Dim of random effects: (2*P)*K
alpha_true = rbind(cbind(c(0.2,3),c(2,1),c(10,-0.5)), #p=1
                   cbind(c(0.5,1),c(1,2),c(2,-0.5)), #p=2
                   cbind(c(3,0.1),c(0.5,4),c(-2,2)))

# Dim of covariance of random effects: (2P)*(2P)
Sigma_a_true = rWishart(K,2*P,.2*diag(2*P))

# Dim of fixed effect: Q_X*P
beta_true = cbind(c(1,0.5),c(3,-0.1),c(4,-0.5))

# Dim of class coefficients: Q_Z*K
gamma_true = cbind(c(1,0.5),c(2,-0.5),c(0,0))
#A_true = rmvnorm(N,rep(0,2*P),Sigma_a_true)

# Initiate list for saving true parameters
pars_true = list()
pars_true[[1]] = beta_true
pars_true[[2]] = alpha_true
pars_true[[3]] = gamma_true
#pars_true[[4]] = A_true
pars_true[[5]] = Sigma_a_true
pars_true[[6]] = s_sq_eps_true
n_class = K
nloop = 2000
burnin = 1000
thin = 1
sim = TRUE
Nsim = 100

gmm.fit =list()
K_true <- array(NA, dim = c(Nsim,N))

## simulate and fit data
start_time <- Sys.time()
#for(nsim in 1:Nsim)
nsim = 1
idx = 1
while(nsim<=Nsim)
{
  tryCatch({
    print(nsim)
    sim.res = mgmm_simulation(data.sim = data.sim,pars_true = pars,X,Z,
                              P=P,n_class=K,nloop=nloop,burnin=burnin,thin=thin)
    gmm.fit[[nsim]] <- sim.res$gmm.fit
    K_true[nsim,] <- sim.res$pars_true[[7]]
    nsim = nsim+1
  },
  error = function(e){print(paste0("idx = ",idx," ",e))})
  idx = idx+1
  if((nsim-1)%%5==0){
    res_list = list(gmm.fit = gmm.fit,
                    pars_true = pars_true,
                    K_true = K_true
                    #time = end_time-start_time
    )
    save(res_list,file = paste0("results/mgmm_simulation_N",N,"_K",K,"_Outcome",P,".rdata"))
  }
}
end_time <- Sys.time()
print(end_time - start_time)

res_list = list(gmm.fit = gmm.fit,
                pars_true = pars_true,
                K_true = K_true,
                time = end_time-start_time)
save(res_list,file = paste0("results/mgmm_simulation_N",N,"_K",K,"_Outcome",P,".rdata"))

