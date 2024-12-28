#################################################################
#################################################################
## Code to simulate and fit univariate growth mixture models   ##
## Wenyi Lin (wel316@ucsd.edu)						                     ##
#################################################################
#################################################################
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

source("gmm_simulation.R")
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
rng_seed <- 2020
set.seed(rng_seed)
N <- 200  			# subjects
P <- 1    			# outcomes
TT <- 6     			# maximum time points
K <- 3  			# number of latent classes
age0 = rnorm(N, 0, 5) 	# baseline ages

data.sim = data.frame(id = rep(1:N,each=TT), 
                      visit = rep(1:TT,N), 
                      age0 = rep(age0,each=TT),
                      apoe = sample(c(0, 1), size = N, replace = TRUE))
data.sim$time = NA
for(i in 1:N){
  data.sim$time[data.sim$id == i] = c(0,sort(runif(TT-1, 0, 5)))  # TT number of time points
}
data.sim$age = data.sim$age0 + data.sim$time

Z = list() # covariates latent class membership probabilities
X = list() # fixed effects covariates
for(i in 1:N){
  Z[[i]] = cbind(1,rnorm(1,0,1)) 
  X[[i]] = as.matrix(cbind(rnorm(1,0,1),data.sim$age0[data.sim$id == i]) )
}
Q_X = dim(X[[1]])[2]	# number of fixed-effect covariates for linear model
Q_Z = dim(Z[[1]])[2]	# number of fixed-effect covariates for multinomial logistic (class-membership) model 
Z_matrix = do.call(rbind, Z)

### set true parameter values & simulate data
s_sq_eps_true = .1
alpha_true = cbind(c(0.2,3),c(2,1),c(10,-0.5))
Sigma_a_true = rWishart(K,2*P,.2*diag(2*P))
beta_true = cbind(c(1,0.5))
gamma_true = cbind(c(1,0.5),c(2,-0.5),c(0,0))
pars_true = list()
pars_true[[1]] = beta_true
pars_true[[2]] = alpha_true
pars_true[[3]] = gamma_true
pars_true[[5]] = Sigma_a_true
pars_true[[6]] = s_sq_eps_true

##################################
##################################
## simulation and model fitting ##
##################################
##################################
n_class = K
nloop = 2000
burnin = 1000
thin = 1
sim = TRUE
Nsim = 100
pars = pars_true
gmm.fit <- NULL
sim.dat = list()
K_true <- array(NA, dim = c(Nsim,N))

## simulate and fit growth mixture model
set.seed(rng_seed)
start_time <- Sys.time()
nsim = 1
idx = 1
while(nsim<=Nsim)
{
  tryCatch({
    print(nsim)
    sim.res <- gmm_simulation(data.sim = data.sim,pars_true = pars,X,Z,
                              P=P,n_class=K,nloop=nloop,burnin=burnin,thin=thin)
    gmm.fit[[nsim]] <- sim.res$gmm.fit
    sim.dat[[nsim]] <- sim.res$sim.dat
    K_true[nsim,] <- sim.res$K_true
    nsim = nsim+1
  },
  error = function(e){print(paste0("idx = ",idx," ",e))})
  idx = idx+1
  if((nsim-1)%%5==0){
    res_list = list(gmm.fit = gmm.fit,
                    sim.dat = sim.dat,
                    K_true = K_true
                    #time = end_time-start_time
    )
    save(res_list,file = paste0("results/mgmm_simulation_N",N,"_K",K,"_Outcome",P,".rdata"))
  }
}
end_time <- Sys.time()
print(end_time - start_time)

res_list = list(gmm.fit = gmm.fit,
                sim.dat = sim.dat,
                K_true = K_true,
                time = end_time-start_time)
save(res_list,file = paste0("results/mgmm_simulation_N",N,"_K",K,"_Outcome",P,".rdata"))
