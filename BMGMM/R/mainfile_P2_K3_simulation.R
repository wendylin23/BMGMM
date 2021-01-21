#################################################################
#################################################################
## Code to simulate and fit multivariate growth mixture models ##
## Wenyi Lin (wel316@ucsd.edu)						                     ##
## January, 2021          								                     ##
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
#########################
#########################
## fit with gmm.R code ##
#########################
#########################
n_class = K
nloop = 2000
burnin = 1000
thin = 1
sim = TRUE
Nsim = 100

gmm.fit =list()
K_true <- array(NA, dim = c(Nsim,N))

## simulate and fit growth mixture model
set.seed(rng_seed)
start_time <- Sys.time()
nsim = 1
idx = 1
while(nsim<=Nsim)
{
  Y_sparse = sim.dat[[nsim]]$Y_sparse
  time_sparse = sim.dat[[nsim]]$time_sparse
  M_i_sparse = sim.dat[[nsim]]$M_i_sparse
  X = sim.dat[[nsim]]$X
  Z = sim.dat[[nsim]]$Z
  pars_true = sim.dat[[nsim]]$pars_true
  tryCatch({
            print(nsim)
            sim.res =  mgmm(Y_sparse=Y_sparse,time_sparse=time_sparse,M_i_sparse=M_i_sparse,
                            X=X,Z=Z,n_class = n_class,P=P,
                            nloop=nloop,burnin=burnin,thin=thin,sim=TRUE,pars = pars_true) 
            gmm.fit[[nsim]] <- sim.res
            K_true[nsim,] <- sim.dat[[nsim]]$pars_true[[7]]
           },
           error = function(e){print(paste0("idx = ",idx," ",e))})
  idx = idx+1
  nsim = nsim+1
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
