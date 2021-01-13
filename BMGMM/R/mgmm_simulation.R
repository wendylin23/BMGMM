mgmm_simulation <- function(data.sim,pars_true,X,Z,plot=FALSE,
                            P=1,n_class=2,nloop=110,burnin=10,thin=1)
{
  library(MASS)
  library(mvtnorm)
  source("mgmm_var.R")
  beta_true = pars_true[[1]]
  alpha_true = pars_true[[2]]
  gamma_true = pars_true[[3]]
  #A_true = pars_true[[4]]
  Sigma_a_true = pars_true[[5]]
  s_sq_eps_true = pars_true[[6]]
  
  data.sim$K = NA
  K_true = NULL
  pi_true = NULL
  A_true = matrix(0, nrow = N,ncol = 2*P)
  for(i in 1:N){
    pi_i = Z[[i]]%*%gamma_true; pi_i = exp(pi_i)/apply(exp(pi_i),1,sum)
    pi_true = rbind(pi_true,pi_i)
    K_i = which(rmultinom(1, size=1, prob=pi_i)==1)
    K_true = c(K_true,K_i)
    #A_true[i,] = rmvnorm(1,rep(0,2*P),Sigma_a_true[,,K_i])
    A_true[i,] = rmvnorm(1,alpha_true[,K_i],Sigma_a_true[,,K_i])
    for(p in 1:P)
    {
      data.sim$K[data.sim$id == i&data.sim$outcome==p] = K_i
      T_i = cbind(1,data.sim$age[data.sim$id == i&data.sim$outcome==p])
      #A_true[i,(2*p-1):(2*p)] = A_true[i,(2*p-1):(2*p)] + cbind(alpha_true[(2*p-1):(2*p),K_i])
      data.sim$y[data.sim$id == i&data.sim$outcome==p] = X[[i]]%*%beta_true[,p] + T_i%*%(cbind(A_true[i,(2*p-1):(2*p)])) + rnorm(TT,0,s_sq_eps_true^.5)
    }
  }
  pars_true[[4]] = A_true
  pars_true[[7]] = K_true
  data.sim$keep = NA	# randomly delete (1-p)*100% of the observations to create sparse dataset
  keep_p = .75
  for(i in 1:N){
    data.sim$keep[data.sim$id == i] = rbinom(sum(data.sim$id == i&data.sim$outcome==1),1,keep_p)
    X[[i]] = rbind(X[[i]][which(data.sim$keep[data.sim$id == i&data.sim$outcome==1] == 1),])
  }
  data.sim = data.sim[data.sim$keep==1,]
  data.sim$keep = NULL
  
  #########################
  #########################
  ## Plot simulated data ##
  #########################
  #########################
  
  if(plot)
  {
    for(p in 1:P)
    {
      data.t = data.sim[data.sim$outcome==p,]
      spaghetti(y = data.t$y,times = data.t$age,ids = data.t$id,ylab="y",xlab="age", cols = K_true)
    } 
  }

  ###########################
  ###########################
  ## format data for input ##
  ###########################
  ###########################
  Y_sparse=list() # list for saving multivariate outcomes
  time_sparse=list() # list for saving time in random effects
  M_i_sparse = list() # list for saving number of inputs for each visit by outcome
  for(i in 1:N){
    Y_sparse[[i]]=list()
    time_sparse[[i]]=list()
    M_i_sparse[[i]]=rep(0,P)
    for(p in 1:P)
    {
      Y_sparse[[i]][[p]] = data.sim$y[data.sim$id == i&data.sim$outcome==p]
      time_sparse[[i]][[p]] = data.sim$age[data.sim$id == i&data.sim$outcome==p]
      M_i_sparse[[i]][p] = length(Y_sparse[[i]][[p]])
    }
  }
  
  gmm.fit = mgmm(Y_sparse=Y_sparse,time_sparse=time_sparse,M_i_sparse=M_i_sparse,
                 X=X,Z=Z,n_class = n_class,P=P,
                 nloop=nloop,burnin=burnin,thin=thin,sim=TRUE,pars = pars_true) 
  sim.dat = list(Y_sparse = Y_sparse,
                 time_sparse = time_sparse,
                 M_i_sparse = M_i_sparse,
                 X = X,
                 Z = Z,
                 pars_true = pars_true)
  return(list(gmm.fit = gmm.fit,sim.dat = sim.dat,K_true=K_true))
}