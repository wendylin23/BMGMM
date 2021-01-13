mgmm=function(Y_sparse,time_sparse,M_i_sparse,X,Z,P=1,n_class=2,nloop=110,burnin=10,thin=1,sim=TRUE, pars = NULL){
  
  library(splines)
  library(MASS)
  library(Matrix)
  library(mvtnorm)
  library(MCMCpack)
  library(bayesm)
  #library(metaSEM)
  library(BayesLogit)
  library(fastDummies)
  
  #Y_sparse = Y_sparse
  #time_sparse = time_sparse
  #X = X
  #Z = Z
  #n_class = 2
  #nloop = 2100
  #burnin = 100
  #thin = 20
  #sim = TRUE
  #pars = pars_true
  
  N=length(Y_sparse)
  P=length(Y_sparse[[1]])
  id=1:N
  Q_X = dim(X[[1]])[2]
  Q_Z = dim(Z[[1]])[2]
  XX = t(X[[1]])%*%X[[1]]
  ZZ = t(Z[[1]])%*%Z[[1]]
  #M_i = length(Y_sparse[[1]][[1]])
  for(i in 2:N){
    XX = XX + t(X[[i]])%*%X[[i]]
    ZZ = ZZ + t(Z[[i]])%*%Z[[i]]
    #M_i = c(M_i, length(Y_sparse[[i]][[1]]))
  }
  M = sum(unlist(M_i_sparse))
  
  if(sim==FALSE){
    est_s_sq_eps = 1
    est_beta = 1
    est_A = 0
    est_Sigma_a = 1
    est_alpha = 1
    est_K = 1
    est_gamma = 1
    
    s_sq_eps = .1							# error variance
    A=array(0,dim=c(N,2)) 					# random effects for all N subjects
    Sigma_a = diag(rep(1,2*P))					# random effect variance
    alpha = array(0, dim = c(2*P,n_class))	# class-specific means of random effects
    beta = array(0, dim = c(Q_X,P))			# fixed effects	of linear model
    gamma = array(0, dim = c(Q_Z,n_class))	# multinomial coefficients for class-membership
    K = sample(1:n_class,N,replace=TRUE)	# indicators of class membership
  }
  if(sim==TRUE){
    est_beta = 1
    if(est_beta == 1) beta = pars[[1]]
    if(est_beta == 0) beta = array(0, dim = c(Q_X,P))	
    
    est_A = 1
    if(est_A == 1) A = pars[[4]]
    if(est_A == 0) A = array(0,dim=c(N,2*P))
    
    est_alpha = 1
    if(est_alpha == 1) alpha = pars[[2]]
    if(est_alpha == 0) alpha = array(0, dim = c(2*P,n_class))  #change the ncol and nrow
    
    est_Sigma_a = 1
    if(est_Sigma_a == 0) Sigma_a = pars[[5]][,,1]
    if(est_Sigma_a == 1) Sigma_a = .1*diag(rep(1,2*P))
    
    est_s_sq_eps = 1
    if(est_s_sq_eps == 1) s_sq_eps = pars[[6]]
    if(est_s_sq_eps == 0) s_sq_eps = .1	
    
    est_K =1
    if(est_K == 1) K = pars[[7]]
    if(est_K == 0) K = sample(1:n_class,N,replace=TRUE)
    
    est_gamma = 1
    if(est_gamma == 1) gamma = pars[[3]]
    if(est_gamma == 0) gamma = array(0, dim = c(Q_Z,n_class))
    #if(est_gamma == 1) gamma = pars[[3]]
  }		
  
  ## initialize parameters   
  
  cc=10000	# variance on hyper-prior for diffuse normal distributions
  eta=1; Psi = 1*diag(rep(1,2*P)) # hyper-parameters for IW prior on Sigma_b
  eta1=1; eta2=1 # hyper-parameters for IG prior on s_sq_eps
  nu = 100 # variance on hyper-prior for class parameter normal distributions
  
  array_ind=0
  
  A_mean = A
  beta_mean=beta
  alpha_mean=alpha
  gamma_mean=gamma
  Sigma_a_mean=Sigma_a
  s_sq_eps_mean=s_sq_eps
  K_mean = K
  
  # set up sampling arrays
  A_array=NULL
  BETA_array=NULL
  ALPHA_array=NULL
  GAMMA_array=NULL
  S_SQ_EPS_array=NULL
  SIGMA_A_array=NULL
  K_pos_array = NULL
  pi_pos_array = NULL
  K_array = matrix(0, nrow = N, ncol = n_class)
  for (i in 1:N) {
    K_array[i,K[i]] = 1
  }
  lambda = array(rep(diag(rep(1,N)),n_class-1),dim = c(N,N,n_class-1))
  C=matrix(rep(1,N*n_class),ncol = n_class,nrow = N)
  Z_array = matrix(0, nrow = N, ncol = n_class-1)
  for (k in 1:(n_class-1)) {
    z_t = abs(rlogis(N,0,1))
    for(i in 1:N)
    {
      if(K_array[i,k] == 1)
        #Z_array[i,k] = z_t[i] * ifelse(Z[[i]]%*%gamma[,k]>0,1,0)
        Z_array[i,k] = z_t[i]
      else
        #Z_array[i,k] = z_t[i] * ifelse(Z[[i]]%*%gamma[,k]<=0,1,0)
        Z_array[i,k] = -z_t[i]
    }
  }
  
  Z_all = do.call(rbind, Z)
  gamma_array = array(0,dim = c(Q_Z,n_class,nloop))
  #gamma_array = array(0,dim = c(Q_Z,nloop))
  #kappa = ifelse(K==1,0.5,-0.5)
  omega = matrix(rpg(N*n_class,1,0),ncol = n_class,nrow = N)
  #omega = rpg(N,1,0)
  
  # run MCMC 
  for(ii in 1:nloop){
    # if(ii%%thin==0)
    #   print(ii)
    
    #draw A
    if(est_A == 1){
      for(i in 1:N){
        #T_i = bdiagMat(lapply(time_sparse[[i]], function(x) cbind(1,x)))
        T_i = as.matrix(bdiag(lapply(time_sparse[[i]], function(x) cbind(1,x))))
        #X_t = matrix(rep(X[[i]][1,],10),ncol=2,byrow = TRUE)
        Y_i = cbind(unlist(Y_sparse[[i]])) - c(X[[i]]%*%beta)
        Sigma_a_i = solve(t(T_i)%*%T_i/s_sq_eps + solve(Sigma_a))
        mu_a_i = Sigma_a_i%*%(t(T_i)%*%Y_i/s_sq_eps +solve(Sigma_a)%*%cbind(alpha[,K[i]]))
        A[i,] = rmvnorm(1,mu_a_i,Sigma_a_i)
      }	
    }
    
    #draw alpha
    if(est_alpha == 1){
      for(k in 1:n_class){
        N_k = sum(K == k)
        Sigma_alpha_k = solve(N_k*solve(Sigma_a)+diag(rep(1,2*P))/cc)
        mu_alpha_k = Sigma_alpha_k%*%solve(Sigma_a)%*%cbind(apply(A[K==k,],2,sum))
        alpha[,k] = t(rmvnorm(1,mu_alpha_k,Sigma_alpha_k))
      }	
      alpha <- alpha[,order(alpha[1,])]
    }
    
    #draw Sigma_a
    if(est_Sigma_a == 1){
      eta_a = eta + N
      Psi_a = Psi + t(A-t(alpha[,K]))%*%(A-t(alpha[,K]))
      Sigma_a = riwish(eta_a, Psi_a)
    }
    
    #draw beta
    if(est_beta == 1){
      Sigma_beta = solve(XX/s_sq_eps + 1/cc*diag(1,dim(beta)[1]))
      for(p in 1:P)
      {
        mu_beta = cbind(rep(0,dim(beta)[1]))
        for(i in 1:N){
          T_i = cbind(1,time_sparse[[i]][[p]])
          Y_i = cbind(Y_sparse[[i]][[p]]) - T_i%*%(cbind(A[i,(2*p-1):(2*p)]))
          mu_beta = mu_beta +	t(X[[i]])%*%Y_i/s_sq_eps
        }	
        mu_beta = Sigma_beta%*%mu_beta
        beta[,p] = t(rmvnorm(1,mu_beta,Sigma_beta))
      }
    }
    
    #draw s_sq_eps
    if(est_s_sq_eps == 1){
      eta1_eps = 1 + M/2
      eta2_eps = 1
      for(i in 1:N){
        #T_i = bdiagMat(lapply(time_sparse[[i]], function(x) cbind(1,x)))
        T_i = as.matrix(bdiag(lapply(time_sparse[[i]], function(x) cbind(1,x))))
        Y_i = cbind(unlist(Y_sparse[[i]])) - cbind(rep((X[[i]][1,])%*%beta,M_i_sparse[[i]])) - T_i%*%(cbind(A[i,]))
        eta2_eps = eta2_eps + t(Y_i)%*%Y_i/2
      }	
      s_sq_eps = rinvgamma(1,eta1_eps,eta2_eps)			
    }
    
    #draw gamma (Holmes & Held 2006)
    # if(est_gamma == 1){
    #   K_array = as.matrix(dummy_cols(K)[,-1])
    #   for(k in 1:(n_class-1))
    #   {
    #     V = solve(t(Z_all)%*%solve(lambda[,,k])%*%Z_all + diag(1/nu,Q_Z))
    #     L = t(chol(V))
    #     B = V%*%t(Z_all)%*%solve(lambda[,,k])%*%(Z_array[,k]+log(C[,k]))
    #     temp = rnorm(Q_Z,0,1)
    #     gamma_array[,k,ii] = B + L%*%temp
    #     gamma[,k] = B + L%*%temp
    #     for(i in 1:N)
    #     {
    #       m = Z[[i]]%*%gamma_array[,k,ii]
    #       if(n_class>=3)
    #         C[i,k] = sum(exp(Z[[i]]%*%gamma_array[,-c(k),ii]))
    #       else
    #         C[i,k]=1
    #       Z_array[i,k] = sampleLogisticInd(m-log(C[i,k]),1,K_array[i,k]) ## Truncated logistic distribution
    #       R = as.numeric(Z_array[i,k]-m)
    #       lambda[i,i,k] = samplelambda(abs(R)) ## Algorithm A4: sample auxilliary variables \lambda
    #     }
    #   }
    # }
    
    #draw gamma (Polson, Scott and Windle 2013)
    
    # if(est_gamma == 1){
    #     #kappa = ifelse(K==1,0.5,-0.5)
    #     kappa = as.matrix(dummy_cols(K)[,-1])-0.5
    #     V = solve(t(Z_all)%*%diag(omega)%*%Z_all + diag(1/nu,Q_Z))
    #     m = V %*% (t(Z_all)%*% kappa)
    #     gamma_array[,1,ii] = rmvnorm(1, m, V)
    #     gamma[,1] = gamma_array[,1,ii]
    #     for(i in 1:N)
    #     {
    #       omega[i] = rpg(1,1,Z[[i]]%*%gamma_array[,1,ii])
    #     }
    # }
    if(est_gamma == 1){
      #kappa = ifelse(K==1,0.5,-0.5)
      kappa = as.matrix(dummy_cols(K)[,-1])-0.5
      for (k in 1:(n_class-1)) {
        for(i in 1:N)
        {
          if(n_class>=3)
            C[i,k] = sum(exp(Z[[i]]%*%gamma[,-k]))
          else
            C[i,k]=1
          omega[i,k] = rpg(1,1,Z[[i]]%*%gamma[,k]-log(C[i,k]))
        }
        V = solve(t(Z_all)%*%diag(omega[,k])%*%Z_all + diag(1/nu,Q_Z))
        m = V %*% (t(Z_all)%*%(kappa[,k]+diag(omega[,k])%*%log(C[,k])))
        gamma_array[,k,ii] = rmvnorm(1, m, V)
      }
      gamma = gamma_array[,,ii]
    }
    
    #draw K
    pi_pos = matrix(NA,nrow = N,ncol = n_class)
    if(est_K == 1){
      for(i in 1:N){
        pi_i = Z[[i]]%*%gamma;
        T_i = cbind(1,time_sparse[[i]][[1]])
        for(k in 1:n_class){
          pi_i[k] = pi_i[k] - .5*t(A[i,]-alpha[,k])%*%solve(Sigma_a)%*%(A[i,]-alpha[,k])
        }
        pi_i = pi_i - max(pi_i)
        pi_i = exp(pi_i)/apply(exp(pi_i),1,sum)
        K[i] = which(rmultinom(1,1:n_class,pi_i)==1)
        pi_pos[i,] = pi_i
        #tmp = K
        #if(mean(A[K==1,1])<mean(A[K==2,1])){
        #	tmp[K==1]=2
        #	tmp[K==2]=1
        #	K=tmp
        #}
      }
    }
    
    if(ii%%thin==0 &ii>=burnin ){
      array_ind=array_ind+1
      ALPHA_array[[array_ind]]=alpha
      BETA_array[[array_ind]]=beta
      GAMMA_array[[array_ind]]=gamma
      A_array[[array_ind]]=A
      SIGMA_A_array[[array_ind]] = Sigma_a
      S_SQ_EPS_array[[array_ind]] = s_sq_eps
      K_pos_array[[array_ind]]=K
      pi_pos_array[[array_ind]]=pi_pos
      
      alpha_mean=((array_ind-1)*alpha_mean+alpha)/array_ind
      beta_mean=((array_ind-1)*beta_mean+beta)/array_ind
      gamma_mean=((array_ind-1)*gamma_mean+gamma)/array_ind
      A_mean=((array_ind-1)*A_mean+A)/array_ind
      Sigma_a_mean=((array_ind-1)*Sigma_a_mean+Sigma_a)/array_ind
      s_sq_eps_mean=((array_ind-1)*s_sq_eps_mean+s_sq_eps)/array_ind
      K_mean=((array_ind-1)*K_mean+K)/array_ind
      
      # if(sim==TRUE){
      #   print(cbind(beta_mean,pars[[1]]))
      #   print(cbind(alpha_mean,pars[[2]]))
      #   print(cbind(gamma_mean,pars[[3]]))
      #   print(cbind(Sigma_a_mean,pars[[5]]))
      #   print(cbind(s_sq_eps_mean,pars[[6]]))
      # }
      # if(sim==FALSE){
      #   print(cbind(beta_mean))
      #   print(cbind(alpha_mean))
      #   print(cbind(gamma_mean))
      #   print(cbind(Sigma_a_mean))
      #   print(cbind(s_sq_eps_mean))
      # }	
    }
  }	
  
  mcmc_results=list()
  mcmc_results[[1]] = BETA_array
  mcmc_results[[2]] = ALPHA_array
  mcmc_results[[3]] = GAMMA_array
  mcmc_results[[4]] = A_array
  mcmc_results[[5]] = SIGMA_A_array
  mcmc_results[[6]] = S_SQ_EPS_array
  mcmc_results[[7]] = K_pos_array
  mcmc_results[[8]] = pars
  mcmc_results[[9]] = beta_mean
  mcmc_results[[10]] = alpha_mean
  mcmc_results[[11]] = gamma_mean
  mcmc_results[[12]] = A_mean
  mcmc_results[[13]] = Sigma_a_mean
  mcmc_results[[14]] = s_sq_eps_mean
  mcmc_results[[15]] = K_mean
  
  #par(mfrow=c(2,1))
  #plot(A_true[,1],A_mean[,1])
  #plot(A_true[,2],A_mean[,2])
  #cor(A_true[,1],A_mean[,1])
  #cor(A_true[,2],A_mean[,2])
  #sum(round(K_mean) == K_true)
  
  return(mcmc_results)
}

samplelambda <- function(r)
{
  ok = 0
  while (ok!=1) {
    Y = rnorm(1,0,1)
    Y = Y*Y
    Y = 1+(Y-sqrt(Y*(4*r+Y)))/(2*r)
    U = runif(1,0,1)
    if(U<=1/(1+Y))
      lambda = r/Y
    else
      lambda = r*Y
    U = runif(1,0,1)
    if(lambda > 4/3)
      ok = rightmost_interval(U,lambda,ok)
    else
      ok = leftmost_interval(U,lambda,ok)
    end
  }
  return(lambda)
}

rightmost_interval <- function(U,lambda,ok)
{
  Z = 1
  X = exp(-.5*lambda)
  j = 0
  
  j = j + 1
  Z = Z-((j+1)^2)*X^(((j+1)^2) - 1)
  
  if (Z > U) 
    ok = 1
  else
  {
    j = j + 1
    Z = Z+((j+1)^2)*X^(((j+1)^2) - 1)
    if(Z < U)
      ok = 0
  }
  return(ok)
}

leftmost_interval <- function(U,lambda,ok)
{
  H = 0.5*log(2) + 2.5*log(pi) - 2.5*log(lambda) - (pi^2)/(2*lambda) + 0.5*lambda
  lU = log(U)
  Z = 1
  X = exp((-pi^2)/(2*lambda))
  K = lambda/pi^2
  j = 0
  
  j = j + 1
  Z = Z - K*X^((j^2)-1)
  
  if (H + log(Z) > lU) 
    ok = 1
  else
  {
    j = j + 1
    Z = Z + ((j+1)^2)*X^(((j+1)^2) -1)
    if(H + log(Z) < lU)
      ok = 0
  }
  return(ok)
}

logisticInversecdf = function(p,a,b)
{
  a + b*(log(p)-log(1-p))
}

sampleLogisticInd = function(m,v,y)
{
  s = sign(y-0.5)
  U = runif(1,0,1)
  sample = logisticInversecdf(U,m,v)
  while(sign(sample)!=s)
  {
    if(s > 0)
    {
      U = 1-runif(1,0,1)*(1-U)
      sample = logisticInversecdf(U,m,v)
    }
    else
    {
      U = runif(1,0,1)*U
      sample = logisticInversecdf(U,m,v)
    }
  }
  return(sample)
}