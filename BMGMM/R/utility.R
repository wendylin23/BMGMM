########################################
##### Useful functions #################
########################################

# Extracting fitting results from MGMM fitting
extract_mgmm_res <- function(fit.res, true_val, K)
{
  BETA_array <- do.call(rbind,lapply(fit.res[[1]],as.vector))
  ALPHA_array <- do.call(rbind,lapply(fit.res[[2]],as.vector))
  GAMMA_array <- do.call(rbind,lapply(fit.res[[3]],as.vector))[,1:(2*(K-1))]
  SIGMA_A_array <- do.call(rbind,lapply(fit.res[[5]],as.vector))
  K_est <- fit.res[[15]]
  beta_mean <- as.vector(fit.res[[9]])
  alpha_mean <- as.vector(fit.res[[10]])
  gamma_mean <- as.vector(fit.res[[11]])[1:(2*(K-1))]
  sigma_a_mean <- as.vector(fit.res[[13]])
  beta_true <- as.vector(true_val[[1]])
  alpha_true <- as.vector(true_val[[2]])
  gamma_true <- as.vector(true_val[[3]])[1:(2*(K-1))]
  sigma_a_true <- as.vector(true_val[[5]])
  
  beta_est = cbind(beta_true,beta_mean,apply(BETA_array,2,median),
                   apply(BETA_array,2,function(x) quantile(x,0.025)),
                   apply(BETA_array,2,function(x) quantile(x,0.975)),
                   as.integer(beta_true>=apply(BETA_array,2,function(x) quantile(x,0.025))&
                                beta_true<=apply(BETA_array,2,function(x) quantile(x,0.975))))
  alpha_est = cbind(alpha_true,alpha_mean,apply(ALPHA_array,2,median),
                    apply(ALPHA_array,2,function(x) quantile(x,0.025)),
                    apply(ALPHA_array,2,function(x) quantile(x,0.975)),
                    as.integer(alpha_true>=apply(ALPHA_array,2,function(x) quantile(x,0.025))&
                                 alpha_true<=apply(ALPHA_array,2,function(x) quantile(x,0.975))))
  gamma_est = cbind(gamma_true,gamma_mean,apply(GAMMA_array,2,median),
                    apply(GAMMA_array,2,function(x) quantile(x,0.025)),
                    apply(GAMMA_array,2,function(x) quantile(x,0.975)),
                    as.integer(gamma_true>=apply(GAMMA_array,2,function(x) quantile(x,0.025))&
                                 gamma_true<=apply(GAMMA_array,2,function(x) quantile(x,0.975))))
  sigma_a_est = cbind(sigma_a_true,sigma_a_mean,apply(SIGMA_A_array,2,median),
                      apply(SIGMA_A_array,2,function(x) quantile(x,0.025)),
                      apply(SIGMA_A_array,2,function(x) quantile(x,0.975)),
                      as.integer(sigma_a_true>=apply(SIGMA_A_array,2,function(x) quantile(x,0.025))&
                                   sigma_a_true<=apply(SIGMA_A_array,2,function(x) quantile(x,0.975))))
  colnames(beta_est) <- colnames(alpha_est) <- colnames(gamma_est) <- colnames(sigma_a_est) <- c("true","mean","median","lower_ci","upper_ci","coverage")
  
  
  return(list(beta_est = beta_est,
              alpha_est = alpha_est,
              gamma_est = gamma_est,
              sigma_a_est = sigma_a_est,
              K_est = K_est))
}

# Compute log-likelihood from MGMM fitting
loglik_mgmm <- function(pars_list,Y_sparse,time_sparse,M_i_sparse,X,Z,P,n_class,N)
{
  beta = pars_list[[1]]
  alpha = pars_list[[2]]
  gamma = pars_list[[3]]
  A = pars_list[[4]]
  Sigma_a = pars_list[[5]]
  s_sq_eps = pars_list[[6]]
  K_pos = round(pars_list[[7]])
  logl = matrix(0, nrow = N, ncol = n_class)
  C = as.matrix(dummy_cols(factor(K_pos, levels = 1:n_class))[,-1])
  
  for(i in 1:N)
  {
    pi_i = Z[[i]]%*%gamma
    pi_i = exp(pi_i)/apply(exp(pi_i),1,sum)
    for(k in 1:n_class)
    {
      T_i = as.matrix(bdiag(lapply(time_sparse[[i]], function(x) cbind(1,x))))
      l_ai = -1/2 * t(A[i,]-alpha[,k]) %*% solve(Sigma_a[,,k]) %*% (A[i,]-alpha[,k]) + log(det(Sigma_a[,,k])^(-1/2))
      l_Yi = rep(0,P)
      for (p in 1:P){
        Ti_p = lapply(time_sparse[[i]], function(x) cbind(1,x))[[p]]
        Yi_t = cbind(unlist(Y_sparse[[i]][[p]])) - c(X[[i]][[p]]%*%beta[,p] - Ti_p%*%(cbind(A[i,(2*p-1):(2*p)])))
        l_Yi[p] = -1/2*t(Yi_t)%*%Yi_t/s_sq_eps + log((1/s_sq_eps)^(M_i_sparse[[i]][p]/2))
      }
      logl[i,k] = l_ai + sum(l_Yi)
    }
  }
  logl_all = sum(logl * C)
  return(list(logl = logl,
              logl_all = logl_all))
}

# Compute WAIC from MGMM fitting
waic <- function(gmm.fit,sim.dat, Nsim, Niter, var = "unequal",n_class=3)
{
  lppd <- rep(0,Nsim)
  pwaic <- rep(0,Nsim)
  for(nsim in 1:Nsim)
  {
    if(nsim %% 10 ==0)
      print(nsim)
    ## load fitted parameters
    if(!is.null(gmm.fit[[nsim]]))
    {
      BETA_array=gmm.fit[[nsim]][[1]]
      ALPHA_array=do.call(rbind,lapply(gmm.fit[[nsim]][[2]],as.vector))
      A_array=gmm.fit[[nsim]][[4]]
      SIGMA_A_array=do.call(rbind,lapply(gmm.fit[[nsim]][[5]],as.vector))
      S_SQ_EPS_array=do.call(rbind,lapply(gmm.fit[[nsim]][[6]],as.vector))
      K_array=do.call(rbind,lapply(gmm.fit[[nsim]][[7]],as.vector))
      
      ## load simulated data
      Y_sparse = sim.dat[[nsim]]$Y_sparse
      time_sparse = sim.dat[[nsim]]$time_sparse
      M_i_sparse = sim.dat[[nsim]]$M_i_sparse
      X = sim.dat[[nsim]]$X
      Z = sim.dat[[nsim]]$Z
      pars_true = sim.dat[[nsim]]$pars_true
      for(i in 1:N)
      {
        y_i = cbind(unlist(Y_sparse[[i]]))
        if(length(y_i)!=0)
        {
          x_i = X[[i]]
          t_i = as.matrix(bdiag(lapply(time_sparse[[i]], function(x) cbind(1,x))))
          k_i = K_array[,i]
          p_yi <- rep(0,Niter)
          for(niter in 1:Niter)
          {
            k = k_i[niter]
            alpha_k = matrix(ALPHA_array[niter,],nrow = 2*P,ncol = n_class)[,k]
            if(var == "equal")
              sigma_k = matrix(SIGMA_A_array[niter,],2*P,2*P)
            else
              sigma_k = array(SIGMA_A_array[niter,],dim = c(2*P,2*P,n_class))[,,k]
            l_ai = mvtnorm::dmvnorm(A_array[[niter]][i,],alpha_k,sigma_k)
            zeta_ik = c(x_i%*%BETA_array[[niter]]) + t_i%*%(cbind(A_array[[niter]][i,]))
            l_y = dnorm(y_i,zeta_ik,S_SQ_EPS_array[niter])
            p_yi[niter] = l_ai * prod(l_y)
          }
          lppd[nsim] <- lppd[nsim] + log(mean(p_yi, na.rm = TRUE))
          pwaic[nsim] <- pwaic[nsim] + var(log(p_yi), na.rm = TRUE)
        }
      }
    }
    else
    {
      lppd[nsim] = NA
      pwaic[nsim] = NA
    }
  }
  WAIC <- -2*(lppd - pwaic)
  return(WAIC)
}