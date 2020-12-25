\name{mgmm}
\alias{mgmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mgmm(Y_sparse, time_sparse, M_i_sparse, X, Z, P = 1, n_class = 2, nloop = 110, burnin = 10, thin = 1, sim = TRUE, pars = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y_sparse}{
%%     ~~Describe \code{Y_sparse} here~~
}
  \item{time_sparse}{
%%     ~~Describe \code{time_sparse} here~~
}
  \item{M_i_sparse}{
%%     ~~Describe \code{M_i_sparse} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{Z}{
%%     ~~Describe \code{Z} here~~
}
  \item{P}{
%%     ~~Describe \code{P} here~~
}
  \item{n_class}{
%%     ~~Describe \code{n_class} here~~
}
  \item{nloop}{
%%     ~~Describe \code{nloop} here~~
}
  \item{burnin}{
%%     ~~Describe \code{burnin} here~~
}
  \item{thin}{
%%     ~~Describe \code{thin} here~~
}
  \item{sim}{
%%     ~~Describe \code{sim} here~~
}
  \item{pars}{
%%     ~~Describe \code{pars} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Y_sparse, time_sparse, M_i_sparse, X, Z, P = 1, n_class = 2, 
    nloop = 110, burnin = 10, thin = 1, sim = TRUE, pars = NULL) 
{
    library(splines)
    library(MASS)
    library(Matrix)
    library(mvtnorm)
    library(MCMCpack)
    library(bayesm)
    library(BayesLogit)
    library(fastDummies)
    N = length(Y_sparse)
    P = length(Y_sparse[[1]])
    id = 1:N
    Q_X = dim(X[[1]])[2]
    Q_Z = dim(Z[[1]])[2]
    XX = list()
    XX = t(X[[1]]) \%*\% X[[1]]
    ZZ = t(Z[[1]]) \%*\% Z[[1]]
    for (i in 2:N) {
        XX = XX + t(X[[i]]) \%*\% X[[i]]
        ZZ = ZZ + t(Z[[i]]) \%*\% Z[[i]]
    }
    M = sum(unlist(M_i_sparse))
    if (sim == FALSE) {
        est_s_sq_eps = 1
        est_beta = 1
        est_A = 1
        est_Sigma_a = 1
        est_alpha = 1
        est_K = 1
        est_gamma = 1
        s_sq_eps = 0.1
        A = array(0, dim = c(N, 2 * P))
        Sigma_a = array(1, dim = c(2 * P, 2 * P, n_class))
        for (k in 1:n_class) {
            Sigma_a[, , k] = diag(rep(1, 2 * P))
        }
        alpha = array(0, dim = c(2 * P, n_class))
        beta = array(0, dim = c(Q_X, P))
        gamma = array(0, dim = c(Q_Z, n_class))
        K = sample(1:n_class, N, replace = TRUE)
    }
    if (sim == TRUE) {
        est_beta = 1
        if (est_beta == 0) 
            beta = pars[[1]]
        if (est_beta == 1) 
            beta = array(0, dim = c(Q_X, P))
        est_A = 1
        if (est_A == 0) 
            A = pars[[4]]
        if (est_A == 1) 
            A = array(0, dim = c(N, 2 * P))
        est_alpha = 1
        if (est_alpha == 0) 
            alpha = pars[[2]]
        if (est_alpha == 1) 
            alpha = array(0, dim = c(2 * P, n_class))
        est_Sigma_a = 1
        if (est_Sigma_a == 0) 
            Sigma_a = pars[[5]]
        if (est_Sigma_a == 1) {
            Sigma_a = array(1, dim = c(2 * P, 2 * P, n_class))
            for (k in 1:n_class) {
                Sigma_a[, , k] = diag(rep(1, 2 * P))
            }
        }
        est_s_sq_eps = 1
        if (est_s_sq_eps == 0) 
            s_sq_eps = pars[[6]]
        if (est_s_sq_eps == 1) 
            s_sq_eps = 0.1
        est_K = 1
        if (est_K == 0) 
            K = pars[[7]]
        if (est_K == 1) 
            K = sample(1:n_class, N, replace = TRUE)
        est_gamma = 1
        if (est_gamma == 0) 
            gamma = pars[[3]]
        if (est_gamma == 1) 
            gamma = array(0, dim = c(Q_Z, n_class))
    }
    cc = 10000
    eta = 1
    Psi = 1 * diag(rep(1, 2 * P))
    xi = matrix(1, 1, 2 * P)
    eta1 = 1
    eta2 = 1
    nu = 100
    array_ind = 0
    A_mean = A
    beta_mean = beta
    alpha_mean = alpha
    gamma_mean = gamma
    Sigma_a_mean = Sigma_a
    s_sq_eps_mean = s_sq_eps
    K_mean = K
    Z_all = do.call(rbind, Z)
    pi_pos_mean = exp(Z_all \%*\% gamma)/apply(exp(Z_all \%*\% gamma), 
        1, sum)
    A_array = NULL
    BETA_array = NULL
    ALPHA_array = NULL
    GAMMA_array = NULL
    S_SQ_EPS_array = NULL
    SIGMA_A_array = NULL
    K_pos_array = NULL
    pi_pos_array = NULL
    K_array = matrix(0, nrow = N, ncol = n_class)
    for (i in 1:N) {
        K_array[i, K[i]] = 1
    }
    lambda = array(rep(diag(rep(1, N)), n_class - 1), dim = c(N, 
        N, n_class - 1))
    C = matrix(rep(1, N * n_class), ncol = n_class, nrow = N)
    Z_array = matrix(0, nrow = N, ncol = n_class - 1)
    for (k in 1:(n_class - 1)) {
        z_t = abs(rlogis(N, 0, 1))
        for (i in 1:N) {
            if (K_array[i, k] == 1) 
                Z_array[i, k] = z_t[i]
            else Z_array[i, k] = -z_t[i]
        }
    }
    Z_all = do.call(rbind, Z)
    gamma_array = array(0, dim = c(Q_Z, n_class, nloop))
    omega = matrix(rpg(N * n_class, 1, 0), ncol = n_class, nrow = N)
    for (ii in 1:nloop) {
        if (est_A == 1) {
            for (i in 1:N) {
                Sigma_a_k = Sigma_a[, , K[i]]
                T_i = as.matrix(bdiag(lapply(time_sparse[[i]], 
                  function(x) cbind(1, x))))
                Y_i = cbind(unlist(Y_sparse[[i]])) - c(X[[i]] \%*\% 
                  beta)
                Sigma_a_i = solve(t(T_i) \%*\% T_i/s_sq_eps + solve(Sigma_a_k))
                mu_a_i = Sigma_a_i \%*\% (t(T_i) \%*\% Y_i/s_sq_eps + 
                  solve(Sigma_a_k) \%*\% cbind(alpha[, K[i]]))
                A[i, ] = rmvnorm(1, mu_a_i, Sigma_a_i)
            }
        }
        if (est_alpha == 1) {
            for (k in 1:n_class) {
                N_k = sum(K == k)
                Sigma_alpha_k = solve(N_k * solve(Sigma_a[, , 
                  k]) + diag(rep(1, 2 * P))/cc)
                mu_alpha_k = Sigma_alpha_k \%*\% solve(Sigma_a[, 
                  , k]) \%*\% cbind(apply(matrix(A[K == k, ], ncol = 2 * 
                  P), 2, sum))
                alpha[, k] = t(rmvnorm(1, mu_alpha_k, Sigma_alpha_k))
            }
            alpha <- alpha[, order(alpha[1, ])]
        }
        if (est_Sigma_a == 1) {
            df = 2 * P
            for (k in 1:n_class) {
                lambdanew = diag(df)
                invSigma = chol2inv(chol(Sigma_a[, , k]))
                for (j in 1:df) {
                  lambdanew[j, j] = 1/rinvgamma(1, shape = ((eta + 
                    df)/2), scale = ((1/(xi[, j]^2)) + (eta * 
                    invSigma[j, j])))
                }
                eta_a = eta + sum(K == k) + df - 1
                Psi_t = sweep(matrix(A[K == k, ], ncol = 2 * 
                  P), 2, alpha[, k])
                Psi_a = 2 * eta * lambdanew + t(Psi_t) \%*\% (Psi_t)
                Sigma_a[, , k] = riwish(eta_a, Psi_a)
            }
        }
        if (est_beta == 1) {
            Sigma_beta = solve(XX/s_sq_eps + 1/cc * diag(1, dim(beta)[1]))
            for (p in 1:P) {
                mu_beta = cbind(rep(0, dim(beta)[1]))
                for (i in 1:N) {
                  T_i = cbind(1, time_sparse[[i]][[p]])
                  Y_i = cbind(Y_sparse[[i]][[p]]) - T_i \%*\% (cbind(A[i, 
                    (2 * p - 1):(2 * p)]))
                  mu_beta = mu_beta + t(X[[i]]) \%*\% Y_i/s_sq_eps
                }
                mu_beta = Sigma_beta \%*\% mu_beta
                beta[, p] = t(rmvnorm(1, mu_beta, Sigma_beta))
            }
        }
        if (est_s_sq_eps == 1) {
            eta1_eps = 1 + M/2
            eta2_eps = 1
            for (i in 1:N) {
                T_i = as.matrix(bdiag(lapply(time_sparse[[i]], 
                  function(x) cbind(1, x))))
                Y_i = cbind(unlist(Y_sparse[[i]])) - cbind(rep((X[[i]][1, 
                  ]) \%*\% beta, M_i_sparse[[i]])) - T_i \%*\% (cbind(A[i, 
                  ]))
                eta2_eps = eta2_eps + t(Y_i) \%*\% Y_i/2
            }
            s_sq_eps = rinvgamma(1, eta1_eps, eta2_eps)
        }
        if (est_gamma == 1) {
            kappa = as.matrix(dummy_cols(factor(K, levels = 1:n_class))[, 
                -1]) - 0.5
            for (k in 1:(n_class - 1)) {
                for (i in 1:N) {
                  if (n_class >= 3) 
                    C[i, k] = sum(exp(Z[[i]] \%*\% gamma[, -k]))
                  else C[i, k] = 1
                  omega[i, k] = rpg(1, 1, Z[[i]] \%*\% gamma[, 
                    k] - log(C[i, k]))
                }
                V = solve(t(Z_all) \%*\% diag(omega[, k]) \%*\% Z_all + 
                  diag(1/nu, Q_Z))
                m = V \%*\% (t(Z_all) \%*\% (kappa[, k] + diag(omega[, 
                  k]) \%*\% log(C[, k])))
                gamma_array[, k, ii] = rmvnorm(1, m, V)
            }
            gamma = gamma_array[, , ii]
        }
        pi_pos = matrix(NA, nrow = N, ncol = n_class)
        if (est_K == 1) {
            for (i in 1:N) {
                pi_i = Z[[i]] \%*\% gamma
                T_i = cbind(1, time_sparse[[i]][[1]])
                for (k in 1:n_class) {
                  pi_i[k] = pi_i[k] - 0.5 * t(A[i, ] - alpha[, 
                    k]) \%*\% solve(Sigma_a[, , k]) \%*\% (A[i, ] - 
                    alpha[, k]) - 0.5 * log(det(Sigma_a[, , k]))
                }
                pi_i = pi_i - max(pi_i)
                pi_i = exp(pi_i)/apply(exp(pi_i), 1, sum)
                K[i] = which(rmultinom(1, 1:n_class, pi_i) == 
                  1)
                pi_pos[i, ] = pi_i
            }
        }
        if (ii\%\%thin == 0 & ii >= burnin) {
            array_ind = array_ind + 1
            ALPHA_array[[array_ind]] = alpha
            BETA_array[[array_ind]] = beta
            GAMMA_array[[array_ind]] = gamma
            A_array[[array_ind]] = A
            SIGMA_A_array[[array_ind]] = Sigma_a
            S_SQ_EPS_array[[array_ind]] = s_sq_eps
            K_pos_array[[array_ind]] = K
            pi_pos_array[[array_ind]] = pi_pos
            alpha_mean = ((array_ind - 1) * alpha_mean + alpha)/array_ind
            beta_mean = ((array_ind - 1) * beta_mean + beta)/array_ind
            gamma_mean = ((array_ind - 1) * gamma_mean + gamma)/array_ind
            A_mean = ((array_ind - 1) * A_mean + A)/array_ind
            Sigma_a_mean = ((array_ind - 1) * Sigma_a_mean + 
                Sigma_a)/array_ind
            s_sq_eps_mean = ((array_ind - 1) * s_sq_eps_mean + 
                s_sq_eps)/array_ind
            K_mean = ((array_ind - 1) * K_mean + K)/array_ind
            pi_pos_mean = ((array_ind - 1) * pi_pos_mean + pi_pos)/array_ind
        }
    }
    mcmc_results = list()
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
    mcmc_results[[16]] = pi_pos_array
    mcmc_results[[17]] = pi_pos_mean
    return(mcmc_results)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
