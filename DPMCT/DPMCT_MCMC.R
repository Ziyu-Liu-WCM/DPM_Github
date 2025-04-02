
source("DPMCT/likelihoodEta_hyp.R")
source("DPMCT/likelihoodEta_inc.R")
source("DPMCT/likelihoodEta_num.R")
source("DPMCT/fn_DPMCT.R")



DPMCT_MCMC = function(n1, n0, y1, y0,  niter, burnin = 5000, thin = 5, 
                      alpha = 2, return_MCMC_spls = TRUE,
                      likelihood_version = "inc",
                      alphaUpdate = FALSE){
  
  ###########################################################################
  # INPUT: 
  # n1: J length vector, total number of patients in each treatment group
  # n0: J length vector, total number of patients in each control group
  # y1: J length vector, total number of response in each treatment group
  # y0: J length vector, total number of response in each control group
  # niter: number of posterior samples we want 
  # burnin: number of burn-in iterations
  # thin: Length of thinning interval for taking the final samples
  # alpha: The alpha parameter for DPM
  # alphaUpdate: Whether update alpha in MCMC
  ###########################################################################
  
  J <- length(n0)
  
  if(all(sapply(list(n1, n0, y1, y0), length) != J)){
    stop("Length of n1,no,y1,y0 must be equal.")
  }
  
  if(max(n1) >30 | max(n0) > 30){
    message("Forcing likelihood_version to be hyp")
    likelihood_version = "hyp"
  }
  
  #################################################################
  # Set hyperparameters
  #################################################################
  
  a_pi <- 0.1
  b_pi <- 0.1
  step_size <- 1
  mu_prior_mean <- 0
  mu_prior_sd <- 2
  a_tau <- 1
  b_tau <- 1
  a_alpha <- 1
  b_alpha <- 1
  m <- 5

  
  
  #################################################################
  # Initialize Markov chains
  #################################################################
  eta_init <- rep(1, J)
  mu_init <- rep(0, J)
  tau_init <- rep(1, J)
  alpha_init <- 2
  
  #################################################################
  # MCMC burnin
  #################################################################
  
  output_burnin <- DPMCT(n0 = n0, y0 = y0, 
                  n1 = n1, y1 = y1,
                  niter = burnin, eta = eta_init,
                  mu = mu_init, tau = tau_init, 
                  alpha = alpha_init, 
                  mu_prior_mean = mu_prior_mean, mu_prior_sd = mu_prior_sd, 
                  a_tau = a_tau, b_tau = b_tau,
                  a_alpha = a_alpha, b_alpha = b_alpha, 
                  m = m, 
                  a_pi = a_pi, b_pi = b_pi, 
                  step_size = step_size, 
                  likelihood_version = likelihood_version, 
                  alphaUpdate = alphaUpdate)
  
  eta_burned = output_burnin$eta_spls[burnin,]
  mu_burned = output_burnin$mu_spls[burnin,]
  tau_burned = output_burnin$tau_spls[burnin,]
  alpha_burned = output_burnin$alpha_spls[burnin]
  
  #################################################################
  # MCMC iterations
  #################################################################
  
  totalIter <- (niter - 1) * thin + 1
  output_post <- DPMCT(n0 = n0, y0 = y0,
                         n1 = n1, y1 = y1,
                         niter = totalIter, eta = eta_burned,
                         mu = mu_burned, tau = tau_burned,
                         alpha = alpha_burned,
                         mu_prior_mean = mu_prior_mean, mu_prior_sd = mu_prior_sd,
                         a_tau = a_tau, b_tau = b_tau,
                         a_alpha = a_alpha, b_alpha = b_alpha,
                         m = m,
                         a_pi = a_pi, b_pi = b_pi,
                         step_size = step_size,
                         likelihood_version = likelihood_version,
                         alphaUpdate = alphaUpdate)

  postIndex <- seq(from = 1, to = totalIter, by = thin)

  eta_postSpls = output_post$eta_spls[postIndex,]
  mu_postSpls = output_post$mu_spls[postIndex,]
  tau_postSpls = output_post$tau_spls[postIndex,]
  alpha_postSpls = output_post$alpha_spls[postIndex]
  
  if(return_MCMC_spls){
    
    MCMC_spls = list()
    MCMC_spls$eta_postSpls = t(eta_postSpls)
    MCMC_spls$mu_postSpls = t(mu_postSpls)
    MCMC_spls$tau_postSpls = t(tau_postSpls)
    MCMC_spls$alpha_postSpls = alpha_postSpls
    return(MCMC_spls)
    
  } else {
    
    return(eta_postSpls)
    
  }
  
}

# niter = 20
# burnin = 50
# spls <- DPMCT_MCMC(n1, n0, y1, y0, niter, burnin = burnin, thin = 5,
#            alpha = 2, return_MCMC_spls = TRUE,
#            likelihood_version = "inc",
#            alphaUpdate = FALSE)
# eta_post <- spls$eta_postSpls
# apply(eta_post > 0, 1, sum) / dim(eta_post)[2]


