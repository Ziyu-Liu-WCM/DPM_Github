DPMcalc_FWER_for_q <- function(q, scenario_data, niter = 100, burnin = 1000, alpha = 2,
                            likelihood_version ="inc", alphaUpdate = FALSE){
  # scenario_data is sim_data[[1]] for GN scenario
  # We'll loop over the 1000 replicates and see how often
  # we reject at least 1 basket.
  
  nRep <- length(scenario_data)
  # Use foreach to iterate over replicates in parallel
  rejects <- foreach(i = seq_len(nRep), .combine = c, 
                     .export = c()) %dopar% {
                       # Extract replicate i
                       
                       n0 <- scenario_data[[i]]$n0
                       y1 <- scenario_data[[i]]$y1
                       n1 <- scenario_data[[i]]$n1
                       y0 <- scenario_data[[i]]$y0
                       
                       dpm <- DPMCT_MCMC(n1, n0, y1, y0,
                                         niter=niter,
                                         burnin=burnin,
                                         thin=5,
                                         alpha=alpha,
                                         return_MCMC_spls=TRUE,
                                         likelihood_version=likelihood_version,
                                         alphaUpdate=alphaUpdate)
                       PP <- apply(dpm$eta_postSpls, 1, function(x) mean(x>0))
                       
                       return(any(PP > q))
                    }
  
  # 'rejects' is now a logical vector; take mean (TRUE=1, FALSE=0)
  FWER_est <- mean(rejects)
  return(FWER_est)
}


DPMfind_q_for_5pct <- function(scenario_data, niter = 100, burnin = 1000, alpha = 2,
                               likelihood_version ="inc", alphaUpdate = FALSE){
  # We'll do a small search or use uniroot:
  # We want f(q) = calc_FWER_for_q(q) - 0.05 = 0
  # q must be between 0.75 and 1.0 realistically. We'll do uniroot:
  
  f <- function(q) DPMcalc_FWER_for_q(q, scenario_data, niter = niter, burnin = burnin, alpha = alpha,
                                      likelihood_version =likelihood_version, alphaUpdate = alphaUpdate) - 0.05

  # We can attempt a bracket search. We'll try from around 0.75 to 0.999
  # Since the final q is often ~ 0.96-0.99 for 4 baskets.
  out <- uniroot(f, interval=c(0.01, 1.0), tol=1e-3)
  return(out$root)
}



evaluatePerformance <- function(
    # True data info:
  true_ratio,       # numeric vector of length J, the true expEta for each basket
  # Posterior inference info:
  eta_draws,        # J x (# draws), posterior samples of the log ratio, for each basket
  # Decision threshold:
  q
) {
  # We assume:
  #  - A "null" basket means true_ratio[basket] == 1.
  #  - A "reject H0" means mean( eta_draws[basket,] > 0 ) > q.
  #
  # We'll compute:
  #  1) FWER
  #  2) Family-wise power (two subtypes)
  #  3) Basket-wise type I error / power
  #  4) Bias & MSE in terms of eta
  # 
  # Return these in a list.
  
  J <- length(true_ratio)          # number of baskets
  true_eta <- log(true_ratio)      # the "true log ratio" for each basket
  
  # Posterior probability that basket j > 0
  PP <- apply(eta_draws, 1, function(x) mean(x > 0))  # length J
  
  #--- 1) Family-wise Type I error rate (FWER)
  # A "false discovery" happens if we reject any basket that is truly null
  # So let's define an indicator for each basket: reject[j] = (PP[j] > q)
  # Then if a truly null basket is rejected, that replicate's indicator is 1
  reject <- (PP > q)
  
  # which baskets are truly null?
  isNull <- (true_ratio == 1)
  
  # do we have at least one false positive in any null basket?
  anyFalsePos <- any(reject[isNull])
  
  # "any false discovery" => 1 if yes, else 0
  # Because this function is for a single replicate, that is either 0 or 1
  FWER <- as.numeric(anyFalsePos)
  
  #--- 2) Family-wise power
  # Subtype A: "any true alternative is rejected while no true null is rejected"
  #   => check if we have zero false rejections, and at least one truly alt is rejected
  # Subtype B: "all true alternatives are rejected while no true null is rejected"
  
  isAlt <- !isNull  # i.e. ratio != 1
  noFalseRejects   <- !any(reject[isNull])
  anyAltRejected   <- any(reject[isAlt])
  allAltRejected   <- all(reject[isAlt]) || sum(isAlt)==0  # if no alt baskets, trivial
  
  # Subtype A:
  famPowerA <- as.numeric(noFalseRejects && anyAltRejected)
  # Subtype B:
  famPowerB <- as.numeric(noFalseRejects && allAltRejected)
  
  #--- 3) Basket-wise type I error / power
  # If basket j is truly null, then type I error indicator is 1 if we rejected that basket.
  # If basket j is truly alt, then "power" indicator is 1 if we rejected that basket.
  basketTypeI <- rep(0, J)
  basketPow   <- rep(0, J)
  for(j in seq_len(J)) {
    if(isNull[j]) {
      basketTypeI[j] <- as.numeric(reject[j])  # 1 if we wrongly reject
      basketPow[j]   <- NA                    # not defined for a null basket
    } else {
      basketTypeI[j] <- NA
      basketPow[j]   <- as.numeric(reject[j])
    }
  }
  
  #--- 4) Bias & MSE in terms of eta
  # For basket j, the posterior mean of eta_j is mean(eta_draws[j,]), call that postMean_j
  # Then the bias for basket j is E[ postMean_j - true_eta_j ] over replicate draws.
  # On a per-replicate basis, we just compute (postMean_j - true_eta_j).
  # Similarly, MSE_j = E[(postMean_j - true_eta_j)^2].
  postMean_eta <- rowMeans(eta_draws)  # length J
  diff         <- postMean_eta - true_eta
  bias_j       <- diff     # length J
  mse_j        <- diff^2   # length J
  
  # Return everything in a list
  list(
    FWER         = FWER,
    famPowerA    = famPowerA,
    famPowerB    = famPowerB,
    basketTypeI  = basketTypeI,
    basketPower  = basketPow,
    bias         = bias_j,
    mse          = mse_j
  )
}