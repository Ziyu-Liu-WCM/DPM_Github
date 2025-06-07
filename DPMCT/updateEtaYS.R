updateEta <- function(n0, y0, n1, y1, 
                      eta, mu, sigma, alpha, m,
                      a_pi, b_pi, step_size, likelihood_version = "inc"){
  if(likelihood_version == "hyp"){
    likelihoodEta <- likelihoodEta_hyp
  } else if(likelihood_version == "inc") {
    likelihoodEta <- likelihoodEta_inc
  } else if(likelihood_version == "num"){
    likelihoodEta <- likelihoodEta_num
  }
  J <- length(n0)
  
  for (i in 1:length(n0)) {
    old_eta  <- eta[i]
    
    # Remove i from its old cluster by "temporarily" clearing out mu[i], sigma[i].
    # We'll reassign it below.  For checking if the old cluster was a singleton,
    # see how many other j != i share the same (old_mu, old_sigma).
    eta[i]  <- NA
    
    # Check if the old cluster is now empty
    # i.e., no other j has (mu[j], sigma[j]) = (old_mu, old_sigma).

    is_singleton <- sum(eta == old_eta, na.rm=TRUE)>0
    
    # Gather unique clusters among the other data points
    eta_other<-eta[-i]
    uniqueEta<-unique(eta_other)
    K <- length(uniqueEta)
    
    #-----------------------------------
    # Build new proposals for i
    #-----------------------------------
    eta_pro  <- rnorm(m,mu,sigma)
    if (is_singleton) {
      # Recycle i's old parameter as first new proposal:
      eta_pro[1]  <- old_eta
      }
    
    #-------------------------------------
    # Compute unnormalized log probabilities
    #-------------------------------------
    # 1) for existing clusters
    log_prob <- numeric(K + m)
    for (k in seq_len(K)) {
      clusterEta  <- uniqueEta[k] 
      indi<-which(eta==clusterEta)
      nk <- sum(length(indi))
      # unnormalized log probability = log(n_k) + log likelihood
      log_prob[k] <- log(nk) + log(likelihoodEta(n0[i], y0[i], n1[i], y1[i], a_pi, b_pi, clusterEta))
    }
    # 2) for new proposals
    for (mm in seq_len(m)) {
      log_prob[K + mm] <- log(alpha)-log(m) + 
        log(likelihoodEta(n0[i], y0[i], n1[i], y1[i], a_pi, b_pi,eta_pro[mm]))
    }
    
    #-------------------------------------
    # Sample a new cluster for data i
    #-------------------------------------
    lp_max <- max(log_prob)
    w <- exp(log_prob - lp_max)
    w <- w / sum(w)
    pick <- which(cumsum(w) > runif(1))[1]
    eta[i]<-ifelse(pick<=K,uniqueEta[pick],eta_pro[pick - K])
   }
  uniqEta<-unique(eta)
  for(i in 1:length(uniqueEta)){
    etaNow<-uniqEta[i]
    etaProp<-etaNow+rnorm(1,0,step_size)
    indi<-which(eta==etaNow)
    dfTemp<-data.frame(n0,n1,y0,y1)
    obsSubset<-dfTemp[indi,]
    likelihoodCurrent<-dnorm(etaNow, mean = mu, sd = sigma, log = TRUE)
    likelihoodPro<-dnorm(etaProp, mean = mu, sd = sigma, log = TRUE)
    for(j in 1:length(indi)){
      likelihoodCurrent<-likelihoodCurrent+log(likelihoodEta(obsSubset$n0[j], obsSubset$y0[j], obsSubset$n1[j], obsSubset$y1[j], a_pi, b_pi, etaNow))
      likelihoodPro<-likelihoodPro+log(likelihoodEta(obsSubset$n0[j], obsSubset$y0[j], obsSubset$n1[j], obsSubset$y1[j], 
                                                     a_pi, b_pi, etaProp))
    }
    u <- runif(1)
    if (log(u) < (likelihoodPro - likelihoodCurrent)) {
      eta[indi] <- etaProp
    } else {
      eta[indi] <- etaNow
    }
  }
  
  # Return the updated mu and sigma
  print(eta)
  return(eta)
}