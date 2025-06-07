
#source("DPMCT/updateMusigma.R")
#source("DPMCT/updateAlpha.R")

DPMCT <- function(n0, y0, n1, y1, burnin, thin, savedIter, eta, mu, sigma, alpha, a_alpha, b_alpha, m, a_pi, b_pi, step_size, likelihood_version = "inc", alphaUpdate = FALSE){
  
  J <- length(n0)
  chain_alpha <- numeric(burnin+thin*savedIter)
  chain_eta   <- matrix(NA, nrow = savedIter, ncol = J)
  
  for (iter in 1:(burnin+thin*savedIter)) {

# show the iteration number and calculate the time for each iteration
    # Print iteration number
    cat("Iteration:", iter, "\n")
    # Start time for the iteration
     start_time <- Sys.time()

    
    # A) Update alpha
    if(alphaUpdate) alpha <- update_alpha(alpha, mu, sigma, a_alpha, b_alpha)
    # C) Update eta using binomial data
    eta <- updateEta(n0, y0, n1, y1, 
                     eta, mu, sigma, alpha, m,
                     a_pi, b_pi, step_size, likelihood_version = likelihood_version)
    cat("number of clusters",length(unique(eta)), "\n")
    # D) Store samples
    chain_alpha[iter]    <- alpha
    if (iter > burnin && (iter - burnin) %% thin == 0) {
      chain_eta[(iter - burnin) / thin, ] <- eta
      cat("Storing eta for iteration", (iter - burnin) / thin, "\n")
    }
    
  }
  
  return(list(  burnin           = burnin,
                thin             = 5,
                savedIter       = savedIter,
                alpha            = alpha,
                a_alpha = a_alpha,
                b_alpha = b_alpha,
                a_pi = a_pi,
                b_pi = b_pi,
                step_size = step_size,
                mu=mu,
                sigma=sigma,
                m=m,
                likelihood_version = likelihood_version,
                alphaUpdate      = alphaUpdate,
    etaSaved= chain_eta,alphaSaved=chain_alpha))
}

## Single iteration
# eta_updated <- update_eta(n0, y0, n1, y1, eta, a_pi, b_pi, mu, sigma)
# mu_sigma_updated <- update_mu_sigma(eta_updated, mu, sigma, alpha,
#                             mu_prior_mean, mu_prior_sd,
#                             a_sigma, b_sigma, m)
# mu_updated <- mu_sigma_updated$mu
# sigma_updated <- mu_sigma_updated$sigma
# alpha_updated <- update_alpha(alpha, mu_updated, sigma_updated, a_alpha, b_alpha)


































# eta_current <- eta
# mu_current <- mu
# sigma_current <- sigma
# alpha_current <- alpha
# num_iter <- 10000
# 
# eta_samples <- matrix(NA, nrow = num_iter, ncol = J)
# mu_samples <- matrix(NA, nrow = num_iter, ncol = J)
# sigma_samples <- matrix(NA, nrow = num_iter, ncol = J)
# alpha_samples <- numeric(num_iter)
# 
# # for (iter in 1:num_iter) {
# #   eta_current <- update_eta(n0, y0, n1, y1, eta_current, a_pi, b_pi, mu_current, sigma_current)
# #   eta_samples[iter, ] <- eta_current
# # }
# 
# for (iter in 1:num_iter) {
#   # Update eta
#   eta_current <- update_eta(n0, y0, n1, y1, eta_current, a_pi, b_pi, mu_current, sigma_current)
#   
#   # Update mu and sigma
#   mu_sigma_updated <- update_mu_sigma(eta_current, mu_current, sigma_current, alpha_current,
#                                   mu_prior_mean, mu_prior_sd, a_sigma, b_sigma, m)
#   mu_current <- mu_sigma_updated$mu
#   sigma_current <- mu_sigma_updated$sigma
#   
#   # Update alpha
#   alpha_current <- update_alpha(alpha_current, mu_current, sigma_current, a_alpha, b_alpha)
#   
#   # Store samples
#   eta_samples[iter, ] <- eta_current
#   mu_samples[iter, ] <- mu_current
#   sigma_samples[iter, ] <- sigma_current
#   alpha_samples[iter] <- alpha_current
# }
# 
# burn_in <- 9000
# 
# posterior_samples <- eta_samples[(burn_in : num_iter),1]
# 
# posterior_mean <- mean(posterior_samples)
# posterior_mode <- as.numeric(density(posterior_samples)$x[which.max(density(posterior_samples)$y)])
# 
# 
# hist(posterior_samples, breaks = 50, probability = TRUE, 
#      main = "Posterior Distribution with treatment(6/10), control(3/10)", xlab = "Parameter Value", 
#      col = "lightblue", border = "white")
# 
# # Overlay density plot
# lines(density(posterior_samples), col = "blue", lwd = 2)
# 
# # Add vertical lines for mean and mode
# abline(v = posterior_mean, col = "red", lwd = 2, lty = 2)  # Mean
# abline(v = posterior_mode, col = "green", lwd = 2, lty = 2)  # Mode
# 
# # Add legend
# legend("topright", legend = c("Mean", "Mode"), 
#        col = c("red", "green"), lty = 2, lwd = 2, bty = "n")
# 
# # Annotate the mean and mode
# text(posterior_mean, 0.05, labels = paste0("Mean = ", round(posterior_mean, 2)), col = "red", pos = 4)
# text(posterior_mode, 0.04, labels = paste0("Mode = ", round(posterior_mode, 2)), col = "green", pos = 4)

