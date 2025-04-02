#-------------------------------------
# Numerical Integration Version
#-------------------------------------

likelihoodEta_num <- function(n0, y0, n1, y1, a_pi, b_pi, eta) {
  
  # The integrand for the marginal likelihood integrand:
  #    [p^(y0 + y1 + a_pi - 1)] * [(1 - p)^(n0 - y0 + b_pi - 1)] * [(1 - p * exp(eta))^(n1 - y1)]
  #
  # We multiply the result of integrate(...) by the outside constant:
  #    choose(n0, y0) * choose(n1, y1) * exp(eta * y1) / Beta(a_pi, b_pi).
  
  integrand <- function(p) {
    p^(y0 + y1 + a_pi - 1) *
      (1 - p)^(n0 - y0 + b_pi - 1) *
      (1 - p * exp(eta))^(n1 - y1)
  }
  
  # Constant prefactor outside the integral
  const_factor <- choose(n0, y0) * choose(n1, y1) * exp(eta * y1) / beta(a_pi, b_pi)
  
  # Integration limits depend on whether eta <= 0 or eta > 0
  if (eta <= 0) {
    lower_limit <- 0
    upper_limit <- 1
  } else {
    lower_limit <- 0
    upper_limit <- exp(-eta)
  }
  
  # Perform the numerical integration
  val <- integrate(integrand, lower = lower_limit, upper = upper_limit,
                   rel.tol = 1e-14, # you can adjust tolerances as needed
                   abs.tol = 1e-14)
  
  # Multiply by the constant factor
  const_factor * val$value
}

# n1 <- 20
# y1 <- 5
# y0 <- 5
# n0 <- 20
# a_pi <- 1
# b_pi <- 1
# # eta <--2.1
# # check_marginal_binomial(n0, y0, n1, y1, a_pi, b_pi, eta)
# 
# eta_values <- seq(-5, 5, by = 0.01)
# results <- sapply(eta_values, function(eta) likelihoodEta_num(n0, y0, n1, y1, a_pi, b_pi, eta))
# # results <- sapply(eta_values, function(eta) likelihoodEta_inc(n0, y0, n1, y1, a_pi, b_pi, eta))
# # results <- NA
# # for (k in seq_len(length(eta_values))){
# #   results[k] <- likelihoodEta(n0, y0, n1, y1, a_pi, b_pi, eta_values[k])
# # }
# # eta_values[72]
# # likelihoodEta(n0, y0, n1, y1, a_pi, b_pi, eta)
# # eta <- eta_values[72]
# #
# #
# # Plot the likelihood function
# plot(eta_values, results, type = "l", col = "blue", lwd = 2,
#      xlab = expression(eta), ylab = "Likelihood",
#      main = "Likelihood Function with treatment(2/100), control(50/100)",
#      xaxt = "n")
# axis(1, at = seq(min(eta_values), max(eta_values), by = 0.1))
# abline(v = log(y1/y0), col = "red", lty = 2, lwd = 2)


# 
# 
# 
# 
# # Define a sequence of x values over [0, 1]
# x <- seq(0, 1, length.out = 200)
# 
# # Choose alpha and beta parameters
# alpha <- 0.5
# beta  <- 0.1
# 
# # Calculate the Beta density
# y <- dbeta(x, shape1 = alpha, shape2 = beta)
# 
# # Plot the Beta distribution
# plot(
#   x, y,
#   type = "l",       # Draw lines
#   main = "Beta Distribution",
#   xlab = "x",
#   ylab = "Density",
#   lwd = 2           # Line width
# )
# #
# # # Optionally, plot multiple Beta distributions with different parameters
# # # to see how the shape changes:
# # lines(x, dbeta(x, 2, 2), lty = 2)  # Beta(2,2)
# # lines(x, dbeta(x, 5, 1), lty = 3)  # Beta(5,1)
# 
# # Add a legend
# legend(
#   "topright",
#   legend = c("Alpha=2, Beta=5", "Alpha=2, Beta=2", "Alpha=5, Beta=1"),
#   lty = 1:3
# )