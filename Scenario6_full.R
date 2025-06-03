#==================================================================
# Script: DPMCT analysis for Scenario 6 only (single replicate)
#==================================================================

source("DPMCT/DPMCT_MCMC.R")
source("DPMCT/helper_functions.R")
library(doParallel)

# — 1) Set up environment and seed
set.seed(1234)

#— 2) Define Scenario 6 (4×2 matrix: global alternative, non‐fixed η, non‐fixed control)
scenario6_4x2 <- matrix(
  c(0.5, 0.1,
    0.3, 0.15,
    0.6, 0.15,
    0.6, 0.2),
  ncol = 2,
  byrow = TRUE
)

#— 3) Expand the 4×2 block into 40×2 (each of the 4 rows repeated 10 times)
nBask <- 40
scenario6_40 <- do.call(
  rbind,
  lapply(seq_len(4), function(i) {
    matrix(scenario6_4x2[i, ], nrow = 10, ncol = 2, byrow = TRUE)
  })
)

#— 4) Simulation parameters
n1    <- 200    # treatment sample size per basket
n0    <- 200    # control sample size per basket

#— 5) Simulate one replicate’s data.frame (40 baskets)
dfRep <- data.frame(
  y1     = numeric(nBask),
  y0     = numeric(nBask),
  n1     = rep(n1, nBask),
  n0     = rep(n0, nBask),
  pi1    = numeric(nBask),
  pi0    = numeric(nBask),
  expEta = numeric(nBask)
)

for (b in seq_len(nBask)) {
  pi_t      <- scenario6_40[b, 1]
  pi_c      <- scenario6_40[b, 2]
  y1_draw   <- rbinom(1, n1, pi_t)
  y0_draw   <- rbinom(1, n0, pi_c)
  dfRep$y1[b]     <- y1_draw
  dfRep$y0[b]     <- y0_draw
  dfRep$pi1[b]    <- pi_t
  dfRep$pi0[b]    <- pi_c
  dfRep$expEta[b] <- pi_t / pi_c
}


#— 7) MCMC hyperparameters
niter              <- 200
burnin             <- 1000
alpha              <- 2
likelihood_version <- "hyp"
alphaUpdate        <- FALSE
q_threshold        <- 0.95  # use 0.95 as in previous calibration

#— 8) Run DPMCT_MCMC on this single replicate
n1_vec <- dfRep$n1
n0_vec <- dfRep$n0
y1_vec <- dfRep$y1
y0_vec <- dfRep$y0
J      <- length(n1_vec)  # should be 40

dpm_result <- DPMCT_MCMC(
  n1               = n1_vec,
  n0               = n0_vec,
  y1               = y1_vec,
  y0               = y0_vec,
  niter            = niter,
  burnin           = burnin,
  thin             = 5,
  alpha            = alpha,
  return_MCMC_spls = TRUE,
  likelihood_version = likelihood_version,
  alphaUpdate      = alphaUpdate
)

#— 9) Evaluate performance
perf_metrics <- evaluatePerformance(
  true_ratio = dfRep$expEta,
  eta_draws  = dpm_result$eta_postSpls,
  q          = q_threshold
)
perf_metrics$scenario  <- 6
perf_metrics$replicate <- 1

#— 10) Compute clustering‐probability matrix
mu_mat   <- dpm_result$mu_postSpls  # J × nSamples
nSamples <- ncol(mu_mat)

clust_prob <- matrix(0, nrow = J, ncol = J)
for (i in seq_len(J)) {
  for (j in seq_len(J)) {
    clust_prob[i, j] <- mean(mu_mat[i, ] == mu_mat[j, ])
  }
}
perf_metrics$clust_prob <- clust_prob

#— 11) Save results to .RData
out_file <- paste0(
  "Result/DPMCT_scenario6_rep1_",
  "niter", niter, "_",
  "burnin", burnin, "_",
  "alpha", alpha, "_",
  likelihood_version, "_",
  "alUp", alphaUpdate, "_",
  "q", q_threshold, ".RData"
)
save(perf_metrics, file = out_file)
cat("Saved results to:", out_file, "\n")

#— 12) (Optional) Quick inspections
print(str(perf_metrics))
image(clust_prob, main = "Clustering‐Probability Matrix for Scenario 6, Replicate 1")