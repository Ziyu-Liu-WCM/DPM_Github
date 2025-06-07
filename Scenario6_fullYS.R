#==================================================================
# Script: DPMCT analysis for Scenario 6 only (single replicate)
#==================================================================
setwd("~/Library/CloudStorage/OneDrive-WeillCornellMedicine/ZiyuLiuDPM/DPM_Github")
library(ggplot2)
library(reshape2)
source("DPMCT/DPMCT_MCMC.R")
source("DPMCT/helper_functions.R")
source("DPMCT/likelihoodEta_hyp.R")
source("DPMCT/likelihoodEta_inc.R")
source("DPMCT/likelihoodEta_num.R")
source("DPMCT/updateEtaYS.R")
source("DPMCT/fn_DPMCT.R")

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

# scenario6_4x2 <- matrix(
#   c(0.5, 0.25,
#     0.3, 0.15,
#     0.6, 0.3,
#     0.4, 0.2),
#   ncol = 2,
#   byrow = TRUE
# )

#— 3) Expand the 4×2 block into 40×2 (each of the 4 rows repeated 10 times)
nBask <- 20
scenario6_40 <- do.call(
  rbind,
  lapply(seq_len(4), function(i) {
    matrix(scenario6_4x2[i, ], nrow = nBask/4, ncol = 2, byrow = TRUE)
  })
)

#— 4) Simulation parameters
n1    <- 500    # treatment sample size per basket
n0    <- 500    # control sample size per basket

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
#likelihood_version <- "inc"  # use "inc" for incremental likelihood
alphaUpdate        <- TRUE
q_threshold        <- 0.95  # use 0.95 as in previous calibration

#— 8) Run DPMCT_MCMC on this single replicate
n1_vec <- dfRep$n1
n0_vec <- dfRep$n0
y1_vec <- dfRep$y1
y0_vec <- dfRep$y0
J      <- length(n1_vec)  # should be 40

dpm_result <- DPMCT(
  n0               = n0_vec,
  y0               = y0_vec,
  n1               = n1_vec,
  y1               = y1_vec,
  burnin           = burnin,
  thin             = 5,
  savedIter       = niter / 5,  # since thin = 5
  eta              = rep(1, J),  # initial eta
  mu               = 0,
  sigma           = 2,
  alpha            = alpha,
  a_alpha = 10,
  b_alpha = 1,
  m=2,
  a_pi = 1,
  b_pi = 1,
  step_size = 0.01,
  likelihood_version = likelihood_version,
  alphaUpdate      = alphaUpdate
)

plot(dpm_result$alphaSaved, type = "l", main = "Alpha Trace", xlab = "Iteration", ylab = "Alpha")

#— 9) Evaluate performance
perf_metrics <- evaluatePerformance(
  true_ratio = dfRep$expEta,
  eta_draws  = dpm_result$eta,
  q          = q_threshold
)
perf_metrics$scenario  <- 6
perf_metrics$replicate <- 1

#— 10) Compute clustering‐probability matrix
eta_mat   <- t(dpm_result$eta) # J × nSamples)
nSamples <- ncol(eta_mat)

clust_prob <- matrix(0, nrow = J, ncol = J)
for (i in seq_len(J)) {
  for (j in seq_len(J)) {
    clust_prob[i, j] <- mean(eta_mat[i, ] == eta_mat[j, ])
  }
}

# plot heatmap about clust_prob

clust_prob_melted <- melt(clust_prob)
scenario6<-ggplot(clust_prob_melted, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(title = "Clustering Probability Matrix",
       x = "Basket Index",
       y = "Basket Index") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(scenario6, filename = "clust_prob_scenario6.png", width = 8, height = 6)
# perf_metrics$clust_prob <- clust_prob
# 
# #— 11) Save results to .RData
# out_file <- paste0(
#   "Result/DPMCT_scenario6_rep1_",
#   "niter", niter, "_",
#   "burnin", burnin, "_",
#   "alpha", alpha, "_",
#   likelihood_version, "_",
#   "alUp", alphaUpdate, "_",
#   "q", q_threshold, ".RData"
# )
# save(perf_metrics, file = out_file)
# cat("Saved results to:", out_file, "\n")
# 
# #— 12) (Optional) Quick inspections
# print(str(perf_metrics))
# image(clust_prob, main = "Clustering‐Probability Matrix for Scenario 6, Replicate 1")