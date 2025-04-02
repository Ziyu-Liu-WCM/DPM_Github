source("DPMCT/DPMCT_MCMC.R")
source("generate_data_largebaskSamp.R")
source("DPMCT/helper_functions.R")
library(doParallel)



niter <- 200
burnin <- 1000
alpha <- 2
likelihood_version <- 'hyp'
alphaUpdate <- FALSE
nCores <- 20

if(min(nonFXCR_simData[[1]][[1]]$n1) > 30) largeSample = TRUE else largeSample = FALSE
####################################################################
# Non-Fixed Control Rates
###################################################################
set.seed(8888)


start_time <- Sys.time()
####################################################################
# 1) Calibrate threshold q under Global Null For Fixed Control Rate
###################################################################
cl <- makeCluster(nCores)
# Register it for parallel use:
registerDoParallel(cl)

clusterEvalQ(cl, {
  # source("fn_RoBoT.R")  # This file defines run_RoBoT_once, RoBoT_MCMC, ...
  source("DPMCT/DPMCT_MCMC.R")
  source("DPMCT/helper_functions.R")
})

FXCR_GN <- FXCR_simData[[1]]

#####################################################################################
## We don't need to calibrate q if we just want to look at the clustering performance
#####################################################################################
# q_calibrated_FX <- DPMfind_q_for_5pct(FXCR_GN, niter = niter, burnin = burnin, alpha = alpha,
#                                    likelihood_version =likelihood_version, alphaUpdate = alphaUpdate)
q_calibrated_FX <- 0.95
cat("Calibrated q for Fixed Control Rate Scenarios:", q_calibrated_FX, "\n")


# Master results list, one element per scenario
resAll <- vector("list", length(FXCR_simData))



for(scen_idx in seq_along(FXCR_simData)){
  
  scenario_data <- FXCR_simData[[scen_idx]]
  nRep          <- length(scenario_data)
  
  resList <- foreach(i = seq_len(nRep), .export = c()) %dopar% {
                       # Extract replicate i
                       n0 <- scenario_data[[i]]$n0
                       y1 <- scenario_data[[i]]$y1
                       n1 <- scenario_data[[i]]$n1
                       y0 <- scenario_data[[i]]$y0
                       J <- length(n0)
                       
                       dpm <- DPMCT_MCMC(n1, n0, y1, y0,
                                         niter=niter,
                                         burnin=burnin,
                                         thin=5,
                                         alpha=alpha,
                                         return_MCMC_spls=TRUE,
                                         likelihood_version=likelihood_version,
                                         alphaUpdate=alphaUpdate)
                       
                       perf <- evaluatePerformance(
                         true_ratio = scenario_data[[i]]$expEta,
                         eta_draws  = dpm$eta_postSpls,
                         q          = q_calibrated_FX   # hypothetical threshold
                       )
                       
                       
                       # Optionally include scenario & replicate indices in the returned result
                       perf$scenario  <- scen_idx
                       perf$replicate <- i
                       
                       # Compute the clustering probability matrix for dataset k
                       mu_mat   <-dpm$mu_postSpls  # dimension: J x numSamples
                       nSamples <- ncol(mu_mat)
                       
                       clust_prob <- matrix(0, nrow = J, ncol = J)
                       for (row_idx in 1:J) {
                         for (col_idx in 1:J) {
                           clust_prob[row_idx, col_idx] <- mean(mu_mat[row_idx, ] == mu_mat[col_idx, ])
                         }
                       }
                       
                       # Store this J x J matrix
                       perf$clust_prob <- clust_prob
                       
                       return(perf)
                       }
  
  # Store the replicate‐wise results in resAll
  resAll[[scen_idx]] <- resList
}
stopCluster(cl)

cat("Current Working Directory:", getwd(), "\n")
outString <- paste("_niter", niter, "_burnin", burnin, "_alpha", alpha, "_", likelihood_version, "_", alphaUpdate, "_q", round(q_calibrated_FX, digits = 3), sep = "")
if(largeSample) outString <- paste(outString, "_largeSample", sep = "")
outString <- paste(outString, "_40B", sep = "")

outFileName <- paste("Result/DPMCT_FXCR", outString, ".RData", sep = "")
save(resAll, file = outFileName)

cat("Saved File at:", outFileName, "\n")

end_time <- Sys.time()
run_time <- as.numeric(end_time - start_time, units = "mins")
cat("Run Time for Fixed Control Rate Scenarios:", run_time, "Minutes\n")






####################################################################
# Non-Fixed Control Rates
###################################################################
set.seed(8888)



start_time <- Sys.time()
####################################################################
# Calibrate threshold q under Global Null For non-Fixed Control Rate
###################################################################
cl <- makeCluster(nCores)
# Register it for parallel use:
registerDoParallel(cl)

clusterEvalQ(cl, {
  # source("fn_RoBoT.R")  # This file defines run_RoBoT_once, RoBoT_MCMC, ...
  source("DPMCT/DPMCT_MCMC.R")
  source("DPMCT/helper_functions.R")
})

nonFXCR_GN <- nonFXCR_simData[[1]]

# q_calibrated_nonFX <- DPMfind_q_for_5pct(nonFXCR_GN, niter = niter, burnin = burnin, alpha = alpha,
#                                    likelihood_version =likelihood_version, alphaUpdate = alphaUpdate)

q_calibrated_nonFX <- 0.95

cat("Calibrated q for non-Fixed Control Rate Scenarios:", q_calibrated_nonFX, "\n")


# Master results list, one element per scenario
resAll <- vector("list", length(nonFXCR_simData))



for(scen_idx in seq_along(nonFXCR_simData)){
  
  scenario_data <- nonFXCR_simData[[scen_idx]]
  nRep          <- length(scenario_data)
  
  resList <- foreach(i = seq_len(nRep), .export = c()) %dopar% {
    # Extract replicate i
    n0 <- scenario_data[[i]]$n0
    y1 <- scenario_data[[i]]$y1
    n1 <- scenario_data[[i]]$n1
    y0 <- scenario_data[[i]]$y0
    J <- length(n0)
    
    dpm <- DPMCT_MCMC(n1, n0, y1, y0,
                      niter=niter,
                      burnin=burnin,
                      thin=5,
                      alpha=alpha,
                      return_MCMC_spls=TRUE,
                      likelihood_version=likelihood_version,
                      alphaUpdate=alphaUpdate)
    
    perf <- evaluatePerformance(
      true_ratio = scenario_data[[i]]$expEta,
      eta_draws  = dpm$eta_postSpls,
      q          = q_calibrated_nonFX   # hypothetical threshold
    )
    
    # Optionally include scenario & replicate indices in the returned result
    perf$scenario  <- scen_idx
    perf$replicate <- i
    
    # Compute the clustering probability matrix for dataset k
    mu_mat   <-dpm$mu_postSpls  # dimension: J x numSamples
    nSamples <- ncol(mu_mat)
    
    clust_prob <- matrix(0, nrow = J, ncol = J)
    for (row_idx in 1:J) {
      for (col_idx in 1:J) {
        clust_prob[row_idx, col_idx] <- mean(mu_mat[row_idx, ] == mu_mat[col_idx, ])
      }
    }
    
    # Store this J x J matrix
    perf$clust_prob <- clust_prob
    
    return(perf)
  }
  
  # Store the replicate‐wise results in resAll
  resAll[[scen_idx]] <- resList
}
stopCluster(cl)

cat("Current Working Directory:", getwd(), "\n")
outString <- paste("_niter", niter, "_burnin", burnin, "_alpha", alpha, "_", likelihood_version, "_", alphaUpdate, "_q", round(q_calibrated_nonFX, digits = 3), sep = "")
if(largeSample) outString <- paste(outString, "_largeSample", sep = "")
outString <- paste(outString, "_40B", sep = "")

outFileName <- paste("Result/DPMCT_nonFXCR", outString, ".RData", sep = "")
save(resAll, file = outFileName)

cat("Saved File at:", outFileName, "\n")

end_time <- Sys.time()
run_time <- as.numeric(end_time - start_time, units = "mins")
cat("Run Time for non-Fixed Control Rate Scenarios:", run_time, "Minutes\n")