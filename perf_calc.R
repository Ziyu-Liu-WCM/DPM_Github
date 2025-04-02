##############################
#This script is for computing performance, including heatmap, but it's in a mess. 
#Contact me if you really need to run this script sysmetically
##############################



#########
##Heatmap
#########

library(fields)
clust_prob_avg <- array(NA, dim = c(40, 40, 5))

# Loop over the 5 scenarios
for (k in seq_along(resAll)) {
  # Number of replications in scenario k
  nReps <- length(resAll[[k]])
  
  # Create an array to hold all "clust_prob" matrices for scenario k
  # Its dimensions are 4 x 4 x nReps
  clust_prob_all <- array(NA, dim = c(40, 40, nReps))
  
  # Loop over each replication within scenario k
  for (i in seq_len(nReps)) {
    # Extract the 4x4 matrix "clust_prob"
    clust_prob_all[ , , i] <- resAll[[k]][[i]]$clust_prob
  }
  
  # Take the average across the third dimension (replications)
  # This collapses the 3D array into a 2D 4x4 matrix of means
  clust_prob_avg[ , , k] <- apply(clust_prob_all, c(1, 2), mean)
}



J <- 40
image.plot(
  x         = 1:J,
  y         = 1:J,
  z         = clust_prob_avg[,,5],
  main      = "Mean Same-μ Clustering Probability (200 datasets), Updating Alpha",
  xlab      = "Basket",
  ylab      = "Basket",
  legend.lab= "P(same cluster)",
  zlim      = c(0, 1),   # probabilities between 0 and 1
  axes      = FALSE
)

axis(1, at = 1:J, labels = paste("B", 1:J, sep=""))
axis(2, at = 1:J, labels = paste("B", 1:J, sep=""))



### DPM
library(dplyr)

dfCombined <- bind_rows(
  lapply(resAll, function(lst) {
    # 'lst' is the list of replicates for this scenario
    # each element is the named list from evaluatePerformance
    # convert each to a data frame with e.g. as.data.frame or tibble
    lapply(lst, as.data.frame)
  }) %>% unlist(recursive = FALSE)
)



dfSummary <- dfCombined %>%
  group_by(scenario) %>%
  summarize(
    FWER        = mean(FWER, na.rm = TRUE),
    FWP1   = mean(famPowerA, na.rm = TRUE),
    FWP2   = mean(famPowerB, na.rm = TRUE),
    basket_T1E   = mean(basketTypeI, na.rm = TRUE),  # basket-wise type I error
    basket_power   = mean(basketPower, na.rm = TRUE),  # basket-wise power
    bias    = mean(bias, na.rm = TRUE),         # average bias(eta_hat - eta)
    MSE     = mean(mse, na.rm = TRUE)           # average MSE
  )

dfSummary



dfSummary1 <- dfSummary

dfSummary_All <- rbind(dfSummary1, dfSummary)


dfSummary_All$group <- c("FXCR Global Null", "FXCR Global Alternative", "FXCR Mixed Null/Alternative", "nonFXCR Global Null", "nonFXCR Global Alternative(Fixed η)",
                         "nonFXCR Global Alternative(non-Fixed η)", "nonFXCR Mixed Null/Alternative 1", "nonFXCR Mixed Null/Alternative 2")

dfSummary_All<-dplyr::select(dfSummary_All, c('group'), everything())










### RoBoT
# 1. Flatten each sub‐list and row‐bind
resAll_df <- do.call(
  rbind,
  lapply(resAll, function(sublist) {
    # For each variable in this sublist, if it's a vector, compute mean; if scalar, just keep it
    sapply(sublist, function(x) {
      if (length(x) > 1) mean(x) else x
    })
  })
)

resAll_df <- as.data.frame(resAll_df)


resAll_df1 <- resAll_df


resAll_dfAll <- rbind(resAll_df1, resAll_df)


resAll_dfAll$group <- c("FXCR Global Null", "FXCR Global Alternative", "FXCR Mixed Null/Alternative", "nonFXCR Global Null", "nonFXCR Global Alternative(Fixed η)",
                         "nonFXCR Global Alternative(non-Fixed η)", "nonFXCR Mixed Null/Alternative 1", "nonFXCR Mixed Null/Alternative 2")

resAll_dfAll<-dplyr::select(resAll_dfAll, c('group'), everything())


