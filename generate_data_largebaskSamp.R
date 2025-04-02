set.seed(1234)

#--- 8 scenarios, each is originally a 4 x 2 matrix (columns = (treatment_rate, control_rate))
scenarios <- list(
  # 1) Global Null, fixed control
  matrix(c(0.1,0.1,
           0.1,0.1,
           0.1,0.1,
           0.1,0.1), ncol=2, byrow=TRUE),
  # 2) Global Alternative, fixed control
  matrix(c(0.3,0.1,
           0.3,0.1,
           0.3,0.1,
           0.3,0.1), ncol=2, byrow=TRUE),
  # 3) Mixed Null/Alternative, fixed control
  matrix(c(0.1,0.1,
           0.1,0.1,
           0.5,0.1,
           0.5,0.1), ncol=2, byrow=TRUE),
  # 4) Global Null, non‐fixed control
  matrix(c(0.1,0.1,
           0.15,0.15,
           0.15,0.15,
           0.2,0.2), ncol=2, byrow=TRUE),
  # 5) Global Alternative (fixed η), non‐fixed control
  matrix(c(0.3,0.1,
           0.45,0.15,
           0.45,0.15,
           0.6,0.2), ncol=2, byrow=TRUE),
  # 6) Global Alternative (non‐fixed η), non‐fixed control
  matrix(c(0.5,0.1,
           0.3,0.15,
           0.6,0.15,
           0.6,0.2), ncol=2, byrow=TRUE),
  # 7) Mixed Null/Alternative 1, non‐fixed control
  matrix(c(0.3,0.1,
           0.45,0.15,
           0.15,0.15,
           0.2,0.2), ncol=2, byrow=TRUE),
  # 8) Mixed Null/Alternative 2, non‐fixed control
  matrix(c(0.4,0.1,
           0.6,0.15,
           0.15,0.15,
           0.8,0.2), ncol=2, byrow=TRUE)
)

nScen  <- length(scenarios)   # 8 total scenarios
nSim   <- 1                  # replicates per scenario
nBask  <- 40                  # now 40 baskets per scenario
n1     <- 200                 # treatment patients per basket
n0     <- 200                 # control patients per basket

#--- Top‐level list: one element per scenario
simDataList <- vector("list", length = nScen)

for(s in seq_len(nScen)) {
  
  # Original 4 x 2 scenario definition
  scenarioDef <- scenarios[[s]]
  
  # Expand it to 40 x 2:
  #   first 10 rows use row 1's rates,
  #   next 10 rows use row 2's rates,
  #   ...
  #   last 10 rows use row 4's rates
  scenarioDef40 <- do.call(rbind, lapply(seq_len(4), function(i) {
    matrix(scenarioDef[i, ], nrow=10, ncol=2, byrow=TRUE)
  }))
  
  # Each scenario is a list of 'nSim' replicates
  simDataList[[s]] <- vector("list", length = nSim)
  
  for(r in seq_len(nSim)) {
    
    # This replicate's data frame: 40 rows (baskets) with outcome and rate columns
    dfRep <- data.frame(
      y1     = numeric(nBask),
      y0     = numeric(nBask),
      n1     = rep(n1, nBask),
      n0     = rep(n0, nBask),
      pi1    = numeric(nBask),
      pi0    = numeric(nBask),
      expEta = numeric(nBask)   # ratio = pi1 / pi0
    )
    
    # Fill in each of the 40 baskets
    for(b in seq_len(nBask)) {
      pi_t <- scenarioDef40[b, 1]  # treatment rate
      pi_c <- scenarioDef40[b, 2]  # control rate
      
      y1_draw <- rbinom(1, n1, pi_t)
      y0_draw <- rbinom(1, n0, pi_c)
      
      dfRep$y1[b]     <- y1_draw
      dfRep$y0[b]     <- y0_draw
      dfRep$pi1[b]    <- pi_t
      dfRep$pi0[b]    <- pi_c
      dfRep$expEta[b] <- pi_t / pi_c  # true exp(eta)
    }
    
    simDataList[[s]][[r]] <- dfRep
  }
}

#--- Now split into two final data sets:
#    1) The fixed‐control scenarios: 1,2,3
#    2) The non‐fixed‐control scenarios: 4,5,6,7,8
FXCR_simData    <- simDataList[1:3]
nonFXCR_simData <- simDataList[4:8]

# Quick checks
# Look at scenario 2 replicate 1:
# FXCR_simData[[3]][[1]]
#
# Look at scenario 8 replicate 10:
# nonFXCR_simData[[5]][[10]]
