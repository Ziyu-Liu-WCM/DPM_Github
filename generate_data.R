set.seed(1234)

#--- 8 scenarios, each is a 4x2 matrix: columns=(treatment_rate, control_rate), rows=baskets
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

nScen  <- length(scenarios)  # 8 total
nSim   <- 200               # replicates per scenario
nBask  <- 4                  # baskets per scenario
n1     <- 20                 # treatment patients per basket
n0     <- 20                # control patients per basket

#--- Top‐level list: one element per scenario
simDataList <- vector("list", length=nScen)

for(s in seq_len(nScen)) {
  
  # Each scenario is a list of 'nSim' replicates
  simDataList[[s]] <- vector("list", length=nSim)
  
  # 4x2 matrix: row = basket, col = (treatment_rate, control_rate)
  scenarioDef <- scenarios[[s]]
  
  for(r in seq_len(nSim)) {
    
    # This replicate's data frame: 4 rows (baskets), and we add ratio col
    dfRep <- data.frame(
      y1   = numeric(nBask),
      y0   = numeric(nBask),
      n1   = rep(n1, nBask),
      n0   = rep(n0, nBask),
      pi1  = numeric(nBask),
      pi0  = numeric(nBask),
      expEta= numeric(nBask)    # ratio = pi1 / pi0
    )
    
    for(b in seq_len(nBask)) {
      pi_t <- scenarioDef[b, 1]  # treatment rate
      pi_c <- scenarioDef[b, 2]  # control rate
      
      y1_draw <- rbinom(1, n1, pi_t)
      y0_draw <- rbinom(1, n0, pi_c)
      
      dfRep$y1[b]   = y1_draw
      dfRep$y0[b]   = y0_draw
      dfRep$pi1[b]  = pi_t
      dfRep$pi0[b]  = pi_c
      dfRep$expEta[b]= pi_t / pi_c  # true exp(eta)
    }
    
    simDataList[[s]][[r]] <- dfRep
  }
}

#--- Now split into two final data sets:
#    1) The fixed‐control scenarios: 1,2,3
#    2) The non‐fixed‐control scenarios: 4,5,6,7,8
FXCR_simData    <- simDataList[1:3]
nonFXCR_simData <- simDataList[4:8]

# Optional: quick checks
# Look at scenario 2 replicate 1:
# simDataListFixed[[2]][[1]]
# 
# # Look at scenario 8 replicate 999:
# simDataListNonFixed[[5]][[999]]  # (5 means scenario 4->8 => 4=1,5=2,... so scenario8 is index5 here)