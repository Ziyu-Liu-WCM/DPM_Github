signedLogSumExp <- function(signs, logs) {
  # signs: numeric vector of +1 or -1
  # logs: numeric vector (same length), logs[i] = log(abs(term_i))
  
  if (length(signs) == 0) return(-Inf)  # no terms => sum=0 => log(0)=-Inf
  stopifnot(length(signs) == length(logs))
  
  # 1) Sort by descending logs
  idx <- order(logs, decreasing = TRUE)
  sgn <- signs[idx]
  lval <- logs[idx]
  
  # 2) Initialize accumulator
  acc_sign <- 0       # 0 means "no sign yet"
  acc_logVal <- -Inf  # "scale" for the sum
  acc_val <- 0        # sum in scaled normal space
  i = 3
  # 3) Accumulate
  for (i in seq_along(sgn)) {
    if (acc_sign == 0) {
      # First nonempty term
      acc_sign <- sgn[i]
      acc_logVal <- lval[i]
      acc_val <- 1.0
    } else {
      # Scale current term by the difference in exponents
      delta <- sgn[i] * exp(lval[i] - acc_logVal)
      newSum <- acc_sign * acc_val + delta
      
      if (abs(newSum) < .Machine$double.eps) {
        # They nearly cancelled each other out:
        acc_sign <- 0
        acc_val <- 0
        acc_logVal <- -Inf
      } else if (newSum > 0) {
        acc_sign <- 1
        acc_val <- newSum
      } else {
        acc_sign <- -1
        acc_val <- -newSum
      }
    }
  }
  
  # 4) Final check: if everything cancelled
  if (acc_sign == 0 || acc_val == 0) {
    return(-Inf)      # sum is 0 => log(0) = -Inf
  } else if (acc_sign < 0) {
    return(NA) 
  } else {
    # sum = sign * exp(acc_logVal) * acc_val
    # => log(sum) = acc_logVal + log(acc_val)
    return(acc_logVal + log(acc_val))
  }
}

# signs <- c(-1, 1, 1)
# logs <- c(10, 9, 8)
# result <- signedLogSumExp(signs, logs)
# print(result)
# 
# 
# log(-exp(10) + exp(9) + exp(8))
