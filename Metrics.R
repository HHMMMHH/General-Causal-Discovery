library(SID)

calcSID <- function(A_true, A_est) {
  return(structIntervDist(A_true, A_est)$sid)
}

calcSHD <- function(A_true, A_est) {
  return(hammingDist(A_true, A_est))
}

#order performance
Dtop <- function(order, adj) {
  err <- 0
  n <- length(order)
  if (n < 2) return(err) 
  
  for (i in 1:(n - 1)) {
    err <- err + sum(adj[ order[(i + 1):n], order[i] ])
  }
  return(err)
}