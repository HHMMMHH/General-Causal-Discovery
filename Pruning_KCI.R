library(dHSIC)
library(CondIndTests)

pruning_KCI <- function(X, pi, alpha) {
  X <- scale(X)
  d <- ncol(X)
  A <- matrix(0, nrow = d, ncol = d)
  
  for (pos_j in seq_along(pi)) {
    j <- pi[pos_j]
    
    if (pos_j == 1) next
    
    Pre_j <- pi[1:(pos_j - 1)]
    
    for (i in Pre_j) {
      Z <- setdiff(Pre_j, i)
      
      if (length(Z) == 0) {
        test_result <- dhsic.test(X[, i], X[, j], method = "permutation",
                                  kernel = c("gaussian", "gaussian"),
                                  pairwise = FALSE, B = 1000)
        p_value <- test_result$p.value
      } else {
        test_result <- KCI(X[, i], X[, j], X[, Z], GP = FALSE)
        p_value <- test_result$pvalue
      }
      
      if (p_value < alpha) {
        A[i, j] <- 1
      }
    }
  }
  
  A
}
