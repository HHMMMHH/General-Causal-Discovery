k_rbf <- function(x, y, sigma) {
  exp(-sum((x - y)^2) / (2 * sigma^2))
}

compute_kernel_matrix <- function(kernel_func, data, sigma) {
  n <- nrow(data)
  K <- matrix(0, n, n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      K[i, j] <- kernel_func(data[i, , drop = FALSE], data[j, , drop = FALSE], sigma)
    }
  }
  K
}

tilde <- function(K) {
  n <- nrow(K)
  H <- diag(n) - matrix(1, n, n) / n
  H %*% K %*% H
}

##### GCD_star and GCD represent the deflation and no-deflation versions of the code, respectively.
##### For bivariate experiments, data should be standardized using `scale()`. 
##### For multivariate experiments, however, we do not recommend applying `scale()` for global standardization, 
##### as this process may distort the true proportional relationships among variables, potentially affecting the identification of causal orders and structures.
GCD_star <- function(data, sigma1 = 0.4, sigma2 = 0.2, lambda = 0.01) {
  n  <- nrow(data)
  d  <- ncol(data)
  causal_order <- integer(0)
  if (d == 2) {
    current_data <- scale(as.matrix(data))
  } else {
    current_data <- as.matrix(data)
  }
  S <- seq_len(d)
  
  while (length(S) > 0) {
    K1_centered <- R_list <- vector("list", length(S))
    names(K1_centered) <- names(R_list) <- S
    
    for (idx in seq_along(S)) {
      j <- S[idx]
      vec_j <- matrix(current_data[, j], ncol = 1)
      
      K1 <- compute_kernel_matrix(k_rbf, vec_j, sigma1)
      K2_cen <- tilde(compute_kernel_matrix(k_rbf, vec_j, sigma2))
      
      K1_centered[[idx]] <- tilde(K1)
      R_list[[idx]] <- lambda * solve(K2_cen + lambda * diag(n))
    }
    
    test_mat <- matrix(0, nrow = length(S), ncol = length(S), dimnames = list(S, S))
    
    for (i in S) {
      for (j in S) {
        if (i != j) {
          JI <- R_list[[as.character(i)]] %*%
            K1_centered[[as.character(j)]] %*%
            t(R_list[[as.character(i)]])
          I <- K1_centered[[as.character(i)]]
          test_mat[as.character(i), as.character(j)] <- sum(diag(JI %*% I))
        }
      }
    }
    
    sum_vals <- sapply(S, function(i) {
      sum(test_mat[as.character(i), as.character(setdiff(S, i))])
    })
    i_star <- S[which.min(sum_vals)]
    causal_order <- c(causal_order, i_star)
    
    S <- setdiff(S, i_star)
    
    if (length(S) > 0) {
      R_star <- R_list[[as.character(i_star)]]
      for (j in S) {
        res_j <- R_star %*% matrix(current_data[, j], ncol = 1)
        current_data[, j] <- res_j[, 1]
      }
    }
  }
  causal_order
}

GCD <- function(data, sigma1 = 0.4, sigma2 = 0.2, lambda = 0.01) {
  n  <- nrow(data)
  d  <- ncol(data)
  if (d == 2) {
    data <- scale(as.matrix(data))
  }
  causal_order <- integer(0)
  S <- seq_len(d)
  
  K1_centered <-R_list <- vector("list", d)
  names(K1_centered) <- names(R_list) <- S
  
  for (idx in seq_along(S)) {
    j <- S[idx]
    vec_j <- matrix(data[, j], ncol = 1)
    
    K1 <- compute_kernel_matrix(k_rbf, vec_j, sigma1)
    K2_cen <- tilde(compute_kernel_matrix(k_rbf, vec_j, sigma2))

    K1_centered[[idx]] <- tilde(K1)
    R_list[[idx]] <- lambda * solve(K2_cen + lambda * diag(n))
  }
  
  test_mat <- matrix(0, nrow = d, ncol = d, dimnames = list(S, S))
  
  for (i in S) {
    for (j in S) {
      if (i != j) {
        JI <- R_list[[as.character(i)]] %*%
          K1_centered[[as.character(j)]] %*%
          t(R_list[[as.character(i)]])
        I <- K1_centered[[as.character(i)]]
        test_mat[as.character(i), as.character(j)] <- sum(diag(JI %*% I))
      }
    }
  }
  
  while (length(S) > 0) {
    sum_vals <- sapply(S, function(i) {
      sum(test_mat[as.character(i), as.character(setdiff(S, i))])
    })
    i_star <- S[which.min(sum_vals)]
    causal_order <- c(causal_order, i_star)
    
    S <- setdiff(S, i_star)
  }
  
  causal_order
}

