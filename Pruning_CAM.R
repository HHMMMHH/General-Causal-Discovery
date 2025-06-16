library(mgcv)
train_gam <-
  function(X,y,pars = list(numBasisFcts = 10))
  {
    if(!("numBasisFcts" %in% names(pars) ))
    { 
      pars$numBasisFcts = 10
    }
    p <- dim(as.matrix(X))
    if(p[1]/p[2] < 3*pars$numBasisFcts)
    {
      pars$numBasisFcts <- ceiling(p[1]/(3*p[2]))
      cat("changed number of basis functions to    ", pars$numBasisFcts, "    in order to have enough samples per basis function\n")
    }
    dat <- data.frame(as.matrix(y),as.matrix(X))
    coln <- rep("null",p[2]+1)
    for(i in 1:(p[2]+1))
    {
      coln[i] <- paste("var",i,sep="")
    }
    colnames(dat) <- coln
    labs<-"var1 ~ "
    if(p[2] > 1)
    {
      for(i in 2:p[2])
      {
        labs<-paste(labs,"s(var",i,",k = ",pars$numBasisFcts,") + ",sep="")
      }
    }
    labs<-paste(labs,"s(var",p[2]+1,",k = ",pars$numBasisFcts,")",sep="")
    mod_gam <- FALSE
    try(mod_gam <- gam(formula=formula(labs), data=dat),silent = TRUE)
    if(typeof(mod_gam) == "logical")
    {
      cat("There was some error with gam. The smoothing parameter is set to zero.\n")
      labs<-"var1 ~ "
      if(p[2] > 1)
      {
        for(i in 2:p[2])
        {
          labs<-paste(labs,"s(var",i,",k = ",pars$numBasisFcts,",sp=0) + ",sep="")
        }
      }
      labs<-paste(labs,"s(var",p[2]+1,",k = ",pars$numBasisFcts,",sp=0)",sep="")
      mod_gam <- gam(formula=formula(labs), data=dat)
    }
    result <- list()
    result$Yfit <- as.matrix(mod_gam$fitted.values)
    result$residuals <- as.matrix(mod_gam$residuals)
    result$model <- mod_gam 
    result$df <- mod_gam$df.residual     
    result$edf <- mod_gam$edf     
    result$edf1 <- mod_gam$edf1     
    
    # for degree of freedom see mod_gam$df.residual
    # for aic see mod_gam$aic
    return(result)
  }
selGam <-
  function(X,pars = list(cutOffPVal = 0.001, numBasisFcts = 10),output = FALSE,k)
  {
    result <- list()
    p <- dim(as.matrix(X))
    if(p[2] > 1)
    {
      selVec <- rep(FALSE, p[2])
      mod_gam <- train_gam(X[,-k],as.matrix(X[,k]),pars)
      pValVec <- summary.gam(mod_gam$model)$s.pv
      if(output)
      {
        cat("vector of p-values:", pValVec, "\n")
      }
      if(length(pValVec) != length(selVec[-k]))
      {
        show("This should never happen (function selGam).")
      }
      selVec[-k] <- (pValVec < pars$cutOffPVal)
    } else
    {
      selVec <- list()
    }
    return(selVec)
  }

pruning <-
  function(X, G, output = FALSE, pruneMethod = selGam, pruneMethodPars = list(cutOffPVal = 0.001, numBasisFcts = 10)) 
  {
    p <- dim(G)[1]
    finalG <- matrix(0,p,p)
    for(i in 1:p)
    {
      parents <- which(G[,i]==1)
      lenpa <- length(parents)
      
      if(output)
      {
        cat("pruning variable:", i, "\n")
        cat("considered parents:", parents, "\n")
      }
      
      if(lenpa>0)
      {
        Xtmp <- cbind(X[,parents],X[,i])
        selectedPar <- pruneMethod(Xtmp, k = lenpa + 1, pars = pruneMethodPars, output = output)
        finalParents <- parents[selectedPar]
        finalG[finalParents,i] <- 1
      }
    }
    
    return(finalG)
  }
pruning_CAM <- function(X, pi, CUTOFF = 0.001, VERBOSE = TRUE) {
  # Reorder columns according to pi
  X <- X[, pi]
  n <- length(pi)
  
  # Create initial adjacency matrix
  nodes <- as.character(pi)
  dag <- matrix(0, nrow = n, ncol = n)
  rownames(dag) <- nodes
  colnames(dag) <- nodes
  
  # Add edges: pi[i] → pi[j] for all j > i
  for (i in 1:(n - 1)) {
    from_node <- nodes[i]
    to_nodes <- nodes[(i + 1):n]
    dag[from_node, to_nodes] <- 1
  }
  
  # Pruning step
  pruned_dag <- pruning(
    X,
    dag,
    pruneMethod = selGam,
    pruneMethodPars = list(cutOffPVal = CUTOFF, numBasisFcts = 10),
    output = VERBOSE
  )
  
  # Restore adjacency matrix with original node labels (1:n)
  restored_dag <- matrix(0, nrow = n, ncol = n)
  rownames(restored_dag) <- as.character(1:n)
  colnames(restored_dag) <- as.character(1:n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      if (pruned_dag[i, j] == 1) {
        from <- pi[i]
        to   <- pi[j]
        restored_dag[as.character(from), as.character(to)] <- 1
      }
    }
  }
  
  return(restored_dag)
}

