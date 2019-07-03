calcNeighbors2 <- function (fit, lambda, type, level, v) 
{
  if (type[v] == "o") n_cats <- fit$nLev else n_cats <- level[v]
  if (! type[v] %in% c("o","c")) {
    coefs_bin <- as.matrix(coef(fit, s = lambda)[-1, ]) != 0
    n_neighbors <- colSums(coefs_bin)
  }
  if (type[v] == "c") {
    m_neighbors <- matrix(0, ncol = length(fit$lambda), 
                          nrow = n_cats)
    coefs_bin <- vector("list", length = n_cats)
    for (ca in 1:n_cats) {
      coefs_bin[[ca]] <- as.matrix(coef(fit, s = lambda)[[ca]][-1, 
                                                               ]) != 0
    }
    n_neighbors <- colSums(Reduce("+", coefs_bin) != 0)
  }

  if (type[v] == "o") {
    m_neighbors <- matrix(0, ncol = length(fit$lambdaVals), nrow = (n_cats-1))
    coefs <- coef(fit, whichLambda = which(fit$lambdaVals == lambda), matrix=TRUE)
    coefs_bin <- vector("list", length = n_cats-1)
    coefs_new_bin <- vector("list", length = n_cats-1)
    coefs_bin <- as.matrix(coefs[-1,1] != 0)
    new_bin <- vector("logical", length = length(level))
    names(new_bin) <-  names(level)
    for (variab in names(level)){
      if (any(coefs[grepl(variab, rownames(coefs)),][,1] != 0)) new_bin[variab] <- TRUE else new_bin[variab] <- FALSE
    }
    n_neighbors <- sum(new_bin)
  }
  return(n_neighbors)
}


