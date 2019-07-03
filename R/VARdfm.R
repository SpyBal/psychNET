VARdfm <- function(x, p) {
  T <- nrow(x)
  Y <- x[(p + 1):T, ]
  X <- c()
  for (i in 1:p) {
    X <- cbind(X, x[(p + 1 - i):(T - i), ])
  }
  A <- solve(t(X) %*% X) %*% t(X) %*% Y
  res <- Y - X %*% A
  
  return(list(Y = Y, X = X, A = A, res = res))
}