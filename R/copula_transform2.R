copula_transform2 <- function(x){
  cdf <- (rank(x)-.5)/(length(x))
  qnorm(cdf)
}

