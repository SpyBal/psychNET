copula_transform <- function(x){
  cdf <- (rank(jitter(x))-.5)/(length(x))
  qnorm(cdf)
}