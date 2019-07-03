em_converged <- function(loglik, previous_loglik, threshold=1e-4, check_increased=TRUE) {
  
  converged <- FALSE
  decrease <- 0
  
  if (check_increased == TRUE) {
    if (loglik - previous_loglik < -0.001) {
      #            cat("*** Likelihood decreased from ", previous_loglik, " to ", loglik, "\n")
      decrease <- 1
    }
  }
  
  delta_loglik <- abs(loglik - previous_loglik)
  avg_loglik <- (abs(loglik) + abs(previous_loglik) + .Machine$double.eps)/2
  
  if ((delta_loglik/avg_loglik) < threshold) {
    converged <- TRUE
  }
  return(converged)
  # return(list('converged'=converged, 'decrease'=decrease))
}