psychopathROOTS <- function (x, modulus = TRUE) 
{
  K <- x$CALL$pars$n_symptoms
  if (x$CALL$pars$model %in% c("SVARHL","SVARHLX","SVARHLMA")){
    if (x$CALL$pars$model == "SVARHLMA") p <- x$fit$VARp else p <- x$fit$p
  }else{
    p <- x$CALL$pars$lag
  }

  if (x$CALL$pars$model=="DFM"){
    if (p==1){
      message("Roots are calculated for the equivalent VAR(1) process\n")
      A <- unlist(t(x$results$Dir_net[[1]]))
    }else{
      message("Roots are calculated for the factor process\n")
      A <- unlist(x$results$A_fact)
      
    }
  }else{
    A <- unlist(x$results$A)
  }
  companion <- matrix(0, nrow = K * p, ncol = K * p)
  companion[1:K, 1:(K * p)] <- A
  if (p > 1) {
    j <- 0
    for (i in (K + 1):(K * p)) {
      j <- j + 1
      companion[i, j] <- 1
    }
  }
  roots <- eigen(companion)$values
  if (modulus) 
    roots <- Mod(roots)
  return(roots)
}
