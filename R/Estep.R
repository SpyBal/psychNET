Estep <- function(y, A, C, Q, R, initx, initV, W) {
  
  os <- dim(y)[1]
  T <- dim(y)[2]
  ss <- nrow(A)
  
  kf <- K_filter(initx, initV, t(y), A, C, R, Q)
  ks <- K_smoother(A, kf$xitt, kf$xittm, kf$Ptt, kf$Pttm, C, R, W)
  
  xsmooth <- ks$xitT
  Vsmooth <- ks$PtT
  Wsmooth <- ks$PtTm
  
  delta <- matrix(0, os, ss)
  gamma <- matrix(0, ss, ss)
  beta <- matrix(0, ss, ss)
  
  for (t in 1:T) {
    z <- y[,t]; z[is.na(z)] <- 0
    # There seems to be a problem here
    delta <- delta + z %*% t(xsmooth[,t])
    gamma <- gamma + xsmooth[,t] %*% t(xsmooth[,t]) + Vsmooth[,,t]
    if (t > 1) {
      beta <- beta + xsmooth[,t] %*% t(xsmooth[,(t-1)]) + Wsmooth[,,t]
    }
  }
  
  gamma1 <- gamma - xsmooth[, T] %*% t(xsmooth[, T]) - Vsmooth[, , T]
  gamma2 <- gamma - xsmooth[, 1] %*% t(xsmooth[, 1]) - Vsmooth[, , 1]
  x1 <- xsmooth[, 1]
  V1 <- Vsmooth[, , 1]
  
  return(list(beta_t = beta, gamma_t = gamma, delta_t = delta, gamma1_t = gamma1,
              gamma2_t = gamma2, x1 = x1, V1 = V1, loglik_t = kf$loglik, xsmooth = xsmooth))
  
}