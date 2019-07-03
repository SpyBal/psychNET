K_filter <- function(initx, initV, x, A, C, R, Q) {
  
  T <- dim(x)[1]
  N <- dim(x)[2]
  r <- dim(A)[1]
  W <- !is.na(x)    
  y <- t(x)
  
  xittm <- matrix(0, r, (T+1))
  xitt <- matrix(0, r, T)
  
  Pttm <- array(0, c(r, r, (T+1)))
  Ptt <- array(0, c(r, r, T))
  
  xittm[,1] <- initx
  Pttm[,,1] <- initV
  
  logl <- c()
  Ci <- C
  Ri <- R
  for (j in 1:T) {
    #      missing_data <- MissData(y[,j], C, R)
    #      C <- missing_data$C
    #      R <- missing_data$R
    C <- Ci[W[j,],, drop=FALSE]
    R <- Ri[W[j,], W[j,], drop=FALSE]
    if (FALSE) #(all(!W[j,])) #(all(is.na(missing_data$y) == TRUE))
    {
      xitt[,,j] <- A %*% xittm[,,j]
      Ptt[,,j] <- C %*% Pttm[,,j] %*% t(C) + R
    } else
    {
      # Innovation covariance (inverse)
      Icov <- C %*% Pttm[,,j] %*% t(C) + R
      L <- solve(Icov)
      # Innovation residual
      Ires <- as.numeric(na.omit(y[,j])) - C %*% xittm[,j]
      # Optimal Kalman gain
      G <- Pttm[,,j] %*% t(C) %*% L
      # Updated state estimate: predicted + (Kalman gain)*fitted
      xitt[,j] <- xittm[,j] + G %*% Ires
      # Updated covariance estimate
      Ptt[,,j] <- Pttm[,,j] - G %*% C %*% Pttm[,,j]
      # State space variable and covariance predictions E[f_t | t-1]
      xittm[,(j+1)] <- A %*% xitt[,j]
      Pttm[,,(j+1)] <- A %*% Ptt[,,j] %*% t(A) + Q
      
      # Compute log-likelihood with Mahalanobis distance
      d <- length(Ires)
      S <- C %*% Pttm[,,j] %*% t(C) + R
      Sinv <- solve(S)
      if (nrow(R) == 1)
      {
        GG <- t(C) %*% solve(R) %*% C
        detS <- prod(R) %*% det(diag(1, r) + Pttm[,,j] %*% GG)
      } else {
        GG <- t(C) %*% diag(1/diag(R)) %*% C
        detS <- prod(diag(R)) * det(diag(1, r) + Pttm[,,j] %*% GG)
      }
      denom <- (2 * pi)^(d/2) * sqrt(abs(detS))
      mahal <- sum(t(Ires) %*% Sinv %*% Ires)
      logl[j] <- -0.5 * mahal - log(denom)
    }
  }
  loglik <- sum(logl, na.rm=TRUE)
  return(list(xitt = xitt, xittm = xittm, Ptt = Ptt, Pttm = Pttm, loglik = loglik))
}