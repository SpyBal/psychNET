K_smoother <- function(A, xitt, xittm, Ptt, Pttm, C, R, W) {
  T <- dim(xitt)[2]
  r <- dim(A)[1]
  
  Pttm <- Pttm[,,(1:(dim(Pttm)[3] - 1)), drop = FALSE]
  xittm <- xittm[,(1:(dim(xittm)[2] - 1)), drop = FALSE]
  
  # Whereas J is of constant dimension, L and K dimensions may vary
  # depending on existence of NAs
  J <- array(0, c(r, r, T))
  L <- list()
  K <- list()
  
  for (i in 1:(T-1)) {
    J[,,i] <- Ptt[,,i] %*% t(A) %*% solve(Pttm[,,(i+1)], tol = 1e-32)
  }
  
  Ci <- C 
  Ri <- R
  for (i in 1:T) {
    # Only keep entries for non-missing data
    C <- Ci[W[i,],, drop=FALSE]
    R <- Ri[W[i,], W[i,], drop=FALSE]
    L[[i]] <- solve(C %*% Pttm[,,i] %*% t(C) + R)
    K[[i]] <- Pttm[,,i] %*% t(C) %*% L[[i]]
  }
  
  xitT <- cbind(matrix(0, r, (T-1)), xitt[,T])
  PtT <- array(0, c(r, r, T))
  PtTm <- array(0, c(r, r, T))
  PtT[,,T] <- Ptt[,,T]
  PtTm[,,T] <- (diag(1, r) - K[[T]] %*% C) %*% A %*% Ptt[,,(T-1)]
  
  for (j in 1:(T-1)) {
    xitT[,(T-j)] <- xitt[,(T-j)] + J[,,(T-j)] %*% (xitT[,(T+1-j)] - xittm[,(T+1-j)])
    PtT[,,(T-j)] <- Ptt[,,(T-j)] + J[,,(T-j)] %*% (PtT[,,(T+1-j)] - Pttm[,,(T+1-j)]) %*% t(J[,,(T-j)])
  }
  
  for (j in 1:(T-2)) {
    PtTm[,,(T-j)] <- Ptt[,,(T-j)] %*% t(J[,,(T-j-1)]) + J[,,(T-j)] %*% (PtTm[,,(T-j+1)] - A %*% Ptt[,,(T-j)]) %*% t(J[,,(T-j-1)])
  }
  
  return(list(xitT = xitT, PtT = PtT, PtTm = PtTm))
}