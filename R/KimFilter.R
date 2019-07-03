KimFilter <- function(x0, P0, y, F, A, R, Q, p)
{
  # Define all containers for further computations. Notations for variables and indices,
  # where appropriate, carefully follow Kim (1994). State vector is denoted as 'x', its
  # covariance as 'P'. Appended letters explicit whether these are updated, approximated
  # or smoothed.
  
  if (is.vector(y) || (ncol(y) <= length(x0)))
  {
    stop("Number of factors should be strictly lower than number of variables. \n
         Increase number of variables or estimate a VAR model instead.")
  }
  
  T <- nrow(y)
  n <- dim(F)[1]
  J <- length(x0)
  s <- dim(p)[1]
  
  ## x:   x^(i,j)_(t|t-1): predicted state vector - (2.6)
  ## xU:  x^(i,j)_(t|t): updated state vector - (2.11)
  ## P:   P^(i,j)_(t|t-1): predicted state covariance - (2.7)
  ## Pu:  P^(i,j)_(t|t): updated state covariance - (2.12)
  ## eta: eta^(i,j)_(t|t-1): conditional forecast error - (2.8)
  ## H:   H^(i,j)_(t): conditional variance of forecast error - (2.9)
  ## K:   K^(i,j)_(t): Kalman gain - (2.10)
  ## lik: f(y_t, S_(t-1)=i, S_t = j | t-1): joint conditional density - (2.16)
  ## loglik: log of (2.16)
  x <- array(NA, c(T,J,s,s))
  xU <- array(NA, c(T,J,s,s))
  P <- array(NA, c(T,J,J,s,s))
  Pu <- array(NA, c(T,J,J,s,s))
  eta <- array(NA, c(T,n,s,s))
  H <- array(NA, c(T,n,n,s,s))
  K <- array(NA, c(T,J,n,s,s))
  lik <- array(NA, c(T,s,s))
  loglik <- array(NA, c(T,s,s))
  ## Pr[S_(t-1) = i, S_t = j | t-1 ]: (2.15)
  ## Pr[S_(t-1) = i, S_t = j | t ]: (2.17)
  ## Pr[S_t = j | t-1 ]: used only for the smoothing part
  ## Pr[S_t = j | t ]: (2.18)
  jointP_fut <- array(NA, c(T,s,s))
  jointP_cur <- array(NA, c((T+1),s,s))
  stateP_fut <- array(NA, c(T,s))
  stateP <- array(NA, c(T,s))
  
  ## x^(j)_(t|t): approximate state vector conditional on S_j - (2.13)
  ## P^(j)_(t|t): approximate state covariance conditional on S_j - (2.14)
  xA <- array(NA, c(T,J,s))
  Pa <- array(0, c(T,J,J,s))
  result <- array(0, c(T,1))
  
  # Some initial conditions to get started
  for (i in 1:s) { xA[1,,i] <- x0 }
  for (i in 1:s) { Pa[1,,,i] <- P0 }
  jointP_cur[1,,] <- matrix(c(0.25,0.25,0.25,0.25), ncol=2)
  
  for (t in 2:T)
  {
    for (j in 1:s)
    {
      for (i in 1:s)
      {
        x[t,,i,j] <- A[,,j] %*% xA[(t-1),,i]
        P[t,,,i,j] <- A[,,j] %*% Pa[(t-1),,,i] %*% t(A[,,j]) + Q
        eta[t,,i,j] <- y[t,] - as(F[,,j], "matrix") %*% x[t,,i,j]
        H[t,,,i,j] <- F[,,j] %*% as(P[t,,,i,j], "matrix") %*% t(F[,,j]) + R
        K[t,,,i,j] <- P[t,,,i,j] %*% t(F[,,j]) %*% solve(H[t,,,i,j])
        xU[t,,i,j] <- x[t,,i,j] + K[t,,,i,j] %*% eta[t,,i,j]
        Pu[t,,,i,j] <- (diag(1,J) - K[t,,,i,j] %*% F[,,j]) %*% P[t,,,i,j]
        jointP_fut[t,i,j] <- p[i,j]*sum(jointP_cur[(t-1),,i]) # is everything alright here?
        lik[t,i,j] <- (2*pi)^(-n/2) * det(H[t,,,i,j])^(-1/2) *
          exp(-1/2*t(eta[t,,i,j]) %*% solve(H[t,,,i,j]) %*% eta[t,,i,j]) *
          jointP_fut[t,i,j]
        loglik[t,i,j] <- log(lik[t,i,j])
        jointP_cur[t,i,j] <- lik[t,i,j]
      }
      # Technically, there should be sum(lik[t,,]) term but it cancels out and is computed later
      stateP[t,j] <- sum(jointP_cur[t,,j])
      stateP_fut[t,j] <- sum(jointP_fut[t,,j])
      # Compute probability-filtered state process and its covariance
      xA[t,,j] <- xU[t,,,j] %*% jointP_cur[t,,j] / stateP[t,j]
      for (i in 1:s)
      {
        Pa[t,,,j] <- Pa[t,,,j] +
          (Pu[t,,,i,j] + (xA[t,,j] - xU[t,,i,j]) %*% t(xA[t,,j] - xU[t,,i,j])) * 
          exp(log(jointP_cur[t,i,j]) - log(stateP[t,j]))
      }
    }
    jointP_cur[t,,] <- exp(log(jointP_cur[t,,]) - log(sum(lik[t,,])))
    stateP[t,] <- exp(log(stateP[t,]) - log(sum(lik[t,,])))
    result[t,1] <- log(sum(lik[t,,]))
  }
  
  return(list("result"=sum(result), "xA"=xA, "Pa"=Pa, "x"=x, "P"=P, "stateP"=stateP, "stateP_fut"=stateP_fut))
  
  }