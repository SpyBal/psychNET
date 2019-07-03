KimSmoother <- function(xA, Pa, A, P, x, p, stateP, stateP_fut)
{
  # Define all containers for further computations. Notations for variables and indices,
  # where appropriate, carefully follow Kim (1994). State vector is denoted as 'x', its
  # covariance as 'P'. Appended letters explicit whether these are updated, approximated
  # or smoothed.
  
  T <- dim(xA)[1]
  J <- dim(xA)[2]
  s <- dim(xA)[3]
  
  ## Pr[S_t = j, S_(t+1) = k | T]: (2.20)
  ## Pr[S_t = j | T]: (2.21)
  jointPs <- array(NA, c(T,s,s))
  ProbS <- array(NA, c(T,s))
  
  ## xS: x^(j,k)_(t|T): inference of x_t based on full sample - (2.24)
  ## Ps: P^(j,k)_(t|T): covariance matrix of x^(j,k)_(t|T) - (2.25)
  ## Ptilde: helper matrix as defined after (2.25) 
  xS <- array(0, c(T,J,s,s))
  Ps <- array(0, c(T,J,J,s,s))
  Ptilde <- array(NA, c(T,J,J,s,s))
  
  ## xAS: x^(j)_(t|T): smoothed and approximated state vector conditional on S_j (2.26)
  ## Pas: P^(j)_(t|T): smoothed and approximated state covariance conditional on S_j (2.27)
  ## xF: x_(t|T): state-weighted [F]inal state vector (2.28)
  ## Pf: P_(t|T): state-weighted [f]inal state covariance
  xAS <- array(0, c(T,J,s))
  Pas <- array(0, c(T,J,J,s))
  xF <- array(0, c(T,J))
  Pf <- array(0, c(T,J,J))
  # Initial conditions for smoothing loop
  ProbS[T,] <- stateP[T,]
  
  for (t in seq(T-1,1,-1))
  {
    for (j in 1:s)
    {
      for (k in 1:s)
      {
        jointPs[t,j,k] <- ProbS[(t+1),k]*stateP[t,j]*p[j,k] / stateP_fut[(t+1),k]
        Ptilde[t,,,j,k] <- Pa[t,,,j] %*% t(A[,,k]) %*% solve(P[(t+1),,,j,k])
        xS[t,,j,k] <- xA[t,,j] + Ptilde[t,,,j,k] %*% (xA[(t+1),,k] - x[(t+1),,j,k])
        Ps[t,,,j,k] <- Pa[t,,,j] +
          Ptilde[t,,,j,k] %*% (Pa[(t+1),,,k] - P[(t+1),,,j,k]) %*% t(Ptilde[t,,,j,k])
        xAS[t,,j] <- xAS[t,,j] + jointPs[t,j,k]*xS[t,,j,k]
        Pas[t,,,j] <- Pas[t,,,j] + jointPs[t,j,k]*(Ps[t,,,j,k] +
                                                     (xAS[t,,j] - xS[t,,j,k]) %*% t(xAS[t,,j] - xS[t,,j,k]))
      }
      ProbS[t,j] <- sum(jointPs[t,j,])
      xAS[t,,j] <- xAS[t,,j] / ProbS[t,j]
      Pas[t,,,j] <- Pas[t,,,j] / ProbS[t,j]
    }
  }
  for (t in 1:T)
  {
    for (j in 1:s)
    {
      xF[t,] <- xF[t,] + xAS[t,,j]*ProbS[t,j]
      Pf[t,,] <- Pf[t,,] + Pas[t,,,j]*ProbS[t,j]
    }
  }
  
  return(list("xF"=xF, "Pf"=Pf, "ProbS"=ProbS))
}