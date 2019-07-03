factorToVAR <- function(lambda, beta, psi, theta){
  if (missing(lambda) | missing(beta) | missing(psi) | missing(theta)){
    stop("'lambda', 'beta', 'psi' and 'theta' may not be missing.")
  }
  # Number of observed:
  nObs <- nrow(lambda)
  
  # Number of latents:
  nLat <- ncol(lambda)
  
  # Compute stationary latent variance:
  vecSigmaPsi <- solve(diag(nLat^2) - kronecker(beta, beta)) %*% c(psi)
  SigmaPsi <- matrix(vecSigmaPsi, nLat, nLat)
  
  # Compute stationary observed variance:
  SigmaY <- lambda %*% SigmaPsi %*% t(lambda) + theta
  
  # Implied beta of VAR model:
  BetaVAR <- lambda %*% beta %*% SigmaPsi %*% t(lambda) %*% solve(SigmaY)
  
  # Implied Psi of VAR Model:
  PsiVAR <- lambda %*% SigmaPsi %*% t(lambda) + theta - BetaVAR %*% SigmaY %*% t(BetaVAR)
  
  # Results object:
  Results <- list(
    beta = BetaVAR,
    psi = PsiVAR,
    PDC = computePDC(BetaVAR, solve(PsiVAR)),
    PCC = computePCC(solve(PsiVAR))
  )
  class(Results) <- "factorToVAR"
  return(Results)
}

