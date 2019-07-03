summary.dfm <- function(x, plot=FALSE) {
  cl <- match.call()
  nf <- dim(x$qml)[2]
  if (plot == TRUE)
  {
    for (i in 1:nf)
    {
      plot.title <- paste0("QML estimated factor ", i)
      plot(x$qml[,i], type='l', main=plot.title, ylab="Value", xlab="Time") 
    }
    boxplot(x$data - t(x$C %*% t(x$qml)), main="Residuals by input variable")
  }
  cat("Observation equation matrix: \n")
  print(x$C)
  
  cat("\nObservation residual covariance matrix: \n")
  print(cov(x$data - t(x$C %*% t(x$qml))))
  
  cat("\nSystem equation transition matrix: \n")
  print(x$A)
  
}