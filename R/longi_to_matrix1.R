longi_to_matrix1  <- function(x){
  nmx <- colnames(x)
  dimx <- dim(x)
  p <- dimx[2]
  t <- dimx[1]
  x <- matrix(x,ncol = p,nrow = t)
  return(x)
}