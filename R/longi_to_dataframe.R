longi_to_dataframe  <- function(x){
  nmx <- colnames(x)
  dimx <- dim(x)
  p <- dimx[2]
  nt <- dimx[1]
  struc <- get.time.repeats(x)
  rr <- struc$repeats
  n <- unique(rr)
  if (length(n) != 1) stop("Data is unbalanced. For multiple subjects, data for each subject must have the same dimensions")
  t <- nt / n
  tau <- struc$time
  time <- rep(tau,unique(rr))
  ID <- rep(1:n,each=length(tau))
  x <- matrix(x,ncol = p,nrow = nt)
  dt_list <- list()
  for (i in 1:n) dt_list[[i]] <- x[seq.int(i,nt,by=n),]
  df <- data.frame(do.call(rbind,dt_list))
  names(df) <- nmx
  df$ID <- ID
  df$TIME <- time
  df$DAY <- rep(1,nrow(df))
  df$BEEP <- time
  return(df)
}


