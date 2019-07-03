data_to_longi <- function(data){
  itm_names <- colnames(data)
  data_type1 <- is.matrix(data)
  data_type2 <- is.data.frame(data)
  GROUPind <- which(colnames(data)=="GROUP")
  TIMEind <- which(colnames(data)=="TIME")
  DAYind <- which(colnames(data)=="DAY")
  BEEPind <- which(colnames(data)=="BEEP")
  if (length(which(colnames(data)=="ID")) == 0){
    data <- cbind(data,rep(1,nrow(data)))
    colnames(data)[ncol(data)] <- "ID"
  } 
  IDind <- which(colnames(data)=="ID")
  no_sympt <- c(IDind, GROUPind, TIMEind, DAYind, BEEPind)
  if (!(data_type1 | data_type2)) stop("The data should be a matrix or a data.frame")
  idvar <- "ID"
  datids <- as.numeric(data[,idvar])
  IDS <- unique(datids)
  n <- length(IDS)
  ns <- ncol(data[,-no_sympt])
  Ti <- sapply(1:n, function(i) nrow(subset(data, datids==IDS[i])))
  check1 <- diff(Ti)
  if (sum(check1 != 0) > 0){
    stop("All units must have the same number of observations")
  }
  TT <- Ti[1]
  lonD <- matrix(NA, nrow = n*TT, ncol = ns )
  d <- data
  for (i in 1:TT){
    d <- lapply(1:n, function(k){
      index <- (data[,"ID"]== unique(data[,"ID"])[k])
      dat <- data[index, ]
      dat[i,]
    })
    dd <- do.call(rbind,d)[,-no_sympt]
    if (n>1) dimnames(dd)[[1]] <- 1:n + (i-1)*n
    lonD[(1:n + (i-1)*n),] <- data.matrix(dd)
  }
  LONGI_OBJ <- as.longitudinal(lonD, repeats = n)
  colnames(LONGI_OBJ) <- itm_names[-no_sympt]
  return(LONGI_OBJ)
}
