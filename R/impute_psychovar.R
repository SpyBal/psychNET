impute_psychovar <- function(res, method="Kalman.struct"){
  x <- res$CALL$pars$data_fixed
  nsym <- res$CALL$pars$no_sympt_dat_ind
  ID <- NULL
  IDS <-  x[["ID"]]
  ncolumns <- ncol(x)
  nsymptoms <- ncolumns-4
  sympt.names <- res$CALL$pars$symptoms
  attrib_vec <- rep(NA,nsymptoms)
  UIDS <- unique(IDS)
  LUIDS <- length(UIDS)
  if (is.null(method)) method <-  res$CALL$pars$impute

  if (! (length(method)==1 | length(method) == nsymptoms) ) stop("Argument impute must be of length one or length equal to the number of symptoms.")
  
  if (length(method)==1){
    method <- rep(method, nsymptoms)
  }
  attrib_vec <- method
  names(attrib_vec) <- sympt.names
  label(x[,-nsym]) = lapply(sympt.names, function(x) attrib_vec[match(x, names(attrib_vec))])

  do.call(rbind, lapply(UIDS, function(y){
    dsub <- subset(x, x$ID==y)[,-nsym]
    dsub_rem <- subset(x, x$ID==y)[,nsym]
    
    cbind(sapply(dsub,  function(z){
      method <- as.character(label(z))
      
      if (method=="Kalman.arima"){
        z <- as.numeric(na.kalman(z,model = "auto.arima"))
      }
      if (method=="Kalman.struct"){
        z <- as.numeric(na.kalman(z,model = "StructTS"))
      }
      if (method=="Interpol.linear"){
        z <- as.numeric(na.interpolation(z, option = "linear"))
      }
      if (method=="Interpol.spline"){
        z <- as.numeric(na.interpolation(z, option = "spline"))
      }
      if (method=="Interpol.stine"){
        z <- as.numeric(na.interpolation(z, option = "stine"))
      }
      if (method=="LOCF"){
        z <- as.numeric(na.locf(z, option = "locf", na.remaining = "rev"))
      }
      if (method=="NOCB"){
        z <- as.numeric(imputeTS::na.locf(z, option = "nocb", na.remaining = "rev"))
      }
      if (method=="MA.simple"){
        z <- as.numeric(imputeTS::na.ma(z, k=4, weighting = "simple"))
      }
      if (method=="MA.linear"){
        z <- as.numeric(imputeTS::na.ma(z, k=4, weighting = "linear"))
      }
      if (method=="MA.exponent"){
        z <- as.numeric(imputeTS::na.ma(z, k=4, weighting = "exponential"))
      }
      if (method=="mean"){
        z <- as.numeric( imputeTS::na.mean(z, option = "mean"))
      }
      if (method=="mode"){
        z <- as.numeric(imputeTS::na.mean(z, option = "mode"))
      }
      if (method=="median"){
        z <- as.numeric(imputeTS::na.mean(z, option = "median"))
      }
      if (method=="random"){
        z <- as.numeric(na.random(z))
      }
      if (method=="Season.int.spline"){
        z <- as.numeric(na.seadec(z, algorithm = "interpolation", option = "spline"))
      }
      if (method=="Season.int.linear"){
        z <- as.numeric(na.seadec(z, algorithm = "interpolation", option = "linear"))
      }
      if (method=="Season.int.stine"){
        z <- as.numeric(na.seadec(z, algorithm = "interpolation", option = "stine"))
      }
      if (method=="Season.LOCF"){
        z <- as.numeric(na.seadec(z, algorithm = "locf", option = "locf", na.remaining = "rev"))
      }
      if (method=="Season.NOCB"){
        z <- as.numeric(na.seadec(z, algorithm = "locf", option = "nocb", na.remaining = "rev"))
      }
      if (method=="Season.mean"){
        z <- as.numeric(na.seadec(z, algorithm = "mean", option = "mean"))
      }
      if (method=="Season.median"){
        z <- as.numeric(na.seadec(z, algorithm = "mean", option = "median"))
      }
      if (method=="Season.mode"){
        z <- as.numeric(na.seadec(z, algorithm = "mean", option = "mode"))
      }
      if (method=="Season.kalman"){
        z <- as.numeric(na.seadec(z, algorithm = "kalman"))
      }
      if (method=="Season.random"){
        z <- as.numeric(na.seadec(z, algorithm = "random"))
      }
      if (method=="Season.MA"){
        z <- as.numeric(na.seadec(z, algorithm = "ma"))
      }
      return(z)
    }),dsub_rem)
  }))
  
}



