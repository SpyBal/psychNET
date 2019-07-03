transform_psychovar <- function(x, method){
  if (missing(method) | is.null(method)) method <- "Copula_skew"
  ID <- NULL
  IDS <-  x[["ID"]]
  ncolumns <- ncol(x)
  nsymptoms <- ncolumns-4
  sympt.names <- colnames(x[,1:nsymptoms])
  attrib_vec <- rep(NA,nsymptoms)
  UIDS <- unique(IDS)
  LUIDS <- length(UIDS)
  
  if (! (length(method)==1 | length(method) == nsymptoms) ) stop("Argument transform must be of length one or length equal to the number of symptoms.")
  
  if (length(method)==1){
    method <- rep(method, nsymptoms)
  }
  attrib_vec <- method
  names(attrib_vec) <- sympt.names
  label(x) = lapply(names(x), function(x) attrib_vec[match(x, names(attrib_vec))])
  
  do.call(rbind, lapply(UIDS, function(y){
    dsub <- subset(x, ID==y)[,1:nsymptoms]
    dsub_rem <- subset(x, ID==y)[,(nsymptoms+1):ncolumns]
    
    
    cbind(sapply(dsub,  function(z){
      method <- as.character(label(z))
      
      if (method=="log"){
        z <- as.numeric(log(z))
      }
      if (method=="log10"){
        z <- as.numeric(log10(z))
      }
      if (method=="Copula_discr"){
        z <- as.numeric(copula_transform(z))
      }
      if (method=="Copula_skew"){
        z <- as.numeric(copula_transform2(z))
      }
      if (method=="Zero.mean"){
        z <- as.numeric(scale(z,scale = FALSE))
      }
      if (method=="Standardize"){
        z <- as.numeric(scale(z))
      }
      if (method=="Power"){
        z <- z^powerTransform(z)$lambda
      }
      if (method=="Logit"){
        z <- logit(z)
      }
      if (method=="Square.root"){
        z <- sqrt(z)
      }
      if (method=="Power2"){
        z <- z*z
      }
      if (method=="Power3"){
        z <- z*z*z
      }
      if (method=="Cube.root"){
        z <- z^(1/3)
      }
      return(z)
    }),dsub_rem)
  }))
  
}



