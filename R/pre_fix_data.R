pre_fix_data <- function(x){
  data <- x$CALL$pars$data
  model <- x$CALL$pars$model
  x$CALL$pars$miss <- any(is.na(data))
  
  cldata <- class(data)
  types <- c("data.frame", "matrix", "longitudinal")
  
  if (!cldata %in% types) stop("The data class must be 'matrix', 'data.frame', 'longitudinal', or 'list' with entries the pre-mentioned data classes")

  VARnames <- colnames(data)
  ncolumns <- ncol(data)
  nrows <- nrow(data)
  
  if (is.null(VARnames)) VARnames <- colnames(data) <- paste("V_",1:ncolumns,sep = "")
    
  if (cldata %in% types[1:2]){
    data <- as.data.frame(data)
    ct <- grepl("^TIME$", toupper(VARnames))
    cd <- grepl("^DAY$", toupper(VARnames))
    cb <- grepl("^BEEP$", toupper(VARnames))
    ci <- grepl("^ID$", toupper(VARnames))
    
  
    if (any(ci==TRUE)){
      x$CALL$pars$ID_index <- which(ci==TRUE)
      x$CALL$pars$ID <- as.numeric(data[,x$CALL$pars$ID_index])
    }else{
      x$CALL$pars$ID_index <- ncolumns + 1 
      x$CALL$pars$ID <- data[["ID"]] <- rep(1,nrows)
    }
    
    if (any(cd==TRUE)){
      x$CALL$pars$DAY_index <- which(cd==TRUE)
      x$CALL$pars$DAY <- as.numeric(data[,x$CALL$pars$DAY_index])
    }else{
      x$CALL$pars$DAY_index <- ncol(data) + 1
      x$CALL$pars$DAY <- data[["DAY"]] <- rep(1,nrows)
    }
    
    if (any(cb==TRUE)){
      x$CALL$pars$BEEP_index <- which(cb==TRUE)
      x$CALL$pars$BEEP <-  data[,x$CALL$pars$BEEP_index]
    }else{
      x$CALL$pars$BEEP_index <- ncol(data) + 1
      x$CALL$pars$BEEP <- data[["BEEP"]] <- do.call(c,lapply(table(x$CALL$pars$ID), function(y) 1:y))
    }
    
    if (any(ct==TRUE)){
      x$CALL$pars$TIME_index <- which(ct==TRUE)
      x$CALL$pars$TIME <-  data[,x$CALL$pars$TIME_index]
    }else{
      x$CALL$pars$TIME_index <- ncol(data) + 1
      x$CALL$pars$TIME <- data[["TIME"]] <- do.call(c,lapply(unique(x$CALL$pars$ID),function(y){
        ind <- which(x$CALL$pars$ID==y)
        beepday2consec(x$CALL$pars$BEEP[ind], x$CALL$pars$DAY[ind])
      } ))
    }
    
  } 
  
  if (cldata == "longitudinal"){
    data <- longi_to_dataframe(data)
    x$CALL$pars$TIME_index <- ncolumns + 2
    x$CALL$pars$TIME <- data[,x$CALL$pars$TIME_index]
    
    x$CALL$pars$BEEP_index <- ncolumns + 4
    x$CALL$pars$BEEP <- data[,x$CALL$pars$BEEP_index]
    
    x$CALL$pars$DAY_index <- ncolumns + 3
    x$CALL$pars$DAY <- data[,x$CALL$pars$DAY_index]
    
    x$CALL$pars$ID_index <- ncolumns + 1
    x$CALL$pars$ID <- data[,x$CALL$pars$ID_index]
  } 
  
  x$CALL$pars$ID_unique <- unique(x$CALL$pars$ID)
  x$CALL$pars$n_ID_unique <- length(x$CALL$pars$ID_unique)
  no_sympt_dat_ind <- c(x$CALL$pars$ID_index, x$CALL$pars$DAY_index,x$CALL$pars$BEEP_index,x$CALL$pars$TIME_index )
  x$CALL$pars$sympt_data <- data[,-no_sympt_dat_ind]
  x$CALL$pars$symptoms <- colnames(x$CALL$pars$sympt_data)
  x$CALL$pars$n_symptoms <- ncol(x$CALL$pars$sympt_data)
  x$CALL$pars$no_sympt_dat_ind <- no_sympt_dat_ind
  x$CALL$pars$other_data <- data[,no_sympt_dat_ind]
  
  if (x$CALL$pars$model== "SMVAR"){
    test_dat <- data 
    whichINTEG <- whichNUMER <-  whichFACTOR <- whichOFACTOR <- NULL
    type_data_each <- split(colnames(test_dat),sapply(test_dat, function(x) paste(class(x), collapse=" ")))
    whichNUMER <- colnames(test_dat) %in% type_data_each$numeric
    whichINTEG <- colnames(test_dat) %in% type_data_each$integer
    whichFACTOR <- colnames(test_dat) %in% type_data_each$factor
    whichOFACTOR <- colnames(test_dat) %in% type_data_each$`ordered factor`
    levels_var <- rep(NA, ncol(test_dat))
    type_var <- rep(NA,ncol(test_dat))
    if (!is.null(whichNUMER)){
      type_var[whichNUMER] <- "g"
      levels_var[whichNUMER] <- 1
    }
    if (!is.null(whichINTEG)){
      type_var[whichINTEG] <- "p"
      if (sum(whichINTEG) == 1) levels_var[whichINTEG] <- length(unique(test_dat[,whichINTEG])) 
      if (sum(whichINTEG) > 1) levels_var[whichINTEG] <- sapply(test_dat[,whichINTEG] ,function(x) length(unique(na.omit(x))))
    }
    if (!is.null(whichFACTOR)){
      type_var[whichFACTOR] <- "c"
      if (sum(whichFACTOR) == 1) levels_var[whichFACTOR] <- length(unique(test_dat[,whichFACTOR])) 
      if (sum(whichFACTOR) > 1) levels_var[whichFACTOR] <- sapply(test_dat[,whichFACTOR] ,function(x) length(unique(na.omit(x))))
    }
    if (!is.null(whichOFACTOR)){
      type_var[whichOFACTOR] <- "o"
      if (sum(whichOFACTOR) == 1) levels_var[whichOFACTOR] <- length(unique(test_dat[,whichOFACTOR])) 
      if (sum(whichOFACTOR) > 1) levels_var[whichOFACTOR] <- sapply(test_dat[,whichOFACTOR] ,function(x) length(unique(na.omit(x))))
    }
    additional_args <- list(type = type_var,level = levels_var)
    x$CALL$pars$dots <- additional_args
  }
  
  
  x$CALL$pars$equidistant <- sapply(x$CALL$pars$ID_unique, function(y){
    tdifferences <- diff(x$CALL$pars$TIME[x$CALL$pars$ID==y])
    utdifferences <- unique(tdifferences)
    if (length(utdifferences)==1) return(TRUE) else return(FALSE) }) 
  
  
  if (any(!x$CALL$pars$equidistant) & (!model %in% c("MLVAR","SMVAR","GVAR") )){
    equi_style <- silver $ bold $ italic $ underline
    message(equi_style("\nMeasurements on the dependent variables are transformed to consecutive by inserting NA's.... \n"))
    x$CALL$pars$ind_nequidistant <- which(!x$CALL$pars$equidistant)
    ids_neq <- x$CALL$pars$ID_unique[x$CALL$pars$ind_nequidistant]
    data <- do.call(rbind,lapply(x$CALL$pars$ID_unique, function(y){
      ind <- x$CALL$pars$ID == y
      time_ind <- x$CALL$pars$TIME[ind]
      mx_time_ind <- max(time_ind)
      datmat <- data.frame(matrix(NA, ncol = x$CALL$pars$n_symptoms, nrow=mx_time_ind))
      colnames(datmat) <- x$CALL$pars$symptoms
      datmat[["ID"]] <- rep(y, mx_time_ind)
      x$CALL$pars$ID_index <- x$CALL$pars$n_symptoms + 1
      
      datmat[["TIME"]] <- 1:mx_time_ind
      x$CALL$pars$TIME_index <- x$CALL$pars$n_symptoms + 2
      
      datmat[["BEEP"]] <- 1:mx_time_ind
      x$CALL$pars$BEEP_index <- x$CALL$pars$n_symptoms + 3
      
      datmat[["DAY"]] <- rep(1,mx_time_ind)
      x$CALL$pars$DAY_index <- x$CALL$pars$n_symptoms + 4
      
      datmat[time_ind,1:x$CALL$pars$n_symptoms] <- x$CALL$pars$sympt_data[time_ind,]
      datmat}))
  }   
  
  x$CALL$pars$data_fixed <- data
  
  if (any(is.na(data))){
    x$CALL$pars$MISSINGS <- TRUE
    equi_style <- silver $ bold $ italic $ underline
    message(equi_style(paste("\nMissing values on the dependent variables are imputed using:",x$CALL$pars$impute," .... \n")))
    if (x$CALL$pars$model== "SMVAR"){
      data <- data.frame(sapply(data, function(x) as.numeric(as.character(x))))
      x$CALL$pars$data_fixed <- data
      
    } 
    data <- impute_psychovar(x, method=x$CALL$pars$impute)
    
  }
  if (!is.null(x$CALL$pars$transform)){
    x$CALL$pars$TRANSFORM <- TRUE
    equi_style <- silver $ bold $ italic $ underline
    message(equi_style(paste("\nTime series are transformed using:",x$CALL$pars$transform," .... \n")))
    data <- transform_psychovar(data, method=x$CALL$pars$transform)
  }
  
  if (any(is.na(data))  )stop("Check your transformation carefully in combination with the imputation method")
  
  if (x$CALL$pars$model== "SMVAR"){
    data <- as.data.frame(data)
    if (x$CALL$pars$model=="SMVAR") data <- data else data <- data.frame(sapply(data, function(x) as.numeric(as.character(x))))
  } 
  x$CALL$pars$data_pre_fixed <- data
  x
}
