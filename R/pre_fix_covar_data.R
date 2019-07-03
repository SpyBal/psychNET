pre_fix_covar_data <- function(x){
  data <- cbind(x$CALL$pars$covariates, x$CALL$pars$other_data)
  model <- x$CALL$pars$model
  cldata <- class(data)
  types <- c("data.frame", "matrix", "longitudinal")
  
  if (!cldata %in% types) stop("The data class must be 'matrix', 'data.frame', 'longitudinal', or 'list' with entries the pre-mentioned data classes")
  
  VARnames <- colnames(data)
  ncolumns <- ncol(data)
  nrows <- nrow(data)
  
  if (is.null(VARnames)) VARnames <- colnames(data) <- paste("symptom_",1:10,sep = "")
  
  if (cldata %in% types[1:2]){
    data <- as.data.frame(data)
    ct <- grepl("TIME", toupper(VARnames))
    cd <- grepl("DAY", toupper(VARnames))
    cb <- grepl("BEEP", toupper(VARnames))
    ci <- grepl("ID", toupper(VARnames))
    
    
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
  x$CALL$pars$other_data <- data[,no_sympt_dat_ind]
  
  
  x$CALL$pars$equidistant <- sapply(x$CALL$pars$ID_unique, function(y){
    tdifferences <- diff(x$CALL$pars$TIME[x$CALL$pars$ID==y])
    utdifferences <- unique(tdifferences)
    if (length(utdifferences)==1) return(TRUE) else return(FALSE) }) 
  
  
  if (any(!x$CALL$pars$equidistant) & (!model %in% c("MLVAR","SMVAR","GVAR") )){
    equi_style <- silver $ bold $ italic $ underline
    message(equi_style("\nMeasurements on covariates are transformed to consecutive by inserting NA's.... \n"))
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
    message(equi_style(paste("\nMissing values on covariates are imputed using:","Kalman.struct"," .... \n")))
    data <- impute_psychovar(data, method="Kalman.struct")
  }
  if (!is.null(x$CALL$pars$transform)){
    x$CALL$pars$TRANSFORM <- TRUE
    message(equi_style(paste("\nTime series are transformed using: 'Copula_skew' .... \n")))
    data <- transform_psychovar(data, method="Copula_skew")
  }
  
  if (any(is.na(data))  )stop("Check your transformation carefully in combination with the imputation method")
  data[,x$CALL$pars$symptoms]
}