fix_data <- function(x){
  data <- x$CALL$pars$data
  types <- c("data.frame", "matrix", "longitudinal", "list")
  cldata <- class(data)
  if (!any(cldata %in% types)) stop("The data class must be 'matrix', 'data.frame', 'longitudinal', or 'list' with entries the abovementioned data classes")
  if (cldata %in% types[1:3]){
    x <- pre_fix_data(x)
    if (!is.null(x$CALL$pars$covariates)){
      x$CALL$pars$covariates <- pre_fix_covar_data(x)
    }
  }else{
    ndata <- length(data)
    nlistx <- lapply(1:ndata, function(i) x)
    for (j in 1:ndata) nlistx[[j]]$CALL$pars$data <- data[[j]]
    x$Multifix_list <- lapply(nlistx,pre_fix_data)
  }
  x
}