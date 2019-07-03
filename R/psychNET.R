psychNET <- function(data, model, lag, criterion, nFact, penalty, lambda1,lambda2, 
                      covariates, impute, transform, ... ){
  start_time <- as.POSIXct(Sys.time())
  if(missing(data)) stop("No 'data' are provided")
  if(missing(model)) stop("Argument 'model' must be specified")
  if(missing(lag) & !model %in% c("SVARHL", "SVARHLX", "SVARHLMA")) lag <- 1
  if(missing(lag) & model %in% c("SVARHL", "SVARHLX", "SVARHLMA")) lag <- NULL
  if(missing(criterion)) criterion <- NULL
  if(missing(nFact)) nFact <- NULL
  if(missing(penalty)) penalty <- NULL
  if(missing(lambda1)) lambda1 <- NULL
  if(missing(lambda2)) lambda2 <- NULL
  if(missing(covariates)) covariates <- NULL
  if(missing(transform)) transform <- NULL
  if(missing(impute)) impute <- "Kalman.struct"
  optimality <- criterion
  penalty.type <- penalty
  three_dot_args <- list(...)
  result <- list()
  result$CALL <- list()
  result$CALL$pars <- list()
  result$CALL$pars$call <- match.call()
  result$CALL$pars$data <- data
  result$CALL$pars$nFact <- nFact
  result$CALL$pars$model <- model
  result$CALL$pars$optimality <- optimality
  result$CALL$pars$lambda1 <- lambda1
  result$CALL$pars$lambda2 <- lambda2
  result$CALL$pars$penalty.type <- penalty.type
  result$CALL$pars$lag <- lag
  result$CALL$pars$covariates <- covariates
  result$CALL$pars$impute <- impute
  result$CALL$pars$transform <- transform
  result$CALL$pars$dots <- three_dot_args
  ## Initial check of the data
  result <- fix_data(result)
  result <- check_arguments(result)
  result$data <- transform_model_data(result)
  result  <- psychofit(result,...)
  class(result) <- "pnt"
  end_time <- as.POSIXct(Sys.time())
  result$CALL$pars$time_elapsed <- difftime(end_time, start_time, units = "mins")
  return(result)
}