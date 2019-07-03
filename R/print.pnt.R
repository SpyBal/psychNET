print.pnt <- function(x, ...){
  
  models_all <- models_in()
  cat("\nCall: ")
  print(x$CALL$pars$call)
  cat("\n")
  cat("Time elapsed:",x$CALL$pars$time_elapsed, " minutes \n")
  cat("\n")
 
  cat("### Data: \n")
  cat("Number of persons:",x$CALL$pars$n_ID_unique, " \n")
  cat("Number of symptoms:",x$CALL$pars$n_symptoms, " \n")
  cat("Time points per person:",table(x$CALL$pars$data_pre_fixed$ID), " \n")
  if (!is.null(x$CALL$pars$transform)) cat("Transformed by:",x$CALL$pars$transform, " \n")
  cat("Missing values:",x$CALL$pars$miss, " \n")
  if (x$CALL$pars$miss) cat("Missing values are imputed by:",x$CALL$pars$impute, " \n")
  cat("Equidistant measurements:",x$CALL$pars$equidistant, " \n")
  if (!any(x$CALL$pars$equidistant) && !x$CALL$pars$model %in% c("SMVAR", "GVAR", "MLVAR")) cat("Consecutive imputed by:",x$CALL$pars$impute, " \n")
  cat("\n")
  
  cat("\n")
  cat("### Model: \n")
  cat("Model fitted:",x$CALL$pars$model, " \n")
  cat("Lag order:",x$CALL$pars$lag, " \n")
  if (x$CALL$pars$model == "DFM")    cat("Number of factors:",x$CALL$pars$nFact, " \n")
  if (x$CALL$pars$model %in% models_all$sparse){
    cat("Penalty type:",x$CALL$pars$penalty.type, " \n")
    cat("Model selection criterion:",x$CALL$pars$optimality, " \n")
  }   
  cat("\n")
  
  cat("\n")
  cat("### Estimates: \n")
  if (x$CALL$pars$model == "DFM" && x$CALL$pars$lag > 1){
    cat("Temporal network:",FALSE, " \n")
  }else{
    cat("Moduli of the roots of the autoregressive companion matrix:", round(psychopathROOTS(x),2),"\n")
    cat("Temporal network:",TRUE, " \n")
    if (x$CALL$pars$model %in% models_all$sparse){
      cat("Sparsity of temporal network:", sapply(1:x$CALL$pars$lag,function(i){
        mat <- x$results$Dir_net[[i]] 
        sum(mat == 0)/ prod(dim(mat))}), " \n")
    }
  }
  if (x$CALL$pars$model %in% c("GVAR","GGVAR","MLVAR") | (x$CALL$pars$model == "DFM" && x$CALL$pars$lag == 1)){
    cat("Contemporaneous network:",TRUE, " \n")
    if (x$CALL$pars$model %in% models_all$sparse){
      cat("Sparsity of contemporaneous network:", sum(x$results$UnDir_net == 0)/ prod(dim(x$results$UnDir_net)), " \n")
    }
  }else{
    cat("Contemporaneous network:",FALSE, " \n")
  }
  cat("\n")
  
}



