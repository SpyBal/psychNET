plot.pnt <- function(x, type, person, community, ...){
  models_all <- models_in()
  plot_list_result <- list()
  if (missing(person)) person <- NULL
  if (missing(type)) type <- "temporal"
  if (is.null(type)) type <- "temporal"
  if (missing(community)) community <- FALSE
  if ((!is.null(person)) && (x$CALL$pars$model != "MLVAR")) stop("The argument 'person' is used only when model='MLVAR")
  if (x$CALL$pars$model == "DFM"){
    return(plot_dfm(x, plot.type = type, communit = community, ...))
  }
  if (x$CALL$pars$model %in% c("VAR","SVAR","SVECM","SVARHL", "SVARHLX", "SVARHLMA","SMVAR")){
    return(plot_var(x, plot.type = type, communit = community,  ...))
  }
  if (x$CALL$pars$model %in% c("GVAR","GGVAR")){
    return(plot_gvar(x, plot.type = type, communit = community, ...))
  }
  if (x$CALL$pars$model == "MLVAR"){
    return(plot_mlvar(x, plot.type = type, individual = person, communit = community, ...))
  }
}


