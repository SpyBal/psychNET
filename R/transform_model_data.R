
transform_model_data <- function(res){
  data_type <- 1
  model <- res$CALL$pars$model
  x <- res$CALL$pars$data_pre_fixed
  if(model=="MLVAR" && data_type != 1){
    if (data_type==2) x <- as.data.frame(x)
    if (data_type==3) x <- longi_to_dataframe(x)
  }
  if (model== "SMVAR"){
    x <- sapply(x, function(i) as.numeric(as.character(i)))
  }
  if(model=="GGVAR" && data_type != 3){
    x <- data_to_longi(x)
  }
  if(model=="GVAR" && (!data_type %in% c(1, 2))){
    x <- longi_to_matrix1(x)
  }
  if(model %in% c("SVARHL","SVARHLX","SVARHLMA") && (data_type != 3 | data_type != 2)){
    x <- as.matrix(x)
  }
  if(model %in% c("SVAR","SVECM") && data_type != 2){
    if (data_type==3) x <- longi_to_matrix1(x)
    if (data_type==1) x <- as.matrix(x)
  }
  return(x)
}