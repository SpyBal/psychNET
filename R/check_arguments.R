check_arguments <- function(res){
  model <- res$CALL$pars$model
  covariates <- res$CALL$pars$covariates
  lambda1 <- res$CALL$pars$lambda1
  lambda2 <- res$CALL$pars$lambda2
  penalty.type <- res$CALL$pars$penalty.type 
  optimality <- res$CALL$pars$optimality
  impute <- res$CALL$pars$impute
  additional_args <- res$CALL$pars$dots
  modls <- models_in()
  
  if (is.null(res$CALL$pars$lag) & !model %in% c("SVARHL", "SVARHLX",  "SVARHLMA")){
    message(red("\n Lag argument cannot be 'NULL'. A lag(1) model will be fitted instead..\n"))
    res$CALL$pars$lag <- 1
  }
  
  if (!(model %in% c(modls$ind,modls$pop))) stop(paste("Model should be one of the following: ",paste(c(modls$ind,modls$pop), collapse = ", ")))
  if ((!is.null(covariates)) & (! model %in% c("SVARHLX","VAR")) ) stop("Covariates are allowed only when model is: VAR or SVARHLX")
  
  doc_style <- magenta $ bold $ italic  
  cit_style <- green $ italic  
  
  if (model %in% c("VAR","MLVAR","DFM")){
    if (model == "VAR"){
      message(doc_style("\n See the documentation of the function ?VAR from the R package vars for details\n"))
      message(cit_style("Bernhard Pfaff (2008). VAR, SVAR and SVEC Models: Implementation Within R Package vars. Journal of
  Statistical Software 27(4). URL http://www.jstatsoft.org/v27/i04/.\nPfaff, B. (2008) Analysis of Integrated and Cointegrated Time Series with R. Second Edition. Springer, New
      York. ISBN 0-387-27960-1\n"))

      if (!is.null(optimality))  message(red("\ncriterion argument is not used for model='VAR'. If you want to specify additional arguments use the ... structure. See ?VAR.\n"))
      if (!is.null(penalty.type)) message(red("\nPenalty is not used for model='VAR'. If you want to specify additional arguments use the ... structure. See ?VAR.\n"))
      if ((!is.null(lambda1) | !is.null(lambda2) )) message(red("\nLambda parameters are not used for model='VAR'. If you want to specify additional arguments use the ... structure. See ?VAR.\n"))
      if (!is.null(res$CALL$pars$nFact)) message(red("\nNumber of factors are not used when model='VAR'. If you want to specify additional arguments use the ... structure. See ?VAR.\n"))
      res$CALL$pars$optimality <- res$CALL$pars$penalty.type <- res$CALL$pars$lambda1 <- res$CALL$pars$lambda1 <- res$CALL$pars$nFact <- NULL
    } 
    
    if (model == "MLVAR"){
      message(doc_style("See the documentation of the function ?mlVAR from the R package mlVAR for details\n "))
      message(cit_style("Sacha Epskamp, Marie K. Deserno and Laura F. Bringmann (2019). mlVAR: Multi-Level Vector Autoregression. R
  package version 0.4.2. https://CRAN.R-project.org/package=mlVAR\n"))


      if (!is.null(optimality))  message(red("\ncriterion argument is not used for model='MLVAR'.\n If you want to specify additional arguments use the ... structure. See ?mlVAR.\n"))
      if (!is.null(penalty.type)) message(red("\nPenalty is not used for model='MLVAR'. If you want to specify additional arguments use the ... structure. See ?mlVAR.\n"))
      if ((!is.null(lambda1) | !is.null(lambda2) )) message(red("\nLambda parameters are not used for model='MLVAR'. If you want to specify additional arguments use the ... structure. See ?mlVAR.\n"))
      if (!is.null(covariates)) message(red("\ncovariates are not used for model='MLVAR'. If you want to specify additional arguments use the ... structure. See ?mlVAR.\n"))
      if (!is.null(res$CALL$pars$nFact)) message(red("\nNumber of factors are not used when model='MLVAR'. If you want to specify additional arguments use the ... structure. See ?mlVAR.\n"))
      res$CALL$pars$optimality <- res$CALL$pars$penalty.type <- res$CALL$pars$lambda1 <- res$CALL$pars$lambda1 <- res$CALL$pars$nFact <- NULL
    }
    
      if (model == "DFM"){
      message(doc_style("See the documentation of the function dfm from the R package dynfactoR for details. The package is available on github.\n "))
      message(cit_style("Bagdziunas R (2019). _dynfactoR: Dynamic factor model estimation for
                    nowcasting_. R package version 0.1.\n"))
      if (!is.null(optimality))  message(red("\ncriterion argument is not used for model= 'DFM'.\n If you want to specify additional arguments use the ... structure."))
      if (!is.null(penalty.type)) message(red("\nPenalty is not used for model='DFM'.\n If you want to specify additional arguments use the ... structure."))
      if ((!is.null(lambda1) | !is.null(lambda2) )) message(red("\nLambda parameters are not used for model='MLVAR'.\n If you want to specify additional arguments use the ... structure."))
      if (!is.null(covariates)) message(red("\ncovariates are not used for model='DFM'.\n If you want to specify additional arguments use the ... structure. 
                                        \nSee ?dfm.\n"))
      res$CALL$pars$optimality <- res$CALL$pars$penalty.type <- res$CALL$pars$lambda1 <- res$CALL$pars$lambda1 <- res$CALL$pars$covariates <- NULL
      if (is.null(res$CALL$pars$nFact)) res$CALL$pars$nFact <- 2
    }
  }
  
  
  if (model %in% c("SVARHL","SVARHLX","SVARHLMA" )){
    message(doc_style("See the documentation of the functions ?sparseVAR ?sparseVARX, and ?sparseVARMA from the R package bigtime for details\n "))
    message(cit_style("Ines Wilms, David S. Matteson, Jacob Bien and Sumanta Basu (2017). bigtime: Sparse Estimation of Large Time
  Series Models. R package version 0.1.0. https://CRAN.R-project.org/package=bigtime\n"))
    
    if (is.null(optimality)){
      message(red("\ncriterion argument has not been specified. It has been set automatically to 'CV'. "))
      res$CALL$pars$optimality <- "CV"
    }else{
      if (optimality != "CV") message(red("\ncriterion argument can be only 'CV' when model= 'SVARHL', 'SVARHLX', or 'SVARHLMA'. It has been set automatically to 'CV'. "))
      res$CALL$pars$optimality <- "CV"
    } 
    
    if (is.null(penalty.type)){
      if (res$CALL$pars$lag > 1){
        message(red("\npenalty argument has not been specified. It has been set automatically to 'HLag'. Other option is 'LASSO' when model= 'SVARHL', 'SVARHLX', or 'SVARHLMA' "))
        res$CALL$pars$penalty.type <- "HLag"
      }else{
        message(red("\npenalty argument has not been specified. It has been set automatically to 'LASSO'. Other option is 'HLag' when model= 'SVARHL', 'SVARHLX', or 'SVARHLMA' "))
        res$CALL$pars$penalty.type <- "LASSO"
      }

    }else{
      if (! penalty.type %in% c("HLag","LASSO")){
        if (res$CALL$pars$lag > 1){
          message(red("\npenalty should be 'HLag' or 'LASSO' when model='SVARHL', 'SVARHLX', or 'SVARHLMA'. It has been set automatically to 'HLag'.\n "))
          res$CALL$pars$penalty.type <- "HLag"
        }else{
          message(red("\npenalty argument has not been specified. It has been set automatically to 'LASSO'. Other option is 'HLag' when model= 'SVARHL', 'SVARHLX', or 'SVARHLMA' "))
          res$CALL$pars$penalty.type <- "LASSO"
        }
      }
    } 
    
    if (!is.null(lambda1)){
      if (length(lambda1) < 2){
        message(red("\nFor model = 'SVARHL', 'SVARHLX', or 'SVARHLMA' : The regularization parameter lambda1 needs to be a vector of length>1 or NULL otherwise. It has been set automatically to NULL\n"))
        res$CALL$pars$lambda1 <- NULL
      } 
    }
    
    if (!is.null(lambda2)){
      message(red("\nThe regularization parameter lambda2 is not used when model='SVARHL', 'SVARHLX', or 'SVARHLMA'. It has been set automatically to NULL\n"))
      res$CALL$pars$lambda2 <- NULL
    }
    
  }
  
  if (model %in% c("SVAR","SVECM" )){
    message(doc_style("\nSee the documentation of the functions ?fitVAR and ?fitVECM from the R package sparsevar for details\n "))
    message(cit_style("Simone Vazzoler, Lorenzo Frattarolo and Monica Billio (2016). sparsevar: A Package for Sparse VAR/VECM
  Estimation. R package version 0.0.10. https://CRAN.R-project.org/package=sparsevar\n"))

    if (is.null(optimality)){
      message(red("\ncriterion argument has not been specified. It has been set automatically to 'CV'.\n "))
      res$CALL$pars$optimality <- "CV"
    }else{
      if (optimality != "CV") message(red("\ncriterion argument can be only 'CV' when model='SVAR' or 'SVECM'. It has been set automatically to 'CV'.\n "))
      res$CALL$pars$optimality <- "CV"
    } 
    if (is.null(penalty.type)){
      message(red("\npenalty argument has not been specified.\n It has been set automatically to 'ENET'. Other options are 'SCAD' and 'MCP' when model='SVAR' or 'SVECM' "))
      res$CALL$pars$penalty.type <- "ENET"
      
    }else{
      if (! penalty.type %in% c("ENET","SCAD","MCP")){
        message(red("\npenalty should be 'ENET', 'SCAD', or 'MCP' when model= 'SVAR' or 'SVECM'. It has been set automatically to 'ENET'.\n "))
        res$CALL$pars$penalty.type <- "ENET"
      }
      
    } 
    
    if (!is.null(lambda1)){
      message(red("\nThe regularization parameter lambda1 is not used when model= 'SVAR' or 'SVECM'. It has been set automatically to NULL\n"))
      res$CALL$pars$lambda1 <- NULL
    }
    
    if (!is.null(lambda2)){
      message(red("\nThe regularization parameter lambda2 is not used when model= 'SVAR' or 'SVECM'. It has been set automatically to NULL\n"))
      res$CALL$pars$lambda2 <- NULL
    }
    
  }
  if (model== "SMVAR"){
    message(doc_style("See the documentation of the function ?mvar from the R package mgm for details\n "))
    message(cit_style("Jonas M. B. Haslbeck, Lourens J. Waldorp (2016). mgm: Structure Estimation for Time-Varying Mixed Graphical
  Models in high-dimensional Data arXiv preprint:1510.06871v2 URL http://arxiv.org/abs/1510.06871v2.\n"))
    
    if (is.null(optimality)){
      message(red("\ncriterion argument has not been specified. It has been set automatically to 'CV'. "))
      res$CALL$pars$optimality <- "CV"
    }else{
      if (!optimality %in% c("CV","EBIC")) message(red("\ncriterion argument can be only 'CV', or 'EBIC' when model='SMVAR'. It has been set automatically to 'CV'. For 'BIC' set criterion='EBIC' and lambdaGam=0\n"))
      res$CALL$pars$optimality <- "CV"
      if (optimality == "EBIC"){
        res$CALL$pars$optimality <- "EBIC"
      } 
    } 
    if (is.null(penalty.type)){
      message(red("\npenalty argument has not been specified. It has been set automatically to 'ENET'. No Other options when model='SMVAR' \n"))
      res$CALL$pars$penalty.type <- "ENET"
    }else{
      if (!penalty.type =="ENET"){
        message(red("\npenalty should be 'ENET', when model= 'SMVAR'. It has been set automatically to 'ENET'.\n "))
        res$CALL$pars$penalty.type <- "ENET"
      }
    } 
    
    if (!is.null(lambda2)){
      message(red("\nThe regularization parameter lambda2 is not used when model='SMVAR' .\n It has been set automatically to NULL\n"))
      res$CALL$pars$lambda1 <- NULL
    }
    res$CALL$pars$dots <- additional_args
    
  }
  
  if (model== "GVAR"){
    message(doc_style("See the documentation of the function ?graphicalVAR from the R package graphicalVAR for details\n "))
    message(cit_style("Sacha Epskamp (2018). graphicalVAR: Graphical VAR for Experience Sampling Data. R package version 0.2.2.
  https://CRAN.R-project.org/package=graphicalVAR\n"))
    
    if (is.null(optimality)){
      message(red("\ncriterion argument has not been specified. It has been set automatically to 'EBIC'. "))
      res$CALL$pars$optimality <- "EBIC"
    }else{
      if (optimality != "EBIC") message(red("\ncriterion argument can be only 'EBIC' when model='GVAR'. It has been set automatically to 'EBIC'. For using standard 'BIC' add gamma=0" ))
      res$CALL$pars$optimality <- "EBIC"
    } 
    if (is.null(penalty.type)){
      message(red("\npenalty argument has not been specified. It has been set automatically to 'LASSO'. No Other options when model='GVAR' "))
      res$CALL$pars$penalty.type <- "LASSO"
    }else{
      if (! penalty.type =="LASSO"){
        message(red("\npenalty should be 'LASSO', when model= 'GVAR'. It has been set automatically to 'LASSO'. "))
        res$CALL$pars$penalty.type <- "LASSO"
      }
    } 
  }
  if (model== "GGVAR"){
    message(doc_style("See function ?sparse.tscgm from the R package SparseTSCGM for details "))
    message(cit_style("Fentaw Abegaz and Ernst Wit (2016). SparseTSCGM: Sparse Time Series Chain Graphical Models. R package
  version 2.5. https://CRAN.R-project.org/package=SparseTSCGM\n"))
    
    if (is.null(optimality)) optimality <- "NULL"
    
    if (res$CALL$pars$lag > 2){
      message(red("\nThe implementation in sparseTSCGM supports up to a second order VAR only."))
      message(red("\nA VAR(2) will be fitted automatically."))
      res$CALL$pars$lag <- 2
    } 
    if (is.null(lambda1) | is.null(lambda2)){
      if (optimality=="NULL" | is.null(optimality)){
        message(red("\ncriterion argument can be one of:\n 'BIC','EBIC','MBIC','GIC','AIC' when model='GGVAR'. It has been set automatically to 'EBIC'. "))
        optimality <- "BIC"
        res$CALL$pars$optimality <- "BIC"
        res$CALL$pars$optimality_pack <- "bic"
      }
    }
    if (!is.null(lambda1) & !is.null(lambda2)){
      if (length(lambda1)==1 & length(lambda2)==1){
        if (optimality != "NULL" | !is.null(optimality)){
          message(red("\ncriterion is not used for a single pair of lambdas. "))
          optimality <- "NULL"
          res$CALL$pars$optimality <- "NULL"
          res$CALL$pars$optimality_pack <- "NULL"
        }
      }
      if (length(lambda1)>1 & length(lambda2)>1){
        if (optimality == "NULL" | is.null(optimality)){
          message(red("\ncriterion argument can be one of: 'BIC','EBIC','MBIC','GIC','AIC' when model='GGVAR'. It has been set automatically to 'EBIC'. "))
          optimality <- "BIC"
          res$CALL$pars$optimality <- "BIC"
          res$CALL$pars$optimality_pack <- "bic"
        }
      }
    }

   
    
    if (! optimality %in% c("NULL","BIC","EBIC","MBIC","GIC","AIC")){
      message(red("\ncriterion argument can be one of: 'NULL','BIC','EBIC','MBIC','GIC','AIC' when model='GGVAR'. It has been set automatically to 'EBIC'. "))
      res$CALL$pars$optimality <- "BIC"
      res$CALL$pars$optimality_pack <- "bic"
    }else{
      if (is.null(optimality) | optimality == "NULL"){
        res$CALL$pars$optimality <- "NULL"
        res$CALL$pars$optimality_pack <- "NULL"
      }
      if (optimality == "BIC"){
        res$CALL$pars$optimality_pack <- "bic"
      }
      if (optimality == "EBIC"){
        res$CALL$pars$optimality_pack <- "bic_ext"
      }
      if (optimality == "MBIC"){
        res$CALL$pars$optimality_pack <- "bic_mod"
      }
      if (optimality == "AIC"){
        res$CALL$pars$optimality_pack <- "aic"
      }
      if (optimality == "GIC"){
        res$CALL$pars$optimality_pack <- "gic"
      }
    }
      
    if (is.null(penalty.type)){
      message(red("\npenalty argument has not been specified. It has been set automatically to 'LASSO'. Other option is 'SCAD' \n "))
      res$CALL$pars$penalty.type <- "LASSO"
    }else{
      if (! penalty.type%in% c("LASSO","SCAD")){
        message(red("\npenalty should be 'LASSO' or 'SCAD' when model= 'GGVAR'. It has been set automatically to 'LASSO'. "))
        res$CALL$pars$penalty.type <- "LASSO"
      }
    } 
  }
  return(res)
  }