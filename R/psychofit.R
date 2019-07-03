psychofit <- function(res,...){
  
  additional_args <- list(...)
  model <- res$CALL$pars$model
  no_sympt_index <- res$CALL$pars$no_sympt_dat_ind
  
  if (model=="VAR"){
    res$CALL$pars$dots <- additional_args
    res$fit <- vars::VAR(y=res$CALL$pars$data_pre_fixed[,-no_sympt_index],
                         p=res$CALL$pars$lag,
                         exogen = res$CALL$pars$covariates, ...)
    
    res$results <- list()
    res$results$A <- Acoef(res$fit)
    res$results$A <- lapply(res$results$A, function(x){
      dimnames(x) <- list(res$CALL$pars$symptoms,res$CALL$pars$symptoms)
      x
    })
    res$results$Covariance <- summary(res$fit)$covres
    res$results$residuals <- residuals(res$fit)
    res$results$Dir_net <- lapply(res$results$A,t)
    res$results$Dir_net <- lapply(res$results$Dir_net, function(x){
      dimnames(x) <- list(res$CALL$pars$symptoms,res$CALL$pars$symptoms)
      x
    })
    res$results$Pval_edge <- lapply(Pcoef(res$fit),t)
    return(res)
  }
  
  if (model=="SVAR"){
    res$CALL$pars$dots <- additional_args
    res$fit <- fitVAR(as.matrix(res$CALL$pars$data_pre_fixed[,-no_sympt_index]),
                      p=res$CALL$pars$lag,
                      penalty = res$CALL$pars$penalty.type,
                      method = tolower(res$CALL$pars$optimality),...)
    res$results <- list()
    res$results$lambda <- res$fit$lambda
    res$results$mu <- res$fit$mu
    res$results$A <- res$fit$A
    res$results$A <- lapply(res$results$A, function(x){
      dimnames(x) <- list(res$CALL$pars$symptoms,res$CALL$pars$symptoms)
      x
    })
    res$results$Covariance <- res$fit$sigma
    res$results$residuals <- res$fit$residuals
    res$results$lambda <- res$fit$lambda
    res$results$Dir_net <- lapply(res$results$A,t)
    res$results$Dir_net <- lapply(res$results$Dir_net, function(x){
      dimnames(x) <- list(res$CALL$pars$symptoms,res$CALL$pars$symptoms)
      x
    })
    return(res)
    
  }
  
  if (model=="SVECM"){
    res$CALL$pars$dots <- additional_args
    res$fit <- fitVECM(as.matrix(res$CALL$pars$data_pre_fixed[,-no_sympt_index]),
                       p=res$CALL$pars$lag,
                       penalty = res$CALL$pars$penalty.type,
                       method = tolower(res$CALL$pars$optimality),
                       logScale = FALSE,...)
    
    res$results <- list()
    res$results$A <- list()
    res$results$A <- res$fit$Pi
    if (!is.list(res$results$A))  res$results$A <- list(res$results$A)
    res$results$A <- lapply(res$results$A, function(x){
      dimnames(x) <- list(res$CALL$pars$symptoms,res$CALL$pars$symptoms)
      x
    })
    res$results$G <- res$fit$G
    res$results$G <- lapply(res$results$G, function(x){
      dimnames(x) <- list(res$CALL$pars$symptoms,res$CALL$pars$symptoms)
      x
    })
    res$results$Covariance <- res$fit$sigma
    res$results$residuals <- res$fit$residuals
    res$results$lambda <- res$fit$lambda
    res$results$Dir_net <- lapply(res$results$A,t)
    res$results$Dir_net <- lapply(res$results$Dir_net, function(x){
      dimnames(x) <- list(res$CALL$pars$symptoms,res$CALL$pars$symptoms)
      x
    })
    return(res)
    
  }
  
  if (model=="GVAR"){
    res$CALL$pars$dots <- additional_args
    if (is.null(res$CALL$pars$lambda1)){
      if (is.null(res$CALL$pars$lambda2)){
        res$fit <- graphicalVAR(res$CALL$pars$data_pre_fixed,
                                lags = 1:res$CALL$pars$lag,
                                vars = res$CALL$pars$symptoms,
                                scale = FALSE,
                                idvar = "ID",
                                beepvar ="BEEP",
                                dayvar = "DAY",
                                centerWithin = FALSE,
                                deleteMissings = TRUE,...)
        
      }else{
        res$fit <- graphicalVAR(res$CALL$pars$data_pre_fixed,
                                lags = 1:res$CALL$pars$lag,
                                verbose = TRUE, #res$CALL$pars$dots$verbose,
                                nLambda = res$CALL$pars$dots$nLambda,
                                lambda_beta = res$CALL$pars$lambda2,
                                vars = res$CALL$pars$symptoms,
                                scale = FALSE,
                                idvar = "ID",
                                beepvar ="BEEP",
                                dayvar = "DAY",
                                centerWithin = FALSE, 
                                deleteMissings = TRUE,...)
      }
      
    }else{
      if (is.null(res$CALL$pars$lambda2)){
        res$fit <- graphicalVAR(res$CALL$pars$data_pre_fixed,
                                lags = 1:res$CALL$pars$lag,
                                verbose = TRUE, #res$CALL$pars$dots$verbose,
                                nLambda = res$CALL$pars$dots$nLambda,
                                lambda_kappa = res$CALL$pars$lambda1,
                                vars = res$CALL$pars$symptoms,
                                scale = FALSE,
                                idvar = "ID",
                                beepvar ="BEEP",
                                dayvar = "DAY",
                                centerWithin = FALSE, 
                                deleteMissings = TRUE,...)
      }else{
        res$fit <- graphicalVAR(res$CALL$pars$data_pre_fixed,
                                lags = 1:res$CALL$pars$lag,
                                verbose = TRUE, #res$CALL$pars$dots$verbose,
                                nLambda = res$CALL$pars$dots$nLambda,
                                lambda_kappa = res$CALL$pars$lambda1,
                                lambda_beta = res$CALL$pars$lambda2,
                                vars = res$CALL$pars$symptoms,
                                gamma = additional_args$gamma,
                                scale = FALSE,
                                idvar = "ID",
                                beepvar ="BEEP",
                                dayvar = "DAY",
                                centerWithin = FALSE, 
                                deleteMissings = TRUE,...)
        
      }
      
    }
    res$results <- list()
    res$results$mu <- res$fit$beta[,1]
    lmats <- res$fit$beta[,-1]
    ind1 <- seq(1,ncol(lmats), by = res$CALL$pars$n_symptoms)
    ind2 <- c(ncol(lmats) / res$CALL$pars$lag,ncol(lmats))
    res$results$A <- lapply(1:res$CALL$pars$lag,function(i) lmats[,ind1[i]:ind2[i]])
    res$results$A <- lapply(res$results$A, function(x){
      dimnames(x) <- list(res$CALL$pars$symptoms,res$CALL$pars$symptoms)
      x
    })
    res$results$Covariance <- solve(res$fit$kappa)
    res$results$Dir_net <- lapply(res$results$A,t)
    res$results$Dir_net <- lapply(res$results$Dir_net, function(x){
      dimnames(x) <- list(res$CALL$pars$symptoms,res$CALL$pars$symptoms)
      x
    })
    res$results$UnDir_net <- res$fit$PCC
    return(res)
  } 
  

  if (model=="SMVAR"){
    if (is.null(res$CALL$pars$lambda1)){
      res$fit <- mvar_psych(data = res$data[,-no_sympt_index],
                      type = res$CALL$pars$dots$type[-no_sympt_index],
                      lambdaSel = res$CALL$pars$optimality,
                      level = res$CALL$pars$dots$level[-no_sympt_index],
                      lags = 1:res$CALL$pars$lag,
                      scale = FALSE,
                      beepvar = res$data[,"BEEP"],
                      dayvar = res$data[,"DAY"],...)
      
    }else{
      res$fit <- mvar_psych(data = res$data[,-no_sympt_index],
                      type = res$CALL$pars$dots$type[-no_sympt_index],
                      level = res$CALL$pars$dots$level[-no_sympt_index],
                      lambdaSel = res$CALL$pars$optimality,
                      lambdaSeq = res$CALL$pars$lambda1,
                      lags = 1:res$CALL$pars$lag,
                      scale = FALSE,
                      beepvar = res$data[,"BEEP"],
                      dayvar = res$data[,"DAY"],...)
      
    }
    res$results <- list()
    res$results$mu <- res$fit$intercepts
    res$results$A <- lapply(1:res$CALL$pars$lag, function(i) res$fit$wadj[,,i])
    res$results$A <- lapply(res$results$A, function(x){
      dimnames(x) <- list(res$CALL$pars$symptoms,res$CALL$pars$symptoms)
      x
    })
    res$results$Dir_net <- lapply(res$results$A,t)
    res$results$Dir_net <- lapply(res$results$Dir_net, function(x){
      dimnames(x) <- list(res$CALL$pars$symptoms,res$CALL$pars$symptoms)
      x
    })
    return(res)
    
  } 
  
  if (model== "SVARHL"){
    res$CALL$pars$penalty.type <- ifelse(res$CALL$pars$penalty.type=="LASSO","L1",res$CALL$pars$penalty.type)
    if (is.null(res$CALL$pars$lambda1)){
      res$fit <- sparseVAR(Y=as.matrix(sapply(res$CALL$pars$data_pre_fixed[,-no_sympt_index],as.numeric)),
                           VARpen = res$CALL$pars$penalty.type, 
                           p= res$CALL$pars$lag, ...)
      
    }else{
      res$fit <- sparseVAR(Y=as.matrix(sapply(res$CALL$pars$data_pre_fixed[,-no_sympt_index],as.numeric)),
                           VARpen = res$CALL$pars$penalty.type,
                           VARlseq = res$CALL$pars$lambda1, 
                           p= res$CALL$pars$lag, ...)
      
    }
    
    res$CALL$pars$lag <- res$fit$p
    res$results <- list()
    res$results$mu <- res$fit$phi0hat
    lmats <- res$fit$Phihat
    ind1 <- seq(1,ncol(lmats), by = res$CALL$pars$n_symptoms)
    if (length(ind1)==1){
      ind2 <- ncol(lmats)
    }else{
      ind2 <- rep(NA, length(ind1))
      test <- ind1 - 1
      ind2[1:(length(ind1)-1)] <- test[-1]
      ind2[length(ind1)] <- ncol(lmats)
    }
    res$results$A <- lapply(1:res$fit$p,function(i) lmats[,ind1[i]:ind2[i]])
    res$results$A <- lapply(res$results$A, function(x){
      dimnames(x) <- list(res$CALL$pars$symptoms,res$CALL$pars$symptoms)
      x
    })
    res$results$Dir_net <- lapply(res$results$A,t)
    res$results$Dir_net <- lapply(res$results$Dir_net, function(x){
      dimnames(x) <- list(res$CALL$pars$symptoms,res$CALL$pars$symptoms)
      x
    })
    return(res)
    
  }
  
  if (model== "SVARHLX"){
    res$CALL$pars$penalty.type <- ifelse(res$CALL$pars$penalty.type=="LASSO","L1",res$CALL$pars$penalty.type)
    if (is.null(res$CALL$pars$lambda1)){
      res$fit <- sparseVARX(Y=as.matrix(sapply(res$CALL$pars$data_pre_fixed[,-no_sympt_index],as.numeric)),
                            X= as.matrix(res$CALL$pars$covariates),
                            p=res$CALL$pars$lag,
                            VARXpen = res$CALL$pars$penalty.type,...)
      
    }else{
      res$fit <- sparseVARX(Y=as.matrix(sapply(res$CALL$pars$data_pre_fixed[,-no_sympt_index],as.numeric)),
                            X= as.matrix(res$CALL$pars$covariates),
                            p=res$CALL$pars$lag,
                            VARXpen = res$CALL$pars$penalty.type,
                            VARXlPhiseq = res$CALL$pars$lambda1,  ...)
      
    }
    res$CALL$pars$lag <- res$fit$p
    res$results <- list()
    res$results$mu <- res$fit$phi0hat
    lmats <- res$fit$Phihat
    ind1 <- seq(1,ncol(lmats), by = res$CALL$pars$n_symptoms)
    if (length(ind1)==1){
      ind2 <- ncol(lmats)
    }else{
      ind2 <- rep(NA, length(ind1))
      test <- ind1 - 1
      ind2[1:(length(ind1)-1)] <- test[-1]
      ind2[length(ind1)] <- ncol(lmats)
    }
    res$results$A <- lapply(1:res$fit$p,function(i) lmats[,ind1[i]:ind2[i]])
    lmats2 <- res$fit$Bhat
    ind11 <- seq(1,ncol(lmats2), by = ncol(res$CALL$pars$covariates))
    if (length(ind11)==1){
      ind22 <- ncol(lmats2)
    }else{
      ind22 <- rep(NA, length(ind11))
      test1 <- ind11 - 1
      ind22[1:(length(ind11)-1)] <- test1[-1]
      ind22[length(ind11)] <- ncol(lmats2)
    }
    res$results$B <- lapply(1:res$fit$s,function(i) lmats2[,ind11[i]:ind22[i]])
    res$results$A <- lapply(res$results$A, function(x){
      dimnames(x) <- list(res$CALL$pars$symptoms,res$CALL$pars$symptoms)
      x
    })
    res$results$B <- lapply(res$results$B, function(x){
      dimnames(x) <- list(res$CALL$pars$symptoms,colnames(res$CALL$pars$covariates))
      x
    })
    res$results$Dir_net <- lapply(res$results$A,t)
    res$results$Dir_net <- lapply(res$results$Dir_net, function(x){
      dimnames(x) <- list(res$CALL$pars$symptoms,res$CALL$pars$symptoms)
      x
    })
    return(res)
    
  }
  
  if (model== "SVARHLMA"){
    res$CALL$pars$penalty.type <- ifelse(res$CALL$pars$penalty.type=="LASSO","L1",res$CALL$pars$penalty.type)
    if (is.null(res$CALL$pars$lambda1)){
      res$fit <- sparseVARMA(Y=as.matrix(sapply(res$CALL$pars$data_pre_fixed[,-no_sympt_index],as.numeric)),
                             VARMAp = res$CALL$pars$lag,
                             VARMApen = res$CALL$pars$penalty.type,...)
      
    }else{
      res$fit <- sparseVARMA(Y=as.matrix(sapply(res$CALL$pars$data_pre_fixed[,-no_sympt_index],as.numeric)),
                             VARMAp = res$CALL$pars$lag,
                             VARMApen = res$CALL$pars$penalty.type,
                             VARMAlPhiseq = res$CALL$pars$lambda1,  ...)

    }
    
    res$CALL$pars$lag <- res$fit$VARMAp
    res$CALL$pars$lagMA <- res$fit$VARMAq
    res$results <- list()
    res$results$mu <- res$fit$VARphi0hat
    lmats <- res$fit$Phihat
    ind1 <- seq(1,ncol(lmats), by = res$CALL$pars$n_symptoms)
    if (length(ind1)==1){
      ind2 <- ncol(lmats)
    }else{
      ind2 <- rep(NA, length(ind1))
      test <- ind1 - 1
      ind2[1:(length(ind1)-1)] <- test[-1]
      ind2[length(ind1)] <- ncol(lmats)
    }
    res$results$A <- lapply(1:res$fit$VARMAp,function(i) lmats[,ind1[i]:ind2[i]])
    lmats2 <- res$fit$Thetahat
    ind11 <- seq(1,ncol(lmats2), by = res$CALL$pars$n_symptoms)
    if (length(ind11)==1){
      ind22 <- ncol(lmats2)
    }else{
      ind22 <- rep(NA, length(ind11))
      test1 <- ind11 - 1
      ind22[1:(length(ind11)-1)] <- test1[-1]
      ind22[length(ind11)] <- ncol(lmats2)
    }
    res$results$B <- lapply(1:res$fit$VARMAq,function(i) lmats2[,ind11[i]:ind22[i]])
    res$results$A <- lapply(res$results$A, function(x){
      dimnames(x) <- list(res$CALL$pars$symptoms,res$CALL$pars$symptoms)
      x
    })
    res$results$B <- lapply(res$results$B, function(x){
      dimnames(x) <- list(res$CALL$pars$symptoms,res$CALL$pars$symptoms)
      x
    })
    res$results$Dir_net <- lapply(res$results$A,t)
    res$results$Dir_net <- lapply(res$results$Dir_net, function(x){
      dimnames(x) <- list(res$CALL$pars$symptoms,res$CALL$pars$symptoms)
      x
    })
    return(res)
    
  }
  
  
  if (model== "DFM"){
    res$fit <- dfm(X=as.matrix(sapply(res$CALL$pars$data_pre_fixed[,-no_sympt_index],as.numeric)),
                   r=res$CALL$pars$nFact,
                   p=res$CALL$pars$lag,
                   rQ= "identity",
                   rC= "upper",
                   max_iter = 100000)
    res$results <- list()
    fact_names <- paste("factor",1:res$CALL$pars$nFact,sep = "")
    res$results$PCA_factors <- res$fit$pca
    colnames(res$results$PCA_factors) <- fact_names
    res$results$QML_factors <- res$fit$qml
    colnames(res$results$QML_factors) <- fact_names
    res$results$TS_factors <- res$fit$twostep
    colnames(res$results$TS_factors) <- fact_names

    lmats <- as.matrix(res$fit$A, ncol= res$CALL$pars$nFact*res$CALL$pars$lag, nrow=res$CALL$pars$nFact)
    ind1 <- seq(1,ncol(lmats), by = res$CALL$pars$nFact)
    if (length(ind1)==1){
      ind2 <- ncol(lmats)
    }else{
      ind2 <- rep(NA, length(ind1))
      test <- ind1 - 1
      ind2[1:(length(ind1)-1)] <- test[-1]
      ind2[length(ind1)] <- ncol(lmats)
    }

    res$results$A_fact <- lapply(1:res$CALL$pars$lag,function(i){
      mats <- as.matrix(lmats[,ind1[i]:ind2[i]])
      dimnames(mats) <- list(fact_names,fact_names)
      mats})
    res$results$Dir_net_fac <- lapply(res$results$A_fact,t)
    res$results$B_fac_symptoms <- as.matrix(res$fit$C)
    dimnames(res$results$B_fac_symptoms) <- list(res$CALL$pars$symptoms, fact_names)
    
    res$results$Obs_Covariance <- res$fit$R
    dimnames(res$results$Obs_Covariance) <- list(res$CALL$pars$symptoms, res$CALL$pars$symptoms)
    
    res$results$System_Covariance <- as.matrix(res$fit$Q)
    dimnames(res$results$System_Covariance) <- list(fact_names, fact_names)
    if (res$CALL$pars$lag==1){
      ftoV <-  factorToVAR(lambda = res$results$B_fac_symptoms, beta = res$results$A_fact[[1]], 
                           psi = res$results$System_Covariance,theta = res$results$Obs_Covariance )
      res$results$Dir_net <- list()
      res$results$Dir_net[[1]] <- ftoV$PDC
      res$results$UnDir_net <- ftoV$PDC
    }
    return(res)
  }
  
  
  if (model== "MLVAR"){
    res$fit <- mlVAR(data=res$CALL$pars$data_pre_fixed,
                     lags=1:res$CALL$pars$lag,
                     vars = res$CALL$pars$symptoms,
                     idvar = "ID",
                     dayvar = "DAY",
                     beepvar = "BEEP",
                     scaleWithin = FALSE,  ...)
    res$results <- list()
    res$results$mu <- res$fit$results$mu$mean
    res$results$mu_subject  <- res$fit$results$mu$subject
    
    res$results$A <- lapply(1:res$CALL$pars$lag, function(i) res$fit$results$Beta$mean[,,i])
    res$results$Covariance <- res$fit$results$Theta$cov$mean
    res$results$Dir_net <- lapply(1:res$CALL$pars$lag, function(i) psychoMLVAR(res$fit, type = "temporal", lag = i, partial = TRUE, SD=FALSE, order = res$CALL$pars$symptoms,
                                                        nonsig = "show"))
    
    res$results$UnDir_net <- psychoMLVAR(res$fit, type = "contemporaneous", partial = TRUE, SD=FALSE, order = res$CALL$pars$symptoms,
                      nonsig = "show")
    dimnames(res$results$UnDir_net) <- list(res$CALL$pars$symptoms,res$CALL$pars$symptoms)
    
    if (res$fit$input$temporal != "fixed"){
      res$results$Subjects <- lapply(1:res$CALL$pars$n_ID_unique, function(x) list())
      res$results$Subjects$Dir_net <- lapply(res$CALL$pars$ID_unique, function(i){
        lapply(1:res$CALL$pars$lag, function(j) psychoMLVAR(res$fit, type = "temporal", 
                                                            lag = j,
                                                            subject = i,
                                                            partial = TRUE, SD=FALSE, 
                                                            order = res$CALL$pars$symptoms,nonsig = "show"))
      })
      
      res$results$Subjects$UnDir_net <- lapply(res$CALL$pars$ID_unique, function(i){
        mat <- psychoMLVAR(res$fit, type = "contemporaneous", 
                           subject = i,
                           partial = TRUE, SD=FALSE, 
                           order = res$CALL$pars$symptoms,
                           nonsig = "show")
        dimnames(mat) <- list(res$CALL$pars$symptoms, res$CALL$pars$symptoms)
        mat
      })
      
    }
    

    return(res)
    
  }
  
  if (model== "GGVAR"){
    res$fit <- sparse.tscgm(res$data, 
                            lam1 = res$CALL$pars$lambda1,
                            lam2 = res$CALL$pars$lambda2,
                            model = ifelse(res$CALL$pars$lag==1,"ar1","ar2"),
                            penalty = tolower(res$CALL$pars$penalty.type),
                            optimality = res$CALL$pars$optimality_pack,...)
    res$results <- list()
    res$results$A  <- lapply(splitmatrix(res$fit$gamma, res$CALL$pars$lag),function(i) t(matrix(i, ncol = res$CALL$pars$n_symptoms,dimnames = list(res$CALL$pars$symptoms,res$CALL$pars$symptoms))))
    res$results$A <- t(res$fit$gamma)
    res$results$Covariance <- solve(res$fit$theta)
    res$results$Dir_net <-  lapply(splitmatrix(res$fit$gamma, res$CALL$pars$lag),function(i) matrix(i, ncol = res$CALL$pars$n_symptoms,dimnames = list(res$CALL$pars$symptoms,res$CALL$pars$symptoms)))
    res$results$UnDir_net <- as.matrix(wi2net(res$fit$theta))
    return(res)
    
  }
  
}

