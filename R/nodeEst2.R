nodeEst2 <- function (y, X, fam, lambdaSeq, lambdaSel, lambdaFolds, lambdaGam, 
          alphac, weights, n, nadj, v, type, level, emp_lev, overparameterize, 
          thresholdCat, obj, Xdat){

  if (missing(alphac)) alphac <- 1
  if (missing(lambdaSeq)) lambdaSeq <- NULL
  
  if (type[v] == "c"){
    fam = "multinomial"
    intercept <- thresholdCat
  } 
  if (type[v] == "g"){
    intercept <- TRUE
    fam = "gaussian"
  } 
  if (type[v] == "p"){
    intercept <- TRUE
    fam = "poisson"
  } 
  if (type[v] == "o") {
    fam = "cumulative"
  }
  
  if (lambdaSel == "EBIC") {
    
    if (fam == "cumulative"){
      nlag <- length(obj$call$lags)
      psi <-  factor(y, levels = sort(unique(y)))
      psi <- ordered(psi)
      Xd <- Xdat[complete.cases(Xdat),]
      level <- sapply(Xd,function(x) length(unique(x)))
      type <- rep(obj$call$type,3)
      
      newX <- matrix(NA, nrow = nrow(Xd), ncol = sum(level), dimnames = list(1:nrow(Xd), as.character(1:sum(level))) )
      for (j in 1:ncol(Xd)){
        colindstart <- which(is.na(newX[1,]))[1]
        colindend <- colindstart + level[j] - 1
        namesj <- paste(names(Xd)[j],unique(Xd[,j]), sep = "")
        if (type[j] %in% c("o","c")) newX[,colindstart:colindend] <-  as.matrix(dummy_cols(Xd[,j])[,-1])
        if (type[j] %in% c("g","p")) newX[,colindstart:colindend] <-  Xd[,j]
        colnames(newX)[colindstart:colindend] <- namesj
        
      }
      
      fit <- ordinalNet(x = newX, y = psi, standardize = FALSE, alpha = alphac, family = fam,lambdaVals = NULL,
                        threshOut = 1e-05,threshIn = 1e-05, maxiterOut = 2000, maxiterIn = 2000)
    
      fit$nulldev <- deviance(polr(psi~1))
      n_lambdas <- length(fit$lambdaVals)
      
      LL_null <- calcLL2(X = newX, y = psi, fit = fit, type = type, level = level, v = v, 
                         weights = rep(1:ncol(newX)), lambda = fit$lambdaVals[1], LLtype = "nullmodel")
      LL_sat <- 1/2 *fit$nulldev + LL_null
      deviance <- (1 - fit$devPct) * fit$nulldev
      LL_lambda_models <- -1/2 * deviance + LL_sat
      n_neighbors <- rep(NA, n_lambdas)
      for (i in 1:n_lambdas) n_neighbors[i] <- calcNeighbors2(fit = fit, lambda = fit$lambdaVals[i], type = type, level = level, v = v)
      EBIC_lambda <- -2 * LL_lambda_models + n_neighbors * log(nadj) + 2 * lambdaGam * n_neighbors * log(ncol(X))
      EBIC_min <- min(EBIC_lambda)
      ind_lambda_min <- which.min(EBIC_lambda)
      lambda_min <- fit$lambdaVals[ind_lambda_min]
      lambad_min_model <- coef(fit, whichLambda = ind_lambda_min)
      outlist <- list(EBIC = EBIC_min, deviance = deviance[which.min(EBIC_lambda)], 
                      lambda = lambda_min, alpha = alphac, model = lambad_min_model, 
                      fitobj = fit)
    }else{
      fit <- glmnet(x = X, y = y, family = fam, alpha = alphac, 
                    weights = weights, lambda = lambdaSeq, intercept = intercept)
      n_lambdas <- length(fit$lambda)
      LL_null <- calcLL(X = X, y = y, fit = fit, type = type, 
                        level = level, v = v, weights = weights, lambda = fit$lambda[1], 
                        LLtype = "nullmodel")
      LL_sat <- 1/2 * fit$nulldev + LL_null
      deviance <- (1 - fit$dev.ratio) * fit$nulldev
      LL_lambda_models <- -1/2 * deviance + LL_sat
      n_neighbors <- rep(NA, n_lambdas)
      for (i in 1:n_lambdas) n_neighbors[i] <- calcNeighbors(fit = fit, lambda = fit$lambda[i], type = type, level = level, v = v)
      EBIC_lambda <- -2 * LL_lambda_models + n_neighbors * log(nadj) + 2 * lambdaGam * n_neighbors * log(ncol(X))
      EBIC_min <- min(EBIC_lambda)
      ind_lambda_min <- which.min(EBIC_lambda)
      lambda_min <- fit$lambda[ind_lambda_min]
      lambad_min_model <- coef(fit, s = lambda_min)
      outlist <- list(EBIC = EBIC_min, deviance = deviance[which.min(EBIC_lambda)], 
                      lambda = lambda_min, alpha = alphac, model = lambad_min_model, 
                      fitobj = fit)
      
    }
    return(outlist)
  }
  if (lambdaSel == "CV") {
    if (fam=="cumulative"){
      nlag <- length(obj$call$lags)
      psi <-  factor(y, levels = sort(unique(y)))
      psi <- ordered(psi)
      Xd <- Xdat[complete.cases(Xdat),]
      level <- sapply(Xd,function(x) length(unique(x)))
      type <- rep(obj$call$type,3)
      
      newX <- matrix(NA, nrow = nrow(Xd), ncol = sum(level), dimnames = list(1:nrow(Xd), as.character(1:sum(level))) )
      for (j in 1:ncol(Xd)){
        colindstart <- which(is.na(newX[1,]))[1]
        colindend <- colindstart + level[j] - 1
        namesj <- paste(names(Xd)[j],unique(Xd[,j]), sep = "")
        if (type[j] %in% c("o","c")) newX[,colindstart:colindend] <-  as.matrix(dummy_cols(Xd[,j])[,-1])
        if (type[j] %in% c("g","p")) newX[,colindstart:colindend] <-  Xd[,j]
        colnames(newX)[colindstart:colindend] <- namesj
        
      }
      fit <- ordinalNetCV(x=newX,y=psi, family= fam, link="logit",alpha = alphac, standardize = FALSE,
                          threshOut = 1e-05,threshIn = 1e-05, maxiterOut = 2000, maxiterIn = 2000,
                          nFolds = lambdaFolds,tuneMethod = "cvMisclass", printProgress = FALSE)
      lambda_min <- fit$lambdaVals[round(median(fit$bestLambdaIndex))]
      lambad_min_model <- coef(fit$fit, whichLambda =round(median(fit$bestLambdaIndex)))
      LL_model <- calcLL2(X = newX, y = psi, fit = fit$fit, type = type, level = level, v = v, 
                          weights = rep(1:ncol(newX)), lambda = lambda_min, LLtype = "model")
      n_neighbors <- calcNeighbors2(fit = fit$fit, lambda = lambda_min, type = type, level = level, v = v)
      EBIC <- -2 * LL_model + n_neighbors * log(nadj) + 2 * lambdaGam * log(ncol(newX))
      outlist <- list(EBIC = EBIC, deviance = NULL, lambda = lambda_min, 
                      alpha = alphac, model = lambad_min_model, fitobj = fit$fit)
    }else{
      fit <- cv.glmnet(x = X, y = y, family = fam, alpha = alphac, 
                       weights = weights, nfolds = lambdaFolds, type.measure = "deviance", 
                       lambda = lambdaSeq, intercept = intercept)
      lambda_min <- fit$lambda.min
      lambad_min_model <- coef(fit, s = lambda_min)
      LL_model <- calcLL(X = X, y = y, fit = fit, type = type, 
                         level = level, v = v, weights = weights, lambda = lambda_min, 
                         LLtype = "model")
      n_neighbors <- calcNeighbors(fit = fit, lambda = lambda_min, 
                                   type = type, level = level, v = v)
      EBIC <- -2 * LL_model + n_neighbors * log(nadj) + 2 * lambdaGam * log(ncol(X))
      outlist <- list(EBIC = EBIC, deviance = NULL, lambda = lambda_min, 
                      alpha = alphac, model = lambad_min_model, fitobj = fit)
    }
    return(outlist)
  }
}
