mvar_psych <- function (data, type, level, lambdaSeq, lambdaSel, lambdaFolds, 
          lambdaGam, alphaSeq, alphaSel, alphaFolds, alphaGam, lags, 
          consec, beepvar, dayvar, weights, threshold, method, binarySign, 
          scale, verbatim, pbar, warnings, saveModels, saveData, overparameterize, 
          thresholdCat, signInfo, ...) 
{
  n <- nrow(data)
  p <- ncol(data)
  data <- as.matrix(data)
  n_var <- n - max(lags)
  n_lags <- length(lags)
  max_lags <- max(lags)
  args <- list(...)
  if (missing(lambdaSeq)) lambdaSeq <- NULL
  if (missing(lambdaSel)) lambdaSel <- "CV"
  if (missing(lambdaFolds)) lambdaFolds <- 10
  if (missing(lambdaGam)) lambdaGam <- 0.25
  if (missing(alphaSeq)) alphaSeq <- 1
  if (missing(alphaSel)) alphaSel <- "CV"
  if (missing(alphaFolds)) alphaFolds <- 10
  if (missing(alphaGam)) alphaGam <- 0.25
  if (missing(lags)) lags <- 1
  if (missing(consec)) consec <- NULL
  if (missing(beepvar)) beepvar <- NULL
  if (missing(dayvar)) dayvar <- NULL
  if (missing(weights)) weights <- rep(1, n)
  if (missing(threshold)) threshold <- "LW"
  if (missing(method)) method <- "glm"
  if (missing(binarySign)) {
    if (!is.null(args$binary.sign)) binarySign <- args$binary.sign else binarySign <- FALSE
  }
  if (!is.null(args$binary.sign)) {
    warning("The argument 'binary.sign' is deprecated Use 'binarySign' instead.")
  }
  if (missing(verbatim)) verbatim <- FALSE
  if (missing(pbar)) pbar <- TRUE
  if (missing(warnings)) warnings <- TRUE
  if (missing(saveModels)) saveModels <- TRUE
  if (missing(saveData)) saveData <- FALSE
  if (missing(overparameterize)) overparameterize <- FALSE
  if (missing(thresholdCat)) thresholdCat <- ifelse(overparameterize, TRUE, TRUE)
  if (missing(signInfo)) signInfo <- TRUE
  if (missing(scale)) scale <- TRUE
  if (verbatim) pbar <- FALSE
  if (verbatim) warnings <- FALSE
  weights_initial <- weights
  if (!is.null(consec) & !is.null(beepvar)) stop("Please specify the consecutiveness of measurements either via consec, or via dayvar and beepvar")
  if (!is.null(consec) & !is.null(dayvar)) stop("Please specify the consecutiveness of measurements either via consec, or via dayvar and beepvar")
  if (!is.null(dayvar)) if (is.null(beepvar)) stop("Argument beepvar not specified.")
  if (!is.null(beepvar)) if (is.null(dayvar)) stop("Argument dayvar not specified.")
  if (!is.null(beepvar) & !is.null(dayvar)) consec <- beepday2consec(beepvar = beepvar, dayvar = dayvar)
  emp_lev <- rep(NA, p)
  type[sapply(1:length(emp_lev), function(i) type[i] == "o" & level[i] <=2)] <- "c"
  ind_cat <- which(type == "c")
  ord_cat <- which(type == "o")
  
  if (length(ind_cat) > 0) for (i in 1:length(ind_cat)) emp_lev[ind_cat][i] <- length(unique(data[, ind_cat[i]]))
  if (length(ord_cat) > 0) for (i in 1:length(ord_cat)) emp_lev[ord_cat][i] <- length(unique(data[, ord_cat[i]]))
  emp_lev[which(!(type %in% c("c","o")))] <- 1

  if (!missing(level)) {
    level_check <- level != emp_lev
    if (sum(level_check) > 0) 
      stop(paste0("Provided levels not equal to levels in data for variables ", paste((1:p)[level_check], collapse = ", ")))
  }
  level <- emp_lev
  if (!missing(weights)) weights <- weights/max(weights)
  nadj <- sum(weights)
  ind_Gauss <- which(type == "g")
  if (scale) for (i in ind_Gauss) data[, i] <- scale(data[, i])
  if (nrow(data) < 2) ("The data matrix has to have at least 2 rows.")
  if (missing(data)) stop("No data provided")
  if (ncol(data) < 2) stop("At least 2 variables required")
  if (any(is.na(data))) stop("No missing values permitted.")
  if (any(!is.finite(as.matrix(data)))) stop("No infinite values permitted.")
  if (any(!(apply(data, 2, class) %in% c("numeric", "integer")))) stop("Only integer and numeric values permitted.")
  if (!is.null(consec)) 
    if (length(consec) != nrow(data)) 
      stop("The length of consec has to be equal to the number of rows of the data matrix.")
  if (!(threshold %in% c("none", "LW", "HW"))) stop("Please select one of the three threshold options \"HW\", \"LW\" and \"none\" ")
  if (is.null(lags)) stop("No lags specified")
  if (any(duplicated(lags))) stop("No duplicates allowed in specified lags.")
  if (any(lags < 1)) stop("Specified lags have to be in {1, 2, ..., n -1}")
  if (any(round(lags) != lags)) stop("Specified lags have to be positive integers")
  if (missing(type)) stop("No type vector provided.")
  if (sum(!(type %in% c("g", "c", "p","o"))) > 0) stop("Only Gaussian 'g', Poisson 'p' or categorical 'c' variables permitted.")
  if (ncol(data) != length(type)) stop("Number of variables is not equal to length of type vector.")
  if (!missing(level)) 
    if (ncol(data) != length(level)) stop("Number of variables is not equal to length of level vector.")
  if ((nrow(data)) != length(weights)) stop("Weights vector has to be equal to the number of observations")
  if ("p" %in% type) {
    ind_Pois <- which(type == "p")
    nPois <- length(ind_Pois)
    v_PoisCheck <- rep(NA, length = nPois)
    for (i in 1:nPois) v_PoisCheck[i] <- sum(data[, ind_Pois[i]] != round(data[, ind_Pois[i]])) > 0
    if (sum(v_PoisCheck) > 0) stop("Only integers permitted for Poisson variables.")
  }
  ind_cat <- which(type == "c")
  ind_binary <- rep(NA, length(ind_cat))
  ind_binary <- as.logical(ind_binary)
  if (length(ind_cat) > 0) {
    for (i in 1:length(ind_cat)) ind_binary[i] <- length(unique(data[, ind_cat[i]])) == 2
  }
  if (sum(ind_binary) > 0) {
    check_binary <- rep(NA, sum(ind_binary))
    for (i in 1:sum(ind_binary)) check_binary[i] <- sum(!(unique(data[, ind_cat[ind_binary][i]]) %in% c(0, 1)))
  }
  if (binarySign) {
    if (sum(check_binary) > 0) 
      stop(paste0("If binarySign = TRUE, all binary variables have to be coded {0,1}. Not satisfied in variable(s) ", 
                  paste(ind_cat[ind_binary][check_binary > 0], collapse = ", ")))
  }
  colnames(data)[1:p] <- paste("V", 1:p, ".", sep = "")
  data <- as.data.frame(data)
  data_lagged <- lagData(data = data, lags = lags, consec = consec)
  data_response <- data_lagged$data_response
  l_data_lags <- data_lagged$l_data_lags
  n_design <- nrow(data_response)
  data_response <- data_response[data_lagged$included, ]
  l_data_lags <- lapply(l_data_lags, function(x) x[data_lagged$included, ])
  weights_design <- weights[data_lagged$included]
  nadj <- sum(weights_design)
  if (!is.null(args$bootstrap)) {
    if (args$bootstrap) {
      data_lagged <- lagData(data = data, lags = lags, consec = consec)
      data_response <- data_lagged$data_response
      l_data_lags <- lapply(data_lagged$l_data_lags, function(x) x)
      data_response <- data_response[args$boot_ind, ]
      l_data_lags <- lapply(l_data_lags, function(x) x[args$boot_ind, ])
      weights_design <- weights[args$boot_ind]
    }
  }
  n_lags <- length(lags)
  for (lag in 1:n_lags) {
    glmnetRequirements_psycho(data = l_data_lags[[lag]], type = type, weights = weights_design)
    
  }
  mvarobj <- list(call = NULL, wadj = NULL, signs = NULL, edgecolor = NULL, 
                  rawlags = list(), intercepts = NULL, nodemodels = NULL)
  mvarobj$call <- list(data = NULL, data_lagged = NULL, type = type, 
                       level = level, lambdaSeq = lambdaSeq, lambdaSel = lambdaSel, 
                       lambdaFolds = lambdaFolds, lambdaGam = lambdaGam, alphaSeq = alphaSeq, 
                       alphaSel = alphaSel, alphaFolds = alphaFolds, alphaGam = alphaGam, 
                       lags = lags, consec = consec, beepvar = beepvar, dayvar = dayvar, 
                       weights = weights_initial, weights_design = weights_design, 
                       threshold = threshold, method = method, binarySign = binarySign, 
                       scale = scale, verbatim = verbatim, pbar = pbar, warnings = warnings, 
                       saveModels = saveModels, saveData = saveData, overparameterize = overparameterize, 
                       thresholdCat = thresholdCat, signInfo = signInfo)
  if (saveData) mvarobj$call$data <- data
  mvarobj$call$data_lagged <- data_lagged
  if (!saveData) {
    mvarobj$call$data_lagged$l_data_lags <- NULL
    mvarobj$call$data_lagged$data_response <- NULL
  }
  mvarobj$call$consec <- consec
  if (pbar == TRUE) 
    pb <- txtProgressBar(min = 0, max = p, initial = 0, char = "-", style = 3)
  npar_standard <- rep(NA, p)
  for (v in 1:p) {
    y <- data_response[, v]
    data_input_MM <- do.call(cbind, l_data_lags)
    data_input_MM <- as.data.frame(data_input_MM)
    for (i in which(type == "c")) data_input_MM[, i] <- factor(data_input_MM[, i], levels = sort(unique(data_input_MM[, i])))
    data_v <- data.frame(y, data_input_MM)
    form <- as.formula("y ~ (.)")
    X_standard <- model.matrix( form , data = data_v)[, -1]
    npar_standard[v] <- ncol(X_standard)
    if (overparameterize){
      type_aug <- rep(type, n_lags)
      level_aug <- rep(level, n_lags)
      X_over <- ModelMatrix(data = data_input_MM, type = type_aug, level = level_aug, labels = colnames(data_input_MM), d = 1, v = NULL)
      X <- X_over
    }else{
      X <- X_standard
    }
    n_alpha <- length(alphaSeq)
    if (alphaSel == "CV") {
      l_alphaModels <- list()
      ind <- sample(1:alphaFolds, size = n_design, replace = TRUE)
      v_mean_OOS_deviance <- rep(NA, n_alpha)
      if (n_alpha > 1) {
        for (a in 1:n_alpha) {
          l_foldmodels <- list()
          v_OOS_deviance <- rep(NA, alphaFolds)
          for (fold in 1:alphaFolds) {
            train_X <- X[ind != fold, ]
            train_y <- y[ind != fold]
            test_X <- X[ind == fold, ]
            test_y <- y[ind == fold]
            n_train <- nrow(train_X)
            nadj_train <- sum(weights_design[ind != fold])
            l_foldmodels[[fold]] <- nodeEst2(y = train_y, 
                                            X = train_X, lambdaSeq = lambdaSeq, lambdaSel = lambdaSel, 
                                            lambdaFolds = lambdaFolds, lambdaGam = lambdaGam, 
                                            alphac = alphaSeq[a], weights = weights_design[ind != fold], n = n_train, nadj = nadj_train, 
                                            v = v, type = type, level = level, emp_lev = emp_lev, 
                                            overparameterize = overparameterize, thresholdCat = thresholdCat, obj=mvarobj, Xdat=data_input_MM)
            LL_model <- calcLL2(X = test_X, y = test_y, 
                                fit = l_foldmodels[[fold]]$fitobj, type = type, 
                                level = level, v = v, weights = weights_design[ind == 
                                                                                 fold], lambda = l_foldmodels[[fold]]$lambda, 
                                LLtype = "model")
            LL_saturated <- calcLL2(X = test_X, y = test_y, 
                                    fit = l_foldmodels[[fold]]$fitobj, type = type, 
                                    level = level, v = v, weights = weights_design[ind == 
                                                                                     fold], lambda = l_foldmodels[[fold]]$lambda, 
                                    LLtype = "saturated")
            v_OOS_deviance[fold] <- 2 * (LL_saturated - LL_model)
          }
          v_mean_OOS_deviance[a] <- mean(v_OOS_deviance)
        }
        alpha_select <- alphaSeq[which.min(v_mean_OOS_deviance)]
      }
      else {
        alpha_select <- alphaSeq
      }
      model <- nodeEst2(y = y, X = X, lambdaSeq = lambdaSeq, 
                       lambdaSel = lambdaSel, lambdaFolds = lambdaFolds, 
                       lambdaGam = lambdaGam, alphac = alpha_select, 
                       weights = weights_design, n = n_design, nadj = nadj, 
                       v = v, type = type, level = level, emp_lev = emp_lev, 
                       overparameterize = overparameterize, thresholdCat = thresholdCat, obj=mvarobj, Xdat=data_input_MM)
      mvarobj$nodemodels[[v]] <- model
    }
    if (alphaSel == "EBIC") {
      l_alphaModels <- list()
      EBIC_Seq <- rep(NA, n_alpha)
      for (a in 1:n_alpha) {
        l_alphaModels[[a]] <- nodeEst2(y = y, X = X, lambdaSeq = lambdaSeq, lambdaSel = lambdaSel, 
                                       lambdaFolds = lambdaFolds, lambdaGam = lambdaGam, 
                                       alphac = alphaSeq[a], weights = weights_design, n = n, nadj = nadj, 
                                       v = v, type = type, level = level, 
                                       emp_lev = emp_lev, overparameterize = overparameterize, 
                                       thresholdCat = thresholdCat, obj=mvarobj, Xdat=data_input_MM)
        EBIC_Seq[a] <- l_alphaModels[[a]]$EBIC
      }
      ind_minEBIC_model <- which.min(EBIC_Seq)
      mvarobj$nodemodels[[v]] <- l_alphaModels[[ind_minEBIC_model]]
    }
    if (pbar == TRUE) setTxtProgressBar(pb, v)
  }
  no_lags <- length(lags)
  l_Par <- vector("list", length = p)
  l_preds_dummy <- vector("list", p * no_lags)
  l_Par <- lapply(l_Par, function(x) l_preds_dummy)
  l_intercepts <- vector("list", length = p)
  pred_names <- unlist(lapply(l_data_lags, colnames))
  if (length(pred_names) != (no_lags * p)) stop("Length of predictor variable names does not match calculated number of predictor varibales")
  for (v in 1:p) {
    for (v2 in 1:(p * no_lags)) {
      if (type[v] == "c") {
        n_cats <- level[v]
        for (cat in 1:n_cats) {
          par_part <- mvarobj$nodemodels[[v]]$model[[cat]]
          par_part_ni <- par_part[-1]
          l_intercepts[[v]][[cat]] <- par_part[1]
          d <- 1
          if (threshold == "LW") tau <- sqrt(d) * sqrt(sum(par_part_ni^2)) * sqrt(log(npar_standard[v])/n)
          if (threshold == "HW") tau <- d * sqrt(log(npar_standard[v])/n)
          if (threshold == "none") tau <- 0
          par_part[abs(par_part) < tau] <- 0
          mvarobj$nodemodels[[v]]$tau <- tau
          ind_v2 <- grepl(pred_names[v2], rownames(par_part)[-1], fixed = TRUE)
          l_Par[[v]][[v2]][[cat]] <- par_part[-1, ][ind_v2]
        }
      }else {
        if (type[v] == "o"){
          n_cats <- level[v]
          for (cat in 1:n_cats) {
            par_part <- mvarobj$nodemodels[[v]]$model
            par_part_ni <- par_part[- seq(1,(n_cats-1))]
            l_intercepts[[v]] <- par_part[seq(1,(n_cats-1))]
            d <- 1
            if (threshold == "LW") tau <- sqrt(d) * sqrt(sum(par_part_ni^2)) * sqrt(log(npar_standard[v])/n)
            if (threshold == "HW") tau <- d * sqrt(log(npar_standard[v])/n)
            if (threshold == "none") tau <- 0
            par_part[abs(par_part) < tau] <- 0
            mvarobj$nodemodels[[v]]$tau <- tau
            ind_v2 <- grepl(pred_names[v2], names(par_part)[- seq(1,(n_cats-1))], fixed = TRUE)
            l_Par[[v]][[v2]] <- par_part[-seq(1,(n_cats-1))][ind_v2]
          }
        }else{
          par_part <- mvarobj$nodemodels[[v]]$model
          par_part_ni <- par_part[-1]
          l_intercepts[[v]] <- par_part[1]
          d <- 1
          if (threshold == "LW") tau <- sqrt(d) * sqrt(sum(par_part_ni^2)) * sqrt(log(npar_standard[v])/n_var)
          if (threshold == "HW") tau <- d * sqrt(log(npar_standard[v])/n_var)
          if (threshold == "none") tau <- 0
          par_part[abs(par_part) < tau] <- 0
          mvarobj$nodemodels[[v]]$tau <- tau
          ind_v2 <- grepl(pred_names[v2], rownames(par_part)[-1], fixed = TRUE)
          l_Par[[v]][[v2]] <- par_part[-1, ][ind_v2]
        }
      }
    }
  }
  l_Par_red <- vector("list", length = p)
  l_signs <- vector("list", length = p)
  l_preds_dummy <- vector("list", p * no_lags)
  l_Par_red <- lapply(l_Par_red, function(x) l_preds_dummy)
  l_signs <- lapply(l_signs, function(x) l_preds_dummy)
  type_long <- rep(type, times = no_lags)
  ind_binary_long <- rep(1:p %in% ind_cat[ind_binary], times = no_lags)
  ind_cat_long <- rep(ind_cat, times = no_lags)
  if (binarySign) {
    set_signdefined <- c(which(type_long == "p"), which(type_long == "g"), ind_cat_long[ind_binary_long])
  }else {
    set_signdefined <- c(which(type_long == "p"), which(type_long == "g"))
  }
  for (v in 1:p) {
    for (v2 in 1:(p * no_lags)) {
      if (type[v] == "o") l_Par_red[[v]][[v2]] <- abs(mean(unlist(l_Par[[v]][[v2]]))) else  l_Par_red[[v]][[v2]] <- mean(abs(unlist(l_Par[[v]][[v2]])))
      if (l_Par_red[[v]][[v2]] != 0) {
        if (sum(!(c(v, v2) %in% set_signdefined)) == 0) {
          if (overparameterize) {
            if (type[v] == "c") {
              if (type[v2] == "c") {
                if (l_Par[[v]][[v2]][[1]][1] != 0) sign_sel <- sign(l_Par[[v]][[v2]][[1]][1])
                if (l_Par[[v]][[v2]][[1]][2] != 0) sign_sel <- -sign(l_Par[[v]][[v2]][[1]][2])
              }else {
                sign_sel <- sign(l_Par[[v]][[v2]][2])
              }
            }else {
              sign_sel <- sign(l_Par[[v]][[v2]])
            }
          }else {
            if (type[v] == "c") {
              sign_sel <- sign(l_Par[[v]][[v2]][[2]])
            }else {
              sign_sel <- sign(l_Par[[v]][[v2]])
            }
          }
          l_signs[[v]][[v2]] <- sign_sel
        }else {
          l_signs[[v]][[v2]] <- 0
        }
      }else {
        l_signs[[v]][[v2]] <- NA
      }
    }
  }
  n_lags <- length(lags)
  lag_indicator <- rep(1:n_lags, each = p)
  list_wadj <- list_signs <- list()
  a_wadj <- a_signs <- array(dim = c(p, p, n_lags))
  a_edgecolor <- array("darkgrey", dim = c(p, p, n_lags))
  mvarobj$rawlags <- vector("list", length = n_lags)
  for (lag in 1:n_lags) {
    lag_seq <- (1:(p * n_lags))[lag_indicator == lag]
    l_lag <- vector("list", length = p)
    l_dum <- vector("list", length = length(lag_seq))
    l_lag <- lapply(l_lag, function(x) l_dum)
    wadj <- m_signs <- matrix(NA, p, p)
    for (i in 1:p) {
      k <- 1
      for (j in lag_seq) {
        a_wadj[i, k, lag] <- l_Par_red[[i]][[j]]
        l_lag[[i]][[j]] <- l_Par[[i]][[j]]
        a_signs[i, k, lag] <- l_signs[[i]][[j]]
        if (!is.na(a_signs[i, k, lag])) {
          if (a_signs[i, k, lag] == -1) 
            a_edgecolor[i, k, lag] <- "red"
          if (a_signs[i, k, lag] == 1) 
            a_edgecolor[i, k, lag] <- "darkgreen"
        }
        k <- k + 1
      }
    }
    mvarobj$wadj <- a_wadj
    mvarobj$signs <- a_signs
    mvarobj$edgecolor <- a_edgecolor
    mvarobj$rawlags[[lag]] <- l_lag
  }
  if (!saveModels) {
    mvarobj$nodemodels <- NULL
    mvarobj$rawfactor <- NULL
  }
  mvarobj$intercepts <- l_intercepts
  if (pbar) {
    if (signInfo) 
      cat("\nNote that the sign of parameter estimates is stored separately; see ?mvar")
  }
  else {
    if (signInfo) 
      cat("Note that the sign of parameter estimates is stored separately; see ?mvar")
  }
  class(mvarobj) <- c("mgm", "mvar")
  return(mvarobj)
}