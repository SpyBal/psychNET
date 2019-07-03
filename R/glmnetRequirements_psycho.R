glmnetRequirements_psycho <- function (data, type, weights) 
{
  var_names <- colnames(data)
  n <- nrow(data)
  var_check <- apply(data, 2, var)
  ind_zero_var <- which(var_check == 0)
  if (length(ind_zero_var) > 0) 
    stop(paste0("Please only provide variables with nonzero variance. Variable(s) with zero variance: ", 
                paste(var_names[ind_zero_var], collapse = ", ")))
  if ("c" %in% type) {
    ind_cat <- which(type == "c")
    l_frqu <- list()
    for (i in 1:length(ind_cat)) l_frqu[[i]] <- table(data[, ind_cat[i]])
    v_check <- unlist(lapply(l_frqu, function(x) {
      frq_norm <- x/sum(x)
      ind_min <- which.min(frq_norm)
      check1 <- !(x[ind_min] > 1)
    }))
    if (sum(v_check) > 0) {
      ind_check1 <- ind_cat[v_check == TRUE]
      stop(paste0("At least 2 events required for each category. Requirement not met for variable(s): ", 
                  paste(var_names[ind_cat[ind_check1]], collapse = ", ")))
    }
  }
  if ("c" %in% type) {
    wtable <- function(x, weights) {
      n_level <- length(unique(x))
      v_level <- unique(x)
      n_obs <- length(x)
      v_wfrq <- rep(NA, n_level)
      for (i in 1:n_level) v_wfrq[i] <- sum(rep(1, sum(x == 
                                                         v_level[i])) * weights[x == v_level[i]])/n
      return(v_wfrq)
    }
    ind_cat <- which(type == "c")
    check2 <- rep(NA, length(ind_cat))
    for (i in 1:length(ind_cat)) check2[i] <- min(wtable(data[, 
                                                              ind_cat[i]], weights)) < 10^-5
    if (sum(check2) > 0) 
      stop(paste0("Each category has to have probability > 10^-5. Requirement not met for variable(s): ", 
                  paste(var_names[ind_cat[check2]], collapse = ", ")))
  }
}
