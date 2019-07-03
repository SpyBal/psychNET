lagData <- function (data, lags, consec = NULL) 
{
  data <- as.matrix(data)
  max_lag <- max(lags)
  lags_ext <- 1:max(lags)
  n <- nrow(data)
  p <- ncol(data)
  n_var <- nrow(data) - max(lags)
  n_lags <- length(lags_ext)
  data_response <- data
  if (!is.null(consec)) 
    m_consec <- matrix(NA, nrow = n, ncol = n_lags)
  l_data_lags <- list()
  lag_pos <- 1
  for (lag in lags) {
    lagged_data <- matrix(NA, nrow = n, ncol = p)
    lagged_data[(lag + 1):n, ] <- data[-((n - lag + 1):n), 
                                       ]
    lagged_data <- matrix(lagged_data, ncol = p, nrow = n)
    colnames(lagged_data) <- paste("V", 1:p, ".lag", 
                                   lag, ".", sep = "")
    l_data_lags[[lag_pos]] <- lagged_data
    lag_pos <- lag_pos + 1
  }
  if (!is.null(consec)) {
    for (lag in lags_ext) m_consec[(lag + 1):n, lag] <- consec[-((n - 
                                                                    lag + 1):n)]
    m_consec_check <- cbind(consec, m_consec)
    v_check <- apply(m_consec_check, 1, function(x) {
      if (any(is.na(x))) {
        FALSE
      }
      else {
        check_row <- x[1] - x[-1] == 1:length(x[-1])
        check_row_relevant <- check_row[lags_ext %in% 
                                          lags]
        if (any(check_row_relevant == FALSE)) 
          FALSE
        else TRUE
      }
    })
  }
  else {
    v_check <- rep(TRUE, n)
    v_check[1:n_lags] <- FALSE
  }
  outlist <- list()
  outlist$data_response <- data_response
  outlist$l_data_lags <- l_data_lags
  outlist$included <- v_check
  return(outlist)
}