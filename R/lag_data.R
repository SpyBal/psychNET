lag_data <- function(data, lag){
  ID <- NULL
  if(missing(lag))  lag <- 1
  uniqueIDs <- unique(data$ID)
  data_ID_list <- lapply(uniqueIDs, function(i){
    dt <- subset(data, ID==i)
    IDind <- which(colnames(dt) == "ID")
    dtnew <- dt[,-IDind]
    data_out <- cbind(dtnew,do.call(cbind, lapply(1:lag, function(j){
      lagged <- lagData(dtnew, j)$l_data_lags[[1]]
    })), ID=dt$ID)
    data_out
  })
  data_ID_list
}

