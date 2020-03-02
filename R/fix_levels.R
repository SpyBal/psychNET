fix_levels <- function(x, min.obs, lag){
  if (missing(min.obs)) min.obs <-  2
  if (missing(lag)) lag <-  0
  requirement <- min.obs + lag
  p <-  ncol(x)
  n <- nrow(x)
  for (i in 1:p){
    coli <- x[,i]
    tabci <- table(coli)
    levels <- as.numeric(names(tabci))
    llevels <- length(levels)
    ncomplete <-  sum(tabci)
    check_categories <- tabci < requirement
    index <- which(check_categories)
    lindex <- length(index)
    if (lindex==1){
      if (llevels>2){
        if (index==1){
          old <- levels[1]
          ind <- coli == old & !is.na(coli)
          coli[ind] <- levels[2]
        }
        if (index==llevels){
          old <- levels[llevels]
          ind <- coli == old & !is.na(coli)
          coli[ind] <- levels[llevels-1]
        }
        
        if (index!=llevels & index!=1){
          old <- levels[index]
          ind <- coli == old & !is.na(coli)
          replacement <- ifelse(which.min(tabci[as.character(c(levels[index-1],levels[index+1]))])==1,levels[index-1],levels[index+1])
          coli[ind] <- replacement
        }
      }
      
    }
    
    
    if (lindex>1){
      if (llevels>2){
        if (any(index==1)){
          old <- levels[1]
          ind <- coli == old & !is.na(coli)
          coli[ind] <- levels[2]
        }
        if (any(index==llevels)){
          tabci2 <- table(coli)
          levels2 <- as.numeric(names(tabci2))
          llevels2 <- length(levels2)
          if (llevels2>2){
            index2 <- llevels2
            lindex2 <- length(index2)
            old2 <- levels2[llevels2]
            ind2 <- coli == old2 & !is.na(coli)
            coli[ind2] <- levels2[llevels2-1]
          }
        }
        
        if (any(index!=llevels) & any(index!=1)){
          tabci2 <- table(coli)
          levels2 <- as.numeric(names(tabci2))
          llevels2 <- length(levels2)
          if (llevels2>2){
            index2 <- which(tabci2 < requirement)
            lindex2 <- length(index2)
            if (lindex2==1){
              old2 <- levels2[index2]
              ind2 <- coli == old2 & !is.na(coli)
              replacement <- ifelse(which.min(tabci2[as.character(c(levels2[index2-1],levels2[index2+1]))])==1,levels2[index2-1],levels2[index2+1])
              coli[ind2] <- replacement
            }
            if (lindex2>1) for (ind in index2) coli[coli == levels2[ind] & !is.na(coli)] <- ifelse(which.min(tabci2[as.character(c(levels2[ind-1],levels2[ind+1]))])==1,levels2[ind-1],levels2[ind+1])
          }
        }
      }
    }
    x[,i] <- coli
  }
  return(x)
}
