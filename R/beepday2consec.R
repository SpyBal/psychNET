beepday2consec <-  function (beepvar, dayvar) 
{
  if (!all(dayvar == round(dayvar))) 
    stop("beepvar has to be a vector of non-negative integers")
  if (!all(beepvar == round(beepvar))) 
    stop("beepvar has to be a vector of non-negative integers")
  if (length(beepvar) != length(dayvar)) 
    stop("beepvar has to have the same length as dayvar")
  fillin_beep <- max(beepvar) + 2
  n <- length(dayvar)
  ind_sameday <- dayvar[-1] == dayvar[-n]
  ind_sameday[1] <- TRUE
  consec <- rep(NA, n)
  consec[1] <- 1
  counter <- 1
  for (i in 2:n) {
    beep_diff <- beepvar[i] - beepvar[i - 1]
    day_diff <- dayvar[i] - dayvar[i - 1]
    if (beep_diff == 1 & day_diff == 0) 
      counter <- counter + 1
    else counter <- counter + 2
    consec[i] <- counter
  }
  return(consec)
}