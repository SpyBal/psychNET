psychoMLVAR <- function (x, type = c("temporal", "contemporaneous", "between"), 
                         lag = 1, partial = TRUE, SD = FALSE, subject, order, nonsig = c("default", 
                                                                                         "show", "hide", "dashed"), rule = c("or", "and"), alpha = 0.05, 
                         onlySig = FALSE, layout = "spring", verbose = TRUE, ...) 
{
  rule <- match.arg(rule)
  if (type[[1]] == "fixed") {
    warning("type = 'fixed' is deprecated. type = 'temporal' instead.")
    type <- "temporal"
  }
  if (type[[1]] == "SD") {
    warning("type = 'SD' is deprecated. type = 'temporal' and SD = TRUE instead.")
    type <- "temporal"
    SD <- TRUE
  }
  if (type[[1]] == "subject") {
    warning("type = 'subject' is deprecated. type = 'temporal' and subject = ... instead")
    type <- "temporal"
    if (missing(subject)) 
      stop("'subject' must be assigned")
  }
  if (onlySig) {
    warning("'onlySig' is deprecated. Setting nonsig = 'hide'.")
    nonsig <- "hide"
  }
  type <- match.arg(type)
  nonsig <- match.arg(nonsig)
  if (nonsig == "default") {
    nonsig <- "hide"
    if (!partial) {
      nonsig <- "show"
    }
    if (!missing(subject)) {
      nonsig <- "show"
    }
    if (verbose) {
      message(paste0("'nonsig' argument set to: '", nonsig, 
                     "'"))
    }
  }
  if (missing(order)) {
    order <- x$input$vars
  }
  if (!missing(subject) && SD) {
    stop("'SD' not available for subject.")
  }
  if (is.character(order)) {
    ord <- match(order, x$input$vars)
  }
  else {
    ord <- order
  }
  if (type == "temporal") {
    if (SD) {
      SIG <- matrix(TRUE, length(ord), length(ord))
      NET <- t(x$results$Beta$SD[, , lag])
    }
    else {
      if (missing(subject)) {
        NET <- t(x$results$Beta$mean[, , lag])
        if (any(is.na(x$results$Beta$P))) {
          if (!any(is.na(x$results$Beta$lower)) && !any(is.na(x$results$Beta$upper))) {
            SIG <- t(x$results$Beta$lower[, , lag]) > 
              0 | t(x$results$Beta$upper[, , lag]) < 
              0
          }
          else {
            SIG <- matrix(TRUE, nrow(NET), ncol(NET))
            if (nonsig != "show") {
              warning("No p-values or CI computed. Can not hide non-significant edges.")
            }
          }
        }
        else {
          SIG <- t(x$results$Beta$P[, , lag]) < alpha
        }
      }
      else {
        NET <- t(x$results$Beta$subject[[subject]][, 
                                                   , lag])
        SIG <- matrix(TRUE, nrow(NET), ncol(NET))
        if (nonsig != "show") {
          warning("Can not hide non-significant edges for subject network.")
        }
      }
    }
  }
  if (type == "contemporaneous") {
    sub <- ifelse(partial, "pcor", "cor")
    if (missing(subject)) {
      if (SD) {
        NET <- x$results$Theta[[sub]]$SD
      }
      else {
        NET <- x$results$Theta[[sub]]$mean
      }
      if (nonsig != "show") {
        if (SD) {
          stop("No significance for SD")
        }
        if (!any(is.na(x$results$Theta[[sub]]$lower)) && 
            !any(is.na(x$results$Theta[[sub]]$upper))) {
          SIG <- x$results$Theta[[sub]]$lower > 0 | x$results$Theta[[sub]]$upper < 
            0
        }
        else if (!any(is.na(x$results$Theta[[sub]]$P))) {
          SIG <- x$results$Theta[[sub]]$P < alpha
        }
        else {
          if (partial && !is.null(x$results$Gamma_Theta) && 
              !all(is.nan(x$results$Gamma_Theta$P))) {
            diag(x$results$Gamma_Theta$P) <- 0
            if (rule == "or") {
              SIG <- x$results$Gamma_Theta$P < alpha | 
                t(x$results$Gamma_Theta$P) < alpha
            }
            else {
              SIG <- x$results$Gamma_Theta$P < alpha & 
                t(x$results$Gamma_Theta$P) < alpha
            }
          }
          else {
            SIG <- matrix(TRUE, nrow(NET), ncol(NET))
            if (nonsig != "show") {
              stop("No p-values or CI computed. Can not hide non-significant edges.")
            }
          }
        }
      }
      else {
        SIG <- matrix(TRUE, nrow(NET), ncol(NET))
      }
    }
    else {
      NET <- x$results$Theta[[sub]]$subject[[subject]]
      SIG <- matrix(TRUE, nrow(NET), ncol(NET))
      if (nonsig != "show") {
        warning("Can not hide non-significant edges for subject network.")
      }
    }
    NET <- makeSym(NET)
  }
  if (type == "between") {
    sub <- ifelse(partial, "pcor", "cor")
    if (!missing(subject)) {
      stop("No subject-specific between network possible")
    }
    if (SD) {
      stop("No SD for between-subjects network.")
    }
    NET <- x$results$Omega_mu[[sub]]$mean
    if (nonsig != "show") {
      if (SD) {
        stop("No significance for SD")
      }
      if (!any(is.na(x$results$Omega_mu[[sub]]$lower)) && 
          !any(is.na(x$results$Omega_mu[[sub]]$upper))) {
        SIG <- x$results$Omega_mu[[sub]]$lower > 0 | 
          x$results$Omega_mu[[sub]]$upper < 0
      }
      else if (!any(is.na(x$results$Omega_mu[[sub]]$P))) {
        SIG <- x$results$Omega_mu[[sub]]$P < alpha
      }
      else {
        if (partial && !is.null(x$results$Gamma_Omega_mu) && 
            !all(is.nan(x$results$Gamma_Omega_mu$P))) {
          diag(x$results$Gamma_Omega_mu$P) <- 0
          if (rule == "or") {
            SIG <- x$results$Gamma_Omega_mu$P < alpha | 
              t(x$results$Gamma_Omega_mu$P) < alpha
          }
          else {
            SIG <- x$results$Gamma_Omega_mu$P < alpha & 
              t(x$results$Gamma_Omega_mu$P) < alpha
          }
        }
        else {
          SIG <- matrix(TRUE, nrow(NET), ncol(NET))
          if (nonsig != "show") {
            stop("No p-values or CI computed. Can not hide non-significant edges.")
          }
        }
      }
    }
    else {
      SIG <- matrix(TRUE, nrow(NET), ncol(NET))
    }
    NET <- makeSym(NET)
  }
  if (nonsig == "dashed") {
    lty <- ifelse(!SIG, 2, 1)
  }
  else {
    lty <- 1
  }
  if (nonsig == "hide") {
    NET <- NET * SIG
  }
  if (any(is.na(NET[ord, ord][upper.tri(NET[ord, ord])])) || 
      any(is.na(NET[ord, ord][lower.tri(NET[ord, ord])]))) {
    stop("Network not estimated correctly.")
  }
  NET[ord, ord]
}