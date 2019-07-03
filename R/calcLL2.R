calcLL2 <- function (X, y, fit, type, level, v, weights, lambda, LLtype = "model") 
{
  if (missing(level)) stop("No levels passed to calcLL !")
  n <- nrow(X)
  if (LLtype == "model") {
    if (type[v] == "g") {
      beta_vector <- matrix(coef(fit, s = lambda), ncol = 1)
      predicted_mean <- cbind(rep(1, n), X) %*% as.vector(beta_vector)
      LL_model <- dnorm(y, mean = predicted_mean, sd = 1, 
                        log = TRUE)
      mean_LL_model <- sum(LL_model * weights)
    }
    if (type[v] == "p") {
      beta_vector <- matrix(coef(fit, s = lambda), ncol = 1)
      predicted_mean <- cbind(rep(1, n), X) %*% as.vector(beta_vector)
      LL_model <- dpois(y, exp(predicted_mean), log = TRUE)
      mean_LL_model <- sum(LL_model * weights)
    }
    if (type[v] == "c") {
      n_cats <- level[v]
      m_respdum <- matrix(NA, n, n_cats)
      m_coefs <- matrix(NA, n, n_cats)
      cats <- unique(y)
      LL_n <- rep(NA, n)
      m_LL_parts <- matrix(NA, nrow = n, ncol = n_cats + 
                             1)
      for (catIter in 1:n_cats) {
        m_respdum[, catIter] <- (y == cats[catIter]) * 
          1
        m_coefs[, catIter] <- cbind(rep(1, n), X) %*% 
          matrix(coef(fit, s = lambda)[[catIter]], ncol = 1)
        m_LL_parts[, catIter] <- m_respdum[, catIter] * 
          m_coefs[, catIter]
      }
      m_LL_parts[, n_cats + 1] <- -log(rowSums(exp(m_coefs)))
      LL_n <- rowSums(m_LL_parts)
      mean_LL_model <- sum(LL_n * weights)
    }
    if (type[v] == "o") {
      n_cats <- length(levels(y))
      m_respdum <- matrix(NA, n, (n_cats-1))
      m_coefs <- matrix(NA, n, (n_cats-1))
      cats <- sort(unique(y))
      LL_n <- rep(NA, n)
      m_LL_parts <- matrix(NA, nrow = n, ncol = n_cats)
      for (catIter in 1:(n_cats-1)) {
        m_respdum[, catIter] <- (y == cats[catIter]) * 1
        m_coefs[, catIter] <- cbind(rep(1, n), X) %*% coef(fit, whichLambda = 1,matrix=TRUE)[,catIter]
        m_LL_parts[, catIter] <- m_respdum[, catIter] * m_coefs[, catIter]
      }
      m_LL_parts[, n_cats ] <- -log(rowSums(exp(m_coefs)))
      LL_n <- rowSums(m_LL_parts)
      mean_LL_model <- sum(LL_n * weights[1:length(LL_n)])
    }
  }
  if (LLtype == "nullmodel") {
    if (type[v] == "g") {
      beta_vector <- matrix(coef(fit, s = 1)[1], ncol = 1)
      predicted_mean <- rep(1, n) * as.vector(beta_vector)
      LL_model <- dnorm(y, mean = predicted_mean, sd = 1, 
                        log = TRUE)
      mean_LL_model <- sum(LL_model * weights)
    }
    if (type[v] == "p") {
      beta_vector <- matrix(coef(fit, s = 1)[1], ncol = 1)
      predicted_mean <- rep(1, n) * as.vector(beta_vector)
      LL_model <- dpois(y, exp(predicted_mean), log = TRUE)
      mean_LL_model <- sum(LL_model * weights)
    }
    if (type[v] == "c") {
      n_cats <- level[v]
      m_respdum <- matrix(NA, n, n_cats)
      m_coefs <- matrix(NA, n, n_cats)
      cats <- sort(unique(y))
      LL_n <- rep(NA, n)
      m_LL_parts <- matrix(NA, nrow = n, ncol = n_cats + 
                             1)
      for (catIter in 1:n_cats) {
        m_respdum[, catIter] <- (y == cats[catIter]) * 
          1
        m_coefs[, catIter] <- cbind(rep(1, n), X) %*% 
          matrix(coef(fit, s = 1)[[catIter]], ncol = 1)
        m_LL_parts[, catIter] <- m_respdum[, catIter] * 
          m_coefs[, catIter]
      }
      m_LL_parts[, n_cats + 1] <- -log(rowSums(exp(m_coefs)))
      LL_n <- rowSums(m_LL_parts)
      mean_LL_model <- sum(LL_n * weights)
    }
    if (type[v]=="o"){
      ffit <-  ordinalNet(x = X, y = y, family = "cumulative", link="logit",lambdaVals = 1)
      n_cats <- length(levels(y))
      m_respdum <- matrix(NA, n, (n_cats-1))
      m_coefs <- matrix(NA, n, (n_cats-1))
      cats <- sort(unique(y))
      LL_n <- rep(NA, n)
      m_LL_parts <- matrix(NA, nrow = n, ncol = n_cats )
      for (catIter in 1:(n_cats-1)) {
        m_respdum[, catIter] <- (y == cats[catIter]) * 
          1
        m_coefs[, catIter] <- cbind(rep(1, n), X) %*% 
          coef(ffit, whichLambda = 1, matrix = TRUE)[, catIter]
        m_LL_parts[, catIter] <- m_respdum[, catIter] * 
          m_coefs[, catIter]
      }
      m_LL_parts[, n_cats] <- -log(rowSums(exp(m_coefs)))
      LL_n <- rowSums(m_LL_parts)
      mean_LL_model <- sum(LL_n * weights[1:length(LL_n)])
    }
  }
  if (LLtype == "saturated") {
    if (type[v] == "g") {
      predicted_mean <- y
      LL_model <- dnorm(y, mean = predicted_mean, sd = 1, 
                        log = TRUE)
      mean_LL_model <- sum(LL_model * weights)
    }
    if (type[v] == "p") {
      predicted_mean <- y
      LL_model <- dpois(y, exp(predicted_mean), log = TRUE)
      mean_LL_model <- sum(LL_model * weights)
    }
    if (type[v] == "c") {
      mean_LL_model <- 0
    }
    if (type[v] == "o") {
      mean_LL_model <- 0
    }
  }
  return(mean_LL_model)
}