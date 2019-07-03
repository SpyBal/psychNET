ModelMatrix <- function (data, type, level, labels, d, moderators = NULL, v = NULL, allCats = FALSE) 
{
  if (!is.null(v)) {
    data <- data[, -v]
    type <- type[-v]
    level <- level[-v]
    labels <- labels[-v]
  }
  p <- ncol(data)
  n <- nrow(data)
  if (p != length(type)) 
    stop("Length of type has to match the number of columns in data.")
  if (p != length(level)) 
    stop("Length of level has to match the number of columns in data.")
  if (!(class(level) == "integer" | class(level) == "numeric")) 
    stop("level has to be an integer vector.")
  l_ind_datasets <- list()
  for (j in 1:p) {
    if (type[j] != "c") {
      l_ind_datasets[[j]] <- as.matrix(data[, j])
      colnames(l_ind_datasets[[j]]) <- paste0(labels[j])
    }
    else {
      if (allCats == FALSE) {
        unique_labels <- unique(data[, j])
        unique_labels_sorted <- sort(unique_labels)
        n_labels <- length(unique_labels_sorted)
        ind_matrix <- matrix(NA, nrow = n, ncol = n_labels)
        for (s in 1:n_labels) ind_matrix[, s] <- data[, 
                                                      j] == unique_labels_sorted[s]
      }
      else {
        unique_labels <- 1:level[j]
        n_labels <- level[j]
        ind_matrix <- matrix(NA, nrow = n, ncol = n_labels)
        for (s in 1:n_labels) ind_matrix[, s] <- data[, 
                                                      j] == unique_labels[s]
      }
      cn <- paste0(labels[j], 1:n_labels)
      colnames(ind_matrix) <- cn
      l_ind_datasets[[j]] <- ind_matrix
    }
  }
  l_ind_datasets_nV <- l_ind_datasets
  Xd1 <- do.call(cbind, l_ind_datasets_nV)
  if (!is.null(moderators)) 
    d <- 2
  if (d > 1) {
    l_interactions <- vector("list", length = d)
    if (is.null(moderators)) {
      for (ord in 1:d) l_interactions[[ord]] <- combn((1:p), 
                                                      ord, simplify = FALSE)
    }
    else {
      l_interactions[[1]] <- list()
      for (i in 1:p) l_interactions[[1]][[i]] <- i
      if (v %in% moderators) {
        l_interactions[[2]] <- combn(1:p, 2, simplify = FALSE)
      }
      else {
        ind_mod_MM <- (1:(p + 1) %in% moderators)[-v]
        n_mods <- sum(ind_mod_MM)
        which_mod <- which(ind_mod_MM)
        l_mods <- list()
        for (i in 1:n_mods) l_mods[[i]] <- expand.grid((1:p)[-which_mod[i]], 
                                                       which_mod[i])
        m_mods <- do.call(rbind, l_mods)
        m_mods <- as.matrix(m_mods)
        l_mods_combn <- list()
        for (i in 1:nrow(m_mods)) l_mods_combn[[i]] <- m_mods[i, 
                                                              ]
        l_interactions[[2]] <- l_mods_combn
      }
    }
    l_collect_terms <- list()
    for (ord in 2:d) {
      n_terms <- length(l_interactions[[ord]])
      l_ord_terms <- list()
      for (it in 1:n_terms) {
        l_it_terms <- list()
        inter_it <- l_interactions[[ord]][[it]]
        l_indicator_it <- list()
        for (i in 1:ord) l_indicator_it[[i]] <- 1:level[inter_it[i]]
        all_combs <- expand.grid(l_indicator_it)
        n_combs <- nrow(all_combs)
        for (comb in 1:n_combs) {
          tarmat <- matrix(NA, nrow = n, ncol = ord)
          for (i in 1:ord) tarmat[, i] <- as.matrix(l_ind_datasets[[inter_it[i]]])[, 
                                                                                   all_combs[comb, i]]
          l_it_terms[[comb]] <- apply(tarmat, 1, prod)
        }
        it_data <- do.call(cbind, l_it_terms)
        all_combs_char <- apply(all_combs, 2, function(x) {
          if (length(unique(x)) == 1) {
            out <- rep("", length(x))
          }
          else {
            out <- x
          }
        })
        cn <- rep(NA, n_combs)
        for (comb in 1:n_combs) cn[comb] <- paste0(labels[inter_it], 
                                                   matrix(all_combs_char, ncol = ord)[comb, ], 
                                                   collapse = ":")
        colnames(it_data) <- cn
        l_ord_terms[[it]] <- it_data
      }
      l_collect_terms[[ord]] <- do.call(cbind, l_ord_terms)
    }
    all_HOI_terms <- do.call(cbind, l_collect_terms)
  }
  if (d > 1) {
    X <- cbind(Xd1, all_HOI_terms)
  }
  else {
    X <- Xd1
  }
  return(X)
}