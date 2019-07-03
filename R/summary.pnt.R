summary.pnt <- function(object, ...){
  models_all <- models_in()
  summary_list_result <- list()
  summary_list_result$temporal <- summary_list_result$contemporaneous <- list()
  
  if (object$CALL$pars$model == "DFM" && (object$CALL$pars$lag>1 | object$CALL$pars$nFact==1)){
    stop("Summary is not available when model='DFM' and lag > 1 or nFact=1")
  }else{
    p <- object$CALL$pars$lag
    if (p!=1){
      igraphsDIR <- lapply(1:p, function(i){
        dn <- graph_from_adjacency_matrix(object$results$Dir_net[[i]], mode ="directed",diag=FALSE, weighted = TRUE,
                                          add.colnames = NULL, add.rownames = NA)
        E(dn)$color <- ifelse(E(dn)$weight > 0,'green','red')
        dn})
    }else{
      dn <- graph_from_adjacency_matrix(object$results$Dir_net[[1]], mode ="directed",diag=FALSE, weighted = TRUE,
                                        add.colnames = NULL, add.rownames = NA)
      E(dn)$color <- ifelse(E(dn)$weight > 0,'green','red')
      E(dn)$weight <- E(dn)$weight
      igraphsDIR <- list()
      igraphsDIR[[1]] <- dn
    }
    object$CALL$pars$lag_used <- 1:object$CALL$pars$lag
    
    if (object$CALL$pars$model %in% models_in()$sparse){
      check_empt <- sapply(igraphsDIR, function(x) is.null(E(x)$weight))
      if (sum(check_empt)== object$CALL$pars$lag) return("All temporal graphs are empty.")
      if (any(check_empt)){
        lgs <- 1:object$CALL$pars$lag
        object$CALL$pars$lag_used <- lgs[which(!check_empt)]
        igraphsDIR <-  igraphsDIR[!check_empt]
      } 
    }

    igraphsDIR_glob <- do.call(cbind,lapply(igraphsDIR, function(i) {
      trans <- transitivity(i, type = "global", isolates = "zero")
      rs1 <- reciprocity(i, mode = "default")
      rs2 <- reciprocity(i, mode = "ratio")
      md <- mean_distance(i, directed = TRUE )
      dens <-  edge_density(i, loops = FALSE)
      E(i)$weight <- abs(E(i)$weight)
      diam <- diameter(i, directed = TRUE)
      c(trans, rs1,rs2, md, dens, diam)
    }))
    
    dimnames(igraphsDIR_glob) <- list(c("Transitivity","Resiprocity","Resiprocity_ratio",
                                        "Mean_distance", "Density","Diameter"),paste("lag_",object$CALL$pars$lag_used,sep = ""))
    summary_list_result$temporal$global <- round(igraphsDIR_glob,3)
    
    igraphsDIR_local <- lapply(igraphsDIR, function(g) {
      descripts <- rbind(transitivity(g, type = "local", isolates = "zero"),
            transitivity(g, type = "weighted", isolates = "zero"),
            degree(g, mode = "in"),
            degree(g, mode = "out"),
            degree(g, mode = "total"),
            expectedInf(g, step = 1, directed = TRUE)$step1,
            expectedInf(g, step = 2, directed = TRUE)$step2)
      E(g)$weight <- abs(E(g)$weight)
      descripts <- rbind(descripts, betweenness(g, directed = TRUE))
      descripts <- rbind(descripts, closeness(g, mode = "out"))
      descripts <- rbind(descripts, closeness(g, mode = "in"))
      descripts <- round(descripts,3)
      dimnames(descripts)[[1]] <- c("Transitivity", "Transitivity_weight","Indegree","Outdegree",
                                    "Degree","Step-1 expected influence", "Step-2 expected influence",
                                    "Betweeness", "Outcloseness","Incloseness" )
      descripts})
    
    summary_list_result$temporal$local <- igraphsDIR_local

    if (object$CALL$pars$model %in% c("DFM","GVAR","GGVAR","MLVAR")){
      igraphUnDIR <- graph_from_adjacency_matrix(object$results$UnDir_net, mode ="undirected",diag=FALSE, weighted = TRUE,
                                          add.colnames = NULL, add.rownames = NA)
      if (is.null(E(igraphUnDIR)$weight)){
        summary_list_result$contemporaneous <- "Not available"
      }else{
        E(igraphUnDIR)$color <- ifelse(E(igraphUnDIR)$weight > 0,'green','red')
        E(igraphUnDIR)$weight <- E(igraphUnDIR)$weight
        
        y <- igraphUnDIR
        E(y)$weight <- abs(E(y)$weight)
        
        igraphsUnDIR_glob <- matrix(c(
          transitivity(igraphUnDIR, type = "global", isolates = "zero"), 
          transitivity(igraphUnDIR, type = "average", isolates = "zero"),
          reciprocity(igraphUnDIR, mode = "default"),
          reciprocity( igraphUnDIR, mode = "ratio"),
          mean_distance(igraphUnDIR, directed = FALSE ),
          edge_density(igraphUnDIR, loops = FALSE),
          diameter(y, directed = FALSE)),ncol = 1)
        dimnames(igraphsUnDIR_glob) <- list(c("Transitivity","Transitivity_ave","Resiprocity","Resiprocity_ratio",
                                              "Mean_distance", "Density","Diameter"),"value" )
        summary_list_result$contemporaneous$global <- round(igraphsUnDIR_glob,3)
        
        
        igraphsUnDIR_local <- rbind(
          transitivity(igraphUnDIR, type = "local", isolates = "zero"),
          transitivity(igraphUnDIR, type = "weighted", isolates = "zero"), 
          degree(igraphUnDIR, mode = "total"),
          expectedInf(igraphUnDIR, step = 1, directed = FALSE)$step1,
          expectedInf(igraphUnDIR, step = 2, directed = FALSE)$step2,
          betweenness(y, directed = TRUE),
          closeness(y, mode = "total"))
        igraphsUnDIR_local <- round(igraphsUnDIR_local,3)
        dimnames(igraphsUnDIR_local)[[1]] <- c("Transitivity", "Transitivity_weight",
                                               "Degree","Step-1 expected influence", "Step-2 expected influence",
                                               "Betweeness","Closeness")
        
        summary_list_result$contemporaneous$local <- igraphsUnDIR_local
        
      }
    }else{
      summary_list_result$contemporaneous <- "Not available"
    }
  }
  return(summary_list_result)
}