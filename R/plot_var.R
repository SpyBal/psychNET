plot_var <- function(x, plot.type, communit,...){
  if (missing(plot.type)) plot.type <- "temporal"
  if (missing(communit)) communit <- FALSE
  
  type <- plot.type
  if (type=="temporal"){
    igraphsDIR <- lapply(1:x$CALL$pars$lag, function(i){
      dn <- graph_from_adjacency_matrix(x$results$Dir_net[[i]], mode ="directed",diag=FALSE, weighted = TRUE,
                                        add.colnames = NULL, add.rownames = NA)
      if(gsize(dn)==0){
        temporal <- list(qgraph = paste("The lag",i,"temporal network is empty"), 
                         igraph = paste("The lag",i,"temporal network is empty"))
      }else{
        dn <- delete.vertices(dn, degree(dn)==0)
        E(dn)$color <- ifelse(E(dn)$weight > 0,'green','red')
        E(dn)$weight <- E(dn)$weight
        qgobj <- qgraph(as_adjacency_matrix(dn,attr = "weight"), DoNotPlot=TRUE )
        plot(dn,vertex.label=V(dn)$name,
             vertex.label.color="black",curved=TRUE,
             vertex.size=12,
             edge.arrow.size=0.5,
             edge.curved=0.15,
             vertex.color="deepskyblue2",
             layout = qgobj$layout,
             main=paste("Lag",i,"temporal network from individual", x$CALL$pars$data_pre_fixed$ID[1]),...)
        temporal <- list(qgraph=qgobj, igraph = dn)
        if (communit){
          cm <- cluster_spinglass(dn, stop.temp = 0.000000000001, cool.fact = 0.999,implementation = "neg", update.rule = "config", gamma = 0, gamma.minus = 0)
          temporal <- list(qgraph=qgobj, igraph = dn, communities = cm)
          plot(cm, dn, layout=qgobj$layout,
               main=paste("Spinglass communities in \n lag",i,"temporal network from individual", x$CALL$pars$data_pre_fixed$ID[1]),...)
        }
      }
      temporal})
    return(igraphsDIR)
  }
  if (type=="contemporaneous"){
    stop("Contemporaneous network is not available when model = 'VAR', 'SVAR', 'SVECM', 'SVARHL', 'SVARHLX', 'SVARHLMA',
         'SMVAR'.")
  }
  if (type=="between"){
    stop("Between subjects network is available only when model= 'MLVAR")
  }
  if (type=="both"){
    igraphsDIR <- lapply(1:x$CALL$pars$lag, function(i){
      dn <- graph_from_adjacency_matrix(x$results$Dir_net[[i]], mode ="directed",diag=FALSE, weighted = TRUE,
                                        add.colnames = NULL, add.rownames = NA)
      if(gsize(dn)==0){
        temporal <- list(qgraph = paste("The lag",i,"temporal network is empty"), 
                         igraph = paste("The lag",i,"temporal network is empty"))
      }else{
        dn <- delete.vertices(dn, degree(dn)==0)
        E(dn)$color <- ifelse(E(dn)$weight > 0,'green','red')
        E(dn)$weight <- E(dn)$weight
        qgobj <- qgraph(as_adjacency_matrix(dn,attr = "weight"), DoNotPlot=TRUE )
        plot(dn,vertex.label=V(dn)$name,
             vertex.label.color="black",curved=TRUE,
             vertex.size=12,
             edge.arrow.size=0.5,
             edge.curved=0.15,
             vertex.color="deepskyblue2",
             layout = qgobj$layout,
             main=paste("Lag",i,"temporal network from individual", x$CALL$pars$data_pre_fixed$ID[1]),...)
        temporal <- list(qgraph=qgobj, igraph = dn)
        if (communit){
          cm <- cluster_spinglass(dn, stop.temp = 0.000000000001, cool.fact = 0.999,implementation = "neg", update.rule = "config", gamma = 0, gamma.minus = 0)
          temporal <- list(qgraph=qgobj, igraph = dn, communities = cm)
          plot(cm, dn, layout=qgobj$layout,
               main=paste("Spinglass communities in \n lag",i,"temporal network from individual", x$CALL$pars$data_pre_fixed$ID[1]),...)
        }
        message("Contemporaneous network is not available when model = 'VAR', 'SVAR', 'SVECM', 'SVARHL', 'SVARHLX', 'SVARHLMA',
         'SMVAR'.")  
      }
      temporal})
    return(igraphsDIR)
  }
}