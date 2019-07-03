plot_gvar <- function(x, plot.type, communit,...){
  plot_out <- list()
  

  if (missing(plot.type)) plot.type <- "temporal"
  if (missing(communit)) communit <- FALSE
  qgobj <- NULL
  type <- plot.type
  
  if (type=="temporal"){
    plot_out$temporal$qgraph <- list()
    plot_out$temporal$igraph <- list()
    igraphsDIR <- lapply(1:x$CALL$pars$lag, function(i){

      dn <- graph_from_adjacency_matrix(x$results$Dir_net[[i]], mode ="directed",diag=FALSE, weighted = TRUE,
                                        add.colnames = NULL, add.rownames = NA)
      E(dn)$color <- ifelse(E(dn)$weight > 0,'green','red')
      E(dn)$weight <- E(dn)$weight
      qgobj <- qgraph(as_adjacency_matrix(dn,attr = "weight"), DoNotPlot=TRUE )
      plot_out$temporal$qgraph[[i]] <- qgobj
      plot_out$temporal$igraph[[i]] <- dn
      plot(dn,
           vertex.label=V(dn)$name,
           vertex.label.color="black",curved=TRUE,
           vertex.size=12,
           edge.arrow.size=0.5,
           vertex.color="deepskyblue2",
           edge.curved=0.15,
           layout = qgobj$layout,
           main=ifelse(x$CALL$pars$model=="GVAR",
                       paste("Lag",i,"temporal network from individual",x$CALL$pars$data_pre_fixed$ID[1]),
                       paste("Lag",i,"temporal network")),...)
      if (communit){
        plot_out$temporal$communities[[i]] <- cluster_spinglass(dn, stop.temp = 0.000000000001, cool.fact = 0.999,
                                                                implementation = "neg", update.rule = "config", gamma = 0, gamma.minus = 0)
        plot(plot_out$temporal$communities[[i]], dn, layout=qgobj$layout,
             main=ifelse(x$CALL$pars$model=="GVAR",
                         paste("Spinglass communities in\n lag",i,"temporal network from individual",x$CALL$pars$data_pre_fixed$ID[1]),
                         paste("Spinglass communities in\n lag",i,"temporal network")),...)
      }
      plot_out})
    return(igraphsDIR)
  }
  if (type=="contemporaneous"){
    plot_out$contemporaneous$igraph <-list()
    plot_out$contemporaneous$qgraph <-list()
    cn <- graph_from_adjacency_matrix(x$results$UnDir_net, mode ="undirected",diag=FALSE, weighted = TRUE,
                                      add.colnames = NULL, add.rownames = NA)
    if (x$CALL$pars$model %in% models_in()$sparse){
      check_empt <- is.null(E(cn)$weight)
      if (check_empt) return("Contemporaneous graph is empty.")
    }
    E(cn)$color <- ifelse(E(cn)$weight > 0,'green','red')
    qgobjcn <- qgraph(as_adjacency_matrix(cn,attr = "weight"), DoNotPlot=TRUE )

    plot_out$contemporaneous$qgraph[[1]] <- qgobjcn
    plot_out$contemporaneous$igraph[[1]] <- cn
    plot(cn,
         vertex.label=V(cn)$name,
         vertex.label.color="black",curved=TRUE,
         vertex.size=12,
         edge.arrow.size=0.5,
         vertex.color="deepskyblue2",
         layout = qgobjcn$layout,
         edge.curved=0.15,
         main=ifelse(x$CALL$pars$model=="GVAR",
                     paste("Contemporaneous network from individual", x$CALL$pars$data_pre_fixed$ID[1]),
                     "Contemporaneous network"),...)
    if (communit){
      plot_out$contemporaneous$communities[[1]] <- cluster_spinglass(cn, stop.temp = 0.000000000001, cool.fact = 0.999,
                                                              implementation = "neg", update.rule = "config", gamma = 0, gamma.minus = 0)
      plot(plot_out$contemporaneous$communities[[1]], cn, layout=qgobj$layout,
           main=ifelse(x$CALL$pars$model=="GVAR",
                       paste("Spinglass communities in the \ncontemporaneous network from individual", x$CALL$pars$data_pre_fixed$ID[1]),
                       "Spinglass communities in the \ncontemporaneous network"),...)
    }
    return(plot_out)
  }
  if (type=="between"){
    stop("Between subjects network is available only when model= 'MLVAR")
  }
  
  if (type=="both"){
    plot_out$temporal$qgraph <- list()
    plot_out$temporal$igraph <- list()
    igraphsDIR <- lapply(1:x$CALL$pars$lag, function(i){
      
      dn <- graph_from_adjacency_matrix(x$results$Dir_net[[i]], mode ="directed",diag=FALSE, weighted = TRUE,
                                        add.colnames = NULL, add.rownames = NA)
      E(dn)$color <- ifelse(E(dn)$weight > 0,'green','red')
      E(dn)$weight <- E(dn)$weight
      qgobj <- qgraph(as_adjacency_matrix(dn,attr = "weight"), DoNotPlot=TRUE )
      plot_out$temporal$qgraph[[i]] <- qgobj
      plot_out$temporal$igraph[[i]] <- dn
      plot(dn,
           vertex.label=V(dn)$name,
           vertex.label.color="black",curved=TRUE,
           vertex.size=12,
           edge.arrow.size=0.5,
           vertex.color="deepskyblue2",
           edge.curved=0.15,
           layout = qgobj$layout,
           main=ifelse(x$CALL$pars$model=="GVAR",
                       paste("Lag",i,"temporal network from individual",x$CALL$pars$data_pre_fixed$ID[1]),
                       paste("Lag",i,"temporal network")),...)
      if (communit){
        plot_out$temporal$communities[[i]] <- cluster_spinglass(dn, stop.temp = 0.000000000001, cool.fact = 0.999,
                                                                implementation = "neg", update.rule = "config", gamma = 0, gamma.minus = 0)
        plot(plot_out$temporal$communities[[i]], dn, layout=qgobj$layout,
             main=ifelse(x$CALL$pars$model=="GVAR",
                         paste("Spinglass communities in\n lag",i,"temporal network from individual",x$CALL$pars$data_pre_fixed$ID[1]),
                         paste("Spinglass communities in\n lag",i,"temporal network")),,...)
      }})
    plot_out$contemporaneous$igraph <-list()
    plot_out$contemporaneous$qgraph <-list()
    cn <- graph_from_adjacency_matrix(x$results$UnDir_net, mode ="undirected",diag=FALSE, weighted = TRUE,
                                      add.colnames = NULL, add.rownames = NA)
    E(cn)$color <- ifelse(E(cn)$weight > 0,'green','red')
    E(cn)$weight <- E(cn)$weight
    qgobjcn <- qgraph(as_adjacency_matrix(cn,attr = "weight"), DoNotPlot=TRUE )
    plot_out$contemporaneous$qgraph[[1]] <- qgobjcn
    plot_out$contemporaneous$igraph[[1]] <- cn
    plot(cn,
         vertex.label=V(cn)$name,
         vertex.label.color="black",curved=TRUE,
         vertex.size=12,
         edge.arrow.size=0.5,
         vertex.color="deepskyblue2",
         layout = qgobjcn$layout,
         edge.curved=0.15,
         main=ifelse(x$CALL$pars$model=="GVAR",
                     paste("Contemporaneous network from individual", x$CALL$pars$data_pre_fixed$ID[1]),
                     "Contemporaneous network"),...)
    if (communit){
      plot_out$contemporaneous$communities[[1]] <- cluster_spinglass(cn, stop.temp = 0.000000000001, cool.fact = 0.999,
                                                                     implementation = "neg", update.rule = "config", gamma = 0, gamma.minus = 0)
      plot(plot_out$contemporaneous$communities[[1]], cn, layout=qgobj$layout,
           main=ifelse(x$CALL$pars$model=="GVAR",
                       paste("Spinglass communities in the \ncontemporaneous network from individual", x$CALL$pars$data_pre_fixed$ID[1]),
                       "Spinglass communities in the \ncontemporaneous network"),...)
    }
    return(plot_out)

    
    

    
  }
  
  
  
  
}