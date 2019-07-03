plot_dfm <- function(x, plot.type, communit, ...){
  plot_out <- list()
  if (missing(plot.type)) plot.type <- "temporal"
  if (missing(communit)) communit <- FALSE
  
  type <- plot.type
  if (x$CALL$pars$lag ==1){
    if (type=="temporal"){
      dn <- graph_from_adjacency_matrix(x$results$Dir_net[[1]], mode ="directed",diag=TRUE, weighted = TRUE,
                                        add.colnames = NULL, add.rownames = NA)
      E(dn)$color <- ifelse(E(dn)$weight > 0,'green','red')
      E(dn)$weight <- E(dn)$weight
      plot_out$temporal$igraph <-list()
      plot_out$temporal$qgraph <-list()
      
      qgobj <- qgraph(as_adjacency_matrix(dn,attr = "weight"), DoNotPlot=TRUE )
      plot_out$temporal$qgraph[[1]] <- qgobj
      plot_out$temporal$igraph[[1]] <- dn
      
      
      plot(dn,
           vertex.label=V(dn)$name,
           vertex.label.color="black",curved=TRUE,
           vertex.size=12,
           edge.arrow.size=0.5,
           vertex.color="deepskyblue2",
           layout = qgobj$layout,
           edge.curved=0.15,
           main=paste("Temporal network from individual", x$CALL$pars$data_pre_fixed$ID[1]),...)
      
      if (communit){
        plot_out$temporal$communities[[1]] <- cluster_spinglass(dn, stop.temp = 0.000000000001, cool.fact = 0.999,
                                                                implementation = "neg", update.rule = "config", gamma = 0, gamma.minus = 0)
        plot(plot_out$temporal$communities[[1]], dn, layout=qgobj$layout,
             main=paste("Spinglass communities \n in the temporal network from individual", x$CALL$pars$data_pre_fixed$ID[1]),...)
      }
      return(plot_out)
      
      
    }
    if (type=="contemporaneous"){
      cn <- graph_from_adjacency_matrix(x$results$UnDir_net, mode ="undirected",diag=FALSE, weighted = TRUE,
                                        add.colnames = NULL, add.rownames = NA)
      E(cn)$color <- ifelse(E(cn)$weight > 0,'green','red')
      E(cn)$weight <- E(cn)$weight
      plot_out$contemporaneous$igraph <-list()
      plot_out$contemporaneous$qgraph <-list()
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
           main=paste("Contemporaneous network from individual", x$CALL$pars$data_pre_fixed$ID[1]),...)
      if (communit){
        plot_out$contemporaneous$communities[[1]] <- cluster_spinglass(cn, stop.temp = 0.000000000001, cool.fact = 0.999,
                                                                implementation = "neg", update.rule = "config", gamma = 0, gamma.minus = 0)
        plot(plot_out$contemporaneous$communities[[1]], cn, layout=qgobj$layout,
             main=paste("Spinglass communities \n in the contemporaneous network from individual", x$CALL$pars$data_pre_fixed$ID[1]),...)
      }
      return(plot_out)
      
    }
    if (type=="between"){
      stop("Between subjects network is available only when model= 'MLVAR")
      
    }
    if (type=="both"){
      dn <- graph_from_adjacency_matrix(x$results$Dir_net, mode ="directed",diag=FALSE, weighted = TRUE,
                                        add.colnames = NULL, add.rownames = NA)
      E(dn)$color <- ifelse(E(dn)$weight > 0,'green','red')
      E(dn)$weight <- E(dn)$weight
      qgobj <- qgraph(as_adjacency_matrix(dn,attr = "weight"), DoNotPlot=TRUE )
      plot_out$temporal$qgraph[[1]] <- qgobj
      plot_out$temporal$igraph[[1]] <- dn
      plot(dn,
           vertex.label=V(dn)$name,
           vertex.label.color="black",curved=TRUE,
           vertex.size=12,
           edge.arrow.size=0.5,
           vertex.color="deepskyblue2",
           layout = qgobj$layout,
           edge.curved=0.15,
           main=paste("Temporal network from individual", x$CALL$pars$data_pre_fixed$ID[1]))
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
           main=paste("Contemporaneous network from individual", x$CALL$pars$data_pre_fixed$ID[1]))
      if (communit){
        plot_out$contemporaneous$communities[[1]] <- cluster_spinglass(cn, stop.temp = 0.000000000001, cool.fact = 0.999,
                                                                       implementation = "neg", update.rule = "config", gamma = 0, gamma.minus = 0)
        plot_out$temporal$communities[[1]] <- cluster_spinglass(dn, stop.temp = 0.000000000001, cool.fact = 0.999,
                                                                implementation = "neg", update.rule = "config", gamma = 0, gamma.minus = 0)
        plot(plot_out$temporal$communities[[1]], dn, 
             layout=qgobj$layout,
             main=paste("Spinglass communities \n in the dynamic network from individual", x$CALL$pars$data_pre_fixed$ID[1]),...)
        plot(plot_out$contemporaneous$communities[[1]], 
             cn, 
             layout=qgobj$layout,
             main=paste("Spinglass communities \n in the contemporaneous network from individual", x$CALL$pars$data_pre_fixed$ID[1]),...)
      }
      return(plot_out)
      
    }
  }else{
    if (type=="temporal"){
      plot_out<-list()
      plot_out$temporal$qgraph <- list()
      plot_out$temporal$igraph <- list()
      igraphsDIR <- lapply(1:x$CALL$pars$lag , function(i){
        dn <- graph_from_adjacency_matrix(x$results$Dir_net_fac[[i]], mode ="directed",diag=FALSE, weighted = TRUE,
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
             main=paste("Temporal network of the factors from individual", x$CALL$pars$data_pre_fixed$ID[1]),...)
        plot_out})
      if (communit) message("No communities are calculated for the factor process")
      return(igraphsDIR)
    }
    if (type=="contemporaneous"){
      stop("Contemporaneous network is only available for one lag DFM ")
    }
    if (type=="between"){
      stop("Between subjects network is available only when model= 'MLVAR")
    }
    if (type=="both"){
      igraphsDIR <- lapply(1:x$CALL$pars$lag, function(i){
        plot_out<-list()
        plot_out$temporal$qgraph <- list()
        plot_out$temporal$igraph <- list()
        dn <- graph_from_adjacency_matrix(x$results$Dir_net_fac[[i]], mode ="directed",diag=FALSE, weighted = TRUE,
                                          add.colnames = NULL, add.rownames = NA)
        E(dn)$color <- ifelse(E(dn)$weight > 0,'green','red')
        E(dn)$weight <- E(dn)$weight
        qgobj <- qgraph(as_adjacency_matrix(dn,attr = "weight"), DoNotPlot=TRUE )
        plot_out$temporal$qgraph[[i]] <- qgobj
        plot_out$temporal$igraph[[i]] <- dn
        plot(dn,
             vertex.label=V(dn)$name,
             vertex.label.color="black",
             edge.curved=0.15,
             vertex.size=12,
             edge.arrow.size=0.5,
             vertex.color="deepskyblue2",
             layout = qgobj$layout,
             main=paste(paste("Lag",i,"temporal network of the factors from individual",sep = ""), x$CALL$pars$data_pre_fixed$ID[1]),...)
        plot_out})
      if (communit) message("No communities are calculated for the factor process")
      message("Contemporaneous network not available when model= 'DFM'.")
      return(igraphsDIR)
    }
    
  }
  
}