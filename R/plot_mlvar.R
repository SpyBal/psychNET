plot_mlvar <- function(x, plot.type, individual,communit, ...){
  plot_out <- list()
  if (missing(individual)) individual <- NULL
  if (missing(communit)) communit <- FALSE
  if (missing(plot.type) | is.null(plot.type)) plot.type <- "temporal"
  type <- plot.type
  qgobj <- NULL
  if (is.null(individual)){
    if (type=="temporal"){
      return(lapply(1:x$CALL$pars$lag, function(i){
        plot_out<-list()
        plot_out$temporal$qgraph <- list()
        plot_out$temporal$igraph <- list()
        dn <- graph_from_adjacency_matrix(x$results$Dir_net[[i]], mode ="directed",diag=FALSE, weighted = TRUE,
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
             main=paste("Lag",i,"temporal network"),...)
        if (communit){
          plot_out$temporal$communities[[i]] <- cluster_spinglass(dn, stop.temp = 0.000000000001, cool.fact = 0.999,
                                                                  implementation = "neg", update.rule = "config", gamma = 0, gamma.minus = 0)
          plot(plot_out$temporal$communities[[i]], dn, layout=qgobj$layout,
               main=paste("Spinglass communities in\n lag",i,"temporal network"),...)
        }
        
        plot_out}))
    }
    if (type=="contemporaneous"){
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
           vertex.label.color="black",
           edge.curved=0.15,
           vertex.size=12,
           edge.arrow.size=0.5,
           vertex.color="deepskyblue2",
           layout = qgobjcn$layout,
           main="Contemporaneous network",... )
      if (communit){
        plot_out$contemporaneous$communities[[1]] <- cluster_spinglass(cn, stop.temp = 0.000000000001, cool.fact = 0.999,
                                                                       implementation = "neg", update.rule = "config", gamma = 0, gamma.minus = 0)
        plot(plot_out$contemporaneous$communities[[1]], cn, layout=qgobj$layout,
             main="Spinglass communities in the \ncontemporaneous network",...)
      }
      return(plot_out)
    }
    if (type=="between"){
      plot_out$between$igraph <-list()
      plot_out$between$qgraph <-list()
      bn <- graph_from_adjacency_matrix(x$fit$results$Omega_mu$pcor$mean, mode ="undirected",diag=FALSE, weighted = TRUE,
                                        add.colnames = NULL, add.rownames = NA)
      E(bn)$color <- ifelse(E(bn)$weight > 0,'green','red')
      E(bn)$weight <- E(bn)$weight
      qgobjbn <- qgraph(as_adjacency_matrix(bn,attr = "weight"), DoNotPlot=TRUE )
      plot_out$between$qgraph[[1]] <- qgobjbn
      plot_out$between$igraph[[1]] <- bn
      plot(bn,
           vertex.label=V(bn)$name,
           vertex.label.color="black",
           edge.curved=0.15,
           vertex.size=12,
           vertex.color="deepskyblue2",
           layout = qgobjbn$layout,
           main= "Between subjects network",...)
      if (communit){
        plot_out$between$communities[[1]] <- cluster_spinglass(bn, stop.temp = 0.000000000001, cool.fact = 0.999,
                                                                       implementation = "neg", update.rule = "config", gamma = 0, gamma.minus = 0)
        plot(plot_out$between$communities[[1]], bn, layout=qgobj$layout,
             main="Spinglass communities in the \nbetween subjects network",...)
      }
      
      return(plot_out)
    }
    if (type == "both"){
      
      list1 <- lapply(1:x$CALL$pars$lag, function(i){
        plot_out<-list()
        plot_out$temporal$qgraph <- list()
        plot_out$temporal$igraph <- list()
        dn <- graph_from_adjacency_matrix(x$results$Dir_net[[i]], mode ="directed",diag=FALSE, weighted = TRUE,
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
             main=paste("Lag",i,"temporal network"))
        if (communit){
          plot_out$temporal$communities[[i]] <- cluster_spinglass(dn, stop.temp = 0.000000000001, cool.fact = 0.999,
                                                                  implementation = "neg", update.rule = "config", gamma = 0, gamma.minus = 0)
          plot(plot_out$temporal$communities[[i]], dn, layout=qgobj$layout,
               main=paste("Spinglass communities in\n lag",i,"temporal network"),...)
        }
        plot_out})
      
      
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
           vertex.label.color="black",
           edge.curved=0.15,
           vertex.size=12,
           edge.arrow.size=0.5,
           vertex.color="deepskyblue2",
           layout = qgobjcn$layout,
           main="Contemporaneous network")
      
      plot_out$between$igraph <-list()
      plot_out$between$qgraph <-list()
      bn <- graph_from_adjacency_matrix(x$fit$results$Omega_mu$pcor$mean, mode ="undirected",diag=FALSE, weighted = TRUE,
                                        add.colnames = NULL, add.rownames = NA)
      E(bn)$color <- ifelse(E(bn)$weight > 0,'green','red')
      E(bn)$weight <- E(bn)$weight
      qgobjbn <- qgraph(as_adjacency_matrix(bn,attr = "weight"), DoNotPlot=TRUE )
      plot_out$between$qgraph[[1]] <- qgobjbn
      plot_out$between$igraph[[1]] <- bn
      print(plot(bn,
           vertex.label=V(bn)$name,
           vertex.label.color="black",
           edge.curved=0.15,
           vertex.size=12,
           vertex.color="deepskyblue2",
           layout = qgobjbn$layout,
           main="Between subjects network"))
      
      
      if (communit){
        plot_out$between$communities[[1]] <- cluster_spinglass(bn, stop.temp = 0.000000000001, cool.fact = 0.999,
                                                               implementation = "neg", update.rule = "config", gamma = 0, gamma.minus = 0)
        plot_out$contemporaneous$communities[[1]] <- cluster_spinglass(cn, stop.temp = 0.000000000001, cool.fact = 0.999,
                                                                       implementation = "neg", update.rule = "config", gamma = 0, gamma.minus = 0)
        plot(plot_out$contemporaneous$communities[[1]], cn, layout=qgobj$layout,
             main="Spinglass communities in the \ncontemporaneous network",...)
        plot(plot_out$between$communities[[1]], bn, layout=qgobj$layout,
             main="Spinglass communities in the \nbetween subjects network",...)
      }
      
      return(c(list1,plot_out))
    }
    
  }else{
    
    if (x$fit$input$temporal != "fixed"){
      indexes <- individual
      if (type=="temporal"){
        lapply(indexes, function(j){
          plot_out<-list()
          plot_out$temporal$qgraph <- list()
          plot_out$temporal$igraph <- list()
          result_temp <- x$results$Subjects$Dir_net[[j]]
          lapply(1:x$CALL$pars$lag, function(i){
            dn <- graph_from_adjacency_matrix(result_temp[[i]], mode ="directed",diag=FALSE, weighted = TRUE,
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
                 main=paste("Lag",i,"temporal network from individual",j),...)
            if (communit){
              plot_out$temporal$communities[[i]] <- cluster_spinglass(dn, stop.temp = 0.000000000001, cool.fact = 0.999,
                                                                      implementation = "neg", update.rule = "config", gamma = 0, gamma.minus = 0)
              plot(plot_out$temporal$communities[[i]], dn, layout=qgobj$layout,
                   main=paste("Spinglass communities in\n Lag",i,"temporal network from individual",j),...)
            }
            
            plot_out})})
      }
      if (type=="contemporaneous"){
        plot_out <- list()
        plot_out$Contemporaneous$qgraph <- list()
        plot_out$Contemporaneous$igraph <- list()
        lapply(indexes, function(j){
          result_temp <- x$results$Subjects$UnDir_net[[j]]

          cn <- graph_from_adjacency_matrix(result_temp, mode ="undirected",diag=FALSE, weighted = TRUE,
                                              add.colnames = NULL, add.rownames = NA)
          E(cn)$color <- ifelse(E(cn)$weight > 0,'green','red')
          E(cn)$weight <- E(cn)$weight
          qgobj <- qgraph(as_adjacency_matrix(cn,attr = "weight"), DoNotPlot=TRUE )
          plot_out$Contemporaneous$qgraph[[j]] <- qgobj
          plot_out$Contemporaneous$igraph[[j]] <- cn
          plot(cn,
               vertex.label=V(cn)$name,
               vertex.label.color="black",
               edge.curved=0.15,
               vertex.size=12,
               edge.arrow.size=0.5,
               vertex.color="deepskyblue2",
               layout = qgobj$layout,
               main=paste("Contemporaneous network from individual",j),...)
          if (communit){
            plot_out$contemporaneous$communities[[j]] <- cluster_spinglass(cn, stop.temp = 0.000000000001, cool.fact = 0.999,
                                                                           implementation = "neg", update.rule = "config", gamma = 0, gamma.minus = 0)
            plot(plot_out$contemporaneous$communities[[j]], cn, layout=qgobj$layout,
                 main=paste("Spinglass communities in the \ncontemporaneous network from individual",j),...)
          }
          
          plot_out})
      }
      if (type=="between"){
      stop("Between subjects network is not available for individuals")
      }
      
      if (type=="both"){
        lapply(indexes, function(j){
          plot_out <- list()
          plot_out$temporal$qgraph <- list()
          plot_out$temporal$igraph <- list()
          plot_out$Contemporaneous$qgraph <- list()
          plot_out$Contemporaneous$igraph <- list()
          plot_out$temporal$communities <- list()
          plot_out$Contemporaneous$communities <- list()
          
          result_temp <- x$results$Subjects$Dir_net[[j]]
          lapply(1:x$CALL$pars$lag, function(i){
            dn <- graph_from_adjacency_matrix(result_temp[[i]], mode ="directed",diag=FALSE, weighted = TRUE,
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
                main=paste("Lag",i,"temporal network from individual",j),...)
            if (communit){
              plot_out$temporal$communities[[i]] <- cluster_spinglass(dn, stop.temp = 0.000000000001, cool.fact = 0.999,
                                                                      implementation = "neg", update.rule = "config", gamma = 0, gamma.minus = 0)
              plot(plot_out$temporal$communities[[i]], dn, layout=qgobj$layout,
                   main=paste("Spinglass communities in\n Lag",i,"temporal network from individual",j),...)
            }})
          
          result_temp2 <- x$results$Subjects$UnDir_net[[j]]
          cn <- graph_from_adjacency_matrix(result_temp2, mode ="undirected",diag=FALSE, weighted = TRUE,
                                            add.colnames = NULL, add.rownames = NA)
          E(cn)$color <- ifelse(E(cn)$weight > 0,'green','red')
          E(cn)$weight <- E(cn)$weight
          qgobj <- qgraph(as_adjacency_matrix(cn,attr = "weight"), DoNotPlot=TRUE )
          plot_out$Contemporaneous$qgraph[[1]] <- qgobj
          plot_out$Contemporaneous$igraph[[1]] <- cn
          plot(cn,
               vertex.label=V(cn)$name,
               vertex.label.color="black",
               edge.curved=0.15,
               vertex.size=12,
               edge.arrow.size=0.5,
               vertex.color="deepskyblue2",
               layout = qgobj$layout,
               main=paste("Contemporaneous network from individual",j),...)
          if (communit){
            plot_out$contemporaneous$communities[[j]] <- cluster_spinglass(cn, stop.temp = 0.000000000001, cool.fact = 0.999,
                                                                           implementation = "neg", update.rule = "config", gamma = 0, gamma.minus = 0)
            plot(plot_out$contemporaneous$communities[[j]], cn, layout=qgobj$layout,
                 main=paste("Spinglass communities in the \ncontemporaneous network from individual",j),...)
          }
          plot_out
          })
      }
      
    }else{
      stop("Individual networks are not available when temporal='fixed'.")
    }
  }
}