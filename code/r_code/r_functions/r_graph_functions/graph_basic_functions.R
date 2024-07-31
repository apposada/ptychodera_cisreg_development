# Helper functions

#' Simple wrapper to calculate a number of metrics on all the vertex of a graph
vertex_metrics <- function(g){
  require(igraph)
  
  ##' degree
  message("degree")
  message("outdegree")
  V(g)$outdegree <- igraph::degree(g, V(g), mode = "out")
  message("indegree")
  V(g)$indegree  <- igraph::degree(g, V(g), mode = "in")
  V(g)$rel_outdegree <- V(g)$outdegree/(V(g)$indegree+V(g)$outdegree)
  
  ##' centrality
  message("centrality, this might take a while...")
  V(g)$centr <- closeness(g, V(g), mode = "all", normalized = TRUE) # , weights = E(g)$weight ??
  # message("outcentrality") # deactivated by default
  # V(g)$outcentr <- closeness(g, V(g), mode = "out")
  # message("incentrality")
  # V(g)$incentr <- closeness(g, V(g), mode = "in")
  
  # message("harmonic centrality") # deactivated for older versions of R and igraph
  # V(g)$outharm <- harmonic_centrality(g, V(g), mode = "out")
  
  ##' betweenness
  message("betweenness, this might take a while...")
  V(g)$betw <- betweenness(g, V(g), directed = TRUE)
  
  #' Tagging TF and target genes
  message("tagging TFs")
  V(g)$is_TF <- ifelse(V(g)$outdegree > 0, TRUE, FALSE)
  message("tagging tgs")
  V(g)$is_tg <- ifelse(V(g)$outdegree == 0, TRUE, FALSE)
  
  return(g)
  
}

#' Retrieve the intersection and the exclusions of two sets of things. Good for the TFs, Targets, etc.
OL <- function(a,b){
  y = 
    list(
      a = a[!(a %in% b)],
      a_b = a[a %in% b],
      b = b[!(b %in% a)]
    )
  return(y)
}

rel <- function(x){
  y = (x-min(x))/(max(x)-min(x))
  return(y)
}

#' Filter top nodes from an OL node list based on metrics such as indegree, outdegree, centrality, etc. ;  see function OL()
filter_top_nodes <- function(l, g_a, g_b, q = .9, metric = "indegree"){
  l_top <-l
  
  q_a <- quantile(vertex_attr(g_a, metric, V(g_a)$name %in% l[[1]]), q)
  q_b <- quantile(vertex_attr(g_b, metric, V(g_b)$name %in% l[[3]]), q)
  q_ab_a <- quantile(vertex_attr(g_a, metric, V(g_a)$name %in% l[[2]]), q)
  q_ab_b <- quantile(vertex_attr(g_b, metric, V(g_b)$name %in% l[[2]]), q)
  
  l_top[[1]] <-
    l_top[[1]][
      l_top[[1]] %in%
        V(g_a)$name[vertex_attr(g_a,metric) > q_a]
    ]
  l_top[[3]] <-
    l_top[[3]][
      l_top[[3]] %in%
        V(g_b)$name[vertex_attr(g_b,metric) > q_b]
    ]
  l_top[[2]] <-
    OL(
      l_top[[2]][
        l_top[[2]] %in% V(g_a)$name[vertex_attr(g_a,metric) > q_ab_a]
      ],
      l_top[[2]][
        l_top[[2]] %in% V(g_b)$name[vertex_attr(g_b,metric) > q_ab_b]
      ]
    )[[2]]
  names(l_top) <- paste0(names(l_top),"_top")
  
  return(l_top)
}

#' Creates an `eulerr::`-friendly vector to plot overlaps of two sets.
calc_overlaps <- function(x){
  
  a = length(x[[1]])
  b = length(x[[3]])
  a_and_b = length(x[[2]])
  
  v = c(a,b,a_and_b)
  n = names(x)[c(1,3)]
  
  res = setNames(v, c(n, paste(n, collapse = "&") ) )
  
  return(res)
}

#' We should have a function to turn the igraph object into and adjacency matrix 

#' Subgraph functions

#' Comb graphs for nice plotting
top_edges_per_tf <- function(g,tf,top=3, mode = "out"){
  
  if(mode == "out"){
    tf_edges <- E(g)[.from(tf)]
  } else if(mode == "in"){
    tf_edges <- E(g)[.to(tf)]
  } else {
    Stop("incorrect mode. Check input parameters")
  }
  
  if (length(tf_edges) > top){
    tf_edges_top <- rev(order( E(g)[tf_edges]$weight ))[ 1:top ]    
  } else {
    tf_edges_top <- rev(order( E(g)[tf_edges]$weight ))
  }
  
  edges_top <- E(g)[which(E(g) %in% tf_edges)][tf_edges_top]
  
  y <- which(E(g) %in% edges_top)
  
  return(y)
}

#' Quick wrapper to retrieve subgraphs keeping top interactions, see subgraph_with_top_edges_per_tf() and top_edges_per_tf()
make_subgraph <- function(g,v){
  require(igraph)
  g_ <- induced.subgraph(
    g,
    vids = v,
    impl = "create_from_scratch"
  )
  g_ <- subgraph_with_top_edges_per_tf(g_, mode = "both", top = 2)
  
  return(g_)
}

#' Create a subgraph keeping the top interactions stemming from every TF
subgraph_with_top_edges_per_tf <- function(g, top = 3, delete_isolated = TRUE, mode = "out"){
  
  top_edges <- integer()
  
  if(mode %in% c("out","in")){
    for(tf in V(g)$name[igraph::degree(g,V(g), mode = "out") > 0]){ # > 0 == TFs
      e <- top_edges_per_tf(g=g, tf = tf, top = top, mode = mode)
      top_edges <- c(top_edges,e)
    }
  } else if(mode == "both"){
    # in
    e_in_top <- integer()
    for(tf in V(g)$name[igraph::degree(g,V(g), mode = "out") > 0]){ # > 0 == TFs
      e <- top_edges_per_tf(g=g, tf = tf, top = top, mode = "in")
      e_in_top <- c(e_in_top,e)
    }
    # out
    e_out_top <- integer()
    for(tf in V(g)$name[igraph::degree(g,V(g), mode = "out") > 0]){ # > 0 == TFs
      e <- top_edges_per_tf(g=g, tf = tf, top = top, mode = "out")
      e_out_top <- c(e_out_top,e)
    }
    
    top_edges <- c(e_in_top, e_out_top)
  }
  
  g_ <- subgraph.edges(graph = g, eids = E(g)[top_edges], delete.vertices = delete_isolated)
  
  return(g_)
}

#' Create a subgraph keeping the top interactions arriving to every gene
subgraph_with_top_edges_per_tg <- function(g, top = 3, delete_isolated = TRUE, mode = "in"){
  
  top_edges <- integer()
  
  if(mode %in% c("out","in")){
    for(tf in V(g)$name[igraph::degree(g,V(g), mode = "in") > 0]){ # "in" > 0 == target genes (but also TFs if they receive interactions!)
      e <- top_edges_per_tf(g=g, tf = tf, top = top, mode = mode)
      top_edges <- c(top_edges,e)
    }
  } else if(mode == "both"){
    # in
    e_in_top <- integer()
    for(tf in V(g)$name[igraph::degree(g,V(g), mode = "in") > 0]){ # "in" > 0 == target genes (but also TFs if they receive interactions!)
      e <- top_edges_per_tf(g=g, tf = tf, top = top, mode = "in")
      e_in_top <- c(e_in_top,e)
    }
    # out
    e_out_top <- integer()
    for(tf in V(g)$name[igraph::degree(g,V(g), mode = "in") > 0]){ # "in" > 0 == target genes (but also TFs if they receive interactions!)
      e <- top_edges_per_tf(g=g, tf = tf, top = top, mode = "out")
      e_out_top <- c(e_out_top,e)
    }
    
    top_edges <- c(e_in_top, e_out_top)
  }
  
  g_ <- subgraph.edges(graph = g, eids = E(g)[top_edges], delete.vertices = delete_isolated)
  
  return(g_)
}

#' Parse and tidy the influence output
influ_table_2 <- function(x,tftable) {
  require(ggplotify)
  require(gridExtra)
  
  x <- merge(x,tftable,by.x=1,by.y=1,all.x=T)
  x$TFclass[is.na(x$TFclass)] <- " "
  x$col[is.na(x$col)] <- "gray"
  x <- x[rev(order(x$influence_score)),]
  rownames(x) <- NULL
  
  p <- as_ggplot( as.grob( function(){
    plot(
      main="Main factors",
      x$factor_fc, x$influence_score,
      pch = 21, bg = x$col, col = "black",
      xlab="log2fold change of TF",
      ylab="ANANSE influence score",
      bty="n",
      xlim= c(0,max(x$factor_fc) + 1)
    )
    text(
      x$factor_fc[1:20], x$influence_score[1:20]+0.01,
      paste(x$factor[1:20],x$TFclass[1:20], sep = " "), cex=0.7
    )
  } ) )
  
  res <- list(
    influence_table = x,
    influence_plot = p
  )
  
  return(res)
  
}
