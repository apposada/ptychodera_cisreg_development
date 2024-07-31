#' Wrappers for handling ANANSE outputs
require(igraph)
require(data.table)

#' Load network
LoadNetworkData <- function(x) {
  nw <- fread(
    text = gsub(
      "—", # replace dash line from ananse rows
      "\t", # with a tab, thus "tf" "\t" "tg"
      readLines(x) #taken from: 
    ),
    skip = 1 # skip the colnames
  )
  if(ncol(nw) == 3) colnames (nw) <- c("tf","tg","prob") # simple output `ananse network [...]`
  if(ncol(nw) == 7) colnames(nw) <- c("tf","tg","prob", "tf_expr", "tg_expr", "w_binding", "activity") # full output `ananse network -f [...]`

  # output of function: a table/dataframe in the structure of igraph
  return(nw)
}

#' Generate network
GenerateNetwork <- function(nw, q=0.85){
  q_filt <- nw$prob > quantile(nw$prob, q)
  ngenes <- length(
    unique(
      c(nw$tf,nw$target)
    )
  )
  
  nw_G <- graph.data.frame( # the initial igraph object, the 2 first columns are the "from" and "to" nodes; the rest of the columns are interpreted as edge attributes
    d = nw[q_filt,],
    directed = T
  )
  E(nw_G)$weight <- E(nw_G)$prob
  E(nw_G)$width <- E(nw_G)$prob
  
  # output of function: a filtered network in the shape of an igraph object
  return(nw_G)
}

#' Parse the network
ParseNetwork <- function(graph, list_attr = list(tflist,efflist,classANDcolor,td_hk,gfams,funcat,markerlist,...) ){
  # 1.1 asigna lo q tenga que asignar de datos para igraph (TF class, color..)
  df_attr <- data.frame(
    id = V(graph)$name,
    index = seq(
      1:length(V(graph))
    )
  )
  
  # merge data
  for (i in list_attr){
    df_attr <- merge(
      df_attr,
      i,
      by.x = 1,
      by.y = 1,
      all.x = T
    )  
  }
  
  # cleanup data
  df_attr[is.na(df_attr)] <- ""
  df_attr <- df_attr[order(df_attr$index),] #IMPORTANT
  df_attr <- df_attr[!duplicated(df_attr$id),]
  rownames(df_attr) <- NULL
  
  # slap together with the network
  for (i in colnames(df_attr)[-1]){
    graph <- set_vertex_attr(graph, paste0(i), value = df_attr[[i]])
  }
  #' output of the function: a graph with additional 
  #' attributes to the different nodes and edges
  output <- list(
    graph,
    df_attr
  )
  return(output)
}

#' Network Comparison
CompareNetworks <- function(
  g_a, g_b, influence = NULL,
  name_network_a = "a", name_network_b = "b",
  col_a = "darkorange", col_b = "purple",
  geneset_interest = NULL, q_targets = 0.9,
  tfs, gene_universe = NULL, id2go = NULL
){
  
  #' Common genes
  message("#### Common and Exclusive genes in graphs ####")
  ##' common and exclusive genes; no distinction
  message("Defining exclusive and common targets")
  OL_ab <- OL(a = V(g_a)$name, b = V(g_b)$name)
  names(OL_ab) <- c("genes_exclusive_a","genes_common_ab","genes_exclusive_b")
  
  ##' common and exclusive TFs
  message("Defining exclusive and common TFs")
  tf_ab <- OL( a = V(g_a)$name[V(g_a)$is_TF], b = V(g_b)$name[V(g_b)$is_TF] )
  names(tf_ab) <- c("TFs_exclusive_a","TFs_common_ab","TFs_exclusive_b") 
  
  ##' common and exclusive TGs
  message("Defining exclusive and common targets")
  tg_ab <- OL( a = V(g_a)$name[V(g_a)$is_tg], b = V(g_b)$name[V(g_b)$is_tg] )
  names(tg_ab) <- c("tgs_exclusive_a","tgs_common_ab","tgs_exclusive_b")
  
  message("Defining exclusive and common targets--top")
  tg_ab_top <-filter_top_nodes(l = tg_ab, g_a = g_a, g_b = g_b, q = q_targets, metric = "indegree")
  
  ##' common and exclusive trans-dev
  message("Defining exclusive and common trans-devs")
  td_ab <- OL( a = V(g_a)$name[V(g_a)$TDHK == "td"], b = V(g_b)$name[V(g_b)$TDHK == "td"] )
  names(td_ab) <- c("td_exclusive_a","td_common_ab","td_exclusive_b")
  
  #' Present the metrics in common dataframes
  message("#### Common Metrics ####")
  ##' Centrality
  message("Centrality of ALL common genes across networks")
  centr_ab <- 
    data.frame(
      a = V(g_a)$centr[match(OL_ab$genes_common_ab, V(g_a)$name)],
      b = V(g_b)$centr[match(OL_ab$genes_common_ab, V(g_b)$name)],
      row.names = OL_ab$genes_common_ab
    )
  message("Centrality of TFs across networks")
  centr_ab_tf <- centr_ab[rownames(centr_ab) %in% tf_ab$TFs_common_ab,]
  message("Centrality trans-dev genes across networks")
  centr_ab_td <- centr_ab[rownames(centr_ab) %in% V(g_a)$name[V(g_a)$TDHK == "td"],]
  
  
  ##' Betweenness
  message("Betweenness of ALL common genes across networks")
  betw_ab <- 
    data.frame(
      a = V(g_a)$betw[match(OL_ab$genes_common_ab, V(g_a)$name)],
      b = V(g_b)$betw[match(OL_ab$genes_common_ab, V(g_b)$name)],
      row.names = OL_ab$genes_common_ab
    )
  message("Betweenness of TFs across networks")
  betw_ab_tf <- betw_ab[rownames(betw_ab) %in% tf_ab$TFS_common_ab,]
  message("Centrality trans-dev genes across networks")
  betw_ab_td <- betw_ab[rownames(betw_ab) %in% V(g_a)$name[V(g_a)$TDHK == "td"],]
  
  ##' Relative out-degree
  message("Relative outdegree of ALL common genes across networks")
  reloutdg_ab <- 
    data.frame(
      a = V(g_a)$centr[match(OL_ab$genes_common_ab, V(g_a)$name)],
      b = V(g_b)$centr[match(OL_ab$genes_common_ab, V(g_b)$name)],
      row.names = OL_ab$genes_common_ab
    )
  message("Relative outdegree of TFs across networks")
  reloutdg_ab_tf <- reloutdg_ab[rownames(reloutdg_ab) %in% tf_ab$TFS_common_ab,]
  message("Relative outdegree of trans-dev genes across networks")
  reloutdg_ab_td <- reloutdg_ab[rownames(reloutdg_ab) %in% V(g_a)$name[V(g_a)$TDHK == "td"],]
  
  #' Sub-Graphs
  message("#### Subgraphs of common and exclusive genes ####")
  ##' Subgraph TFs exclusive to a
  g_tfs_a <- make_subgraph(g = g_a, v = tf_ab$TFs_exclusive_a )
  ##' Subgraph common TFs, more prominent in a
  g_tfs_ab_a <- make_subgraph(g = g_a, v = rownames(centr_ab_tf[rel(centr_ab_tf$b)/rel(centr_ab_tf$a) < .95,]) )
  ##' Subgraph common TFs, more prominent in b
  g_tfs_ab_b <- make_subgraph(g = g_b, v = rownames(centr_ab_tf[rel(centr_ab_tf$b)/rel(centr_ab_tf$a) > 1,]) )
  ##' Subgraph TFs exclusive to b
  g_tfs_b <- make_subgraph(g = g_b, v = tf_ab$TFs_exclusive_b )
  
  ##' Subgraph TFs exclusive to a
  g_td_a <- make_subgraph(g = g_a, v = td_ab$TD_exclusive_a )
  ##' Subgraph common TFs, more prominent in a
  g_td_ab_a <- make_subgraph(g = g_a, v = rownames(centr_ab_td[rel(centr_ab_td$b)/rel(centr_ab_td$a) < .95,]) )
  ##' Subgraph common TFs, more prominent in b
  g_td_ab_b <- make_subgraph(g = g_b, v = rownames(centr_ab_td[rel(centr_ab_td$b)/rel(centr_ab_td$a) > 1,]) )
  ##' Subgraph TFs exclusive to b
  g_td_b <- make_subgraph(g = g_b, v = td_ab$TD_exclusive_b )
  
  #' Influence
  message("#### Influence results ####")
  if (is.null(influence) == FALSE){
    inf_res = influ_table_2(x = influence, tftable = tfs)
    g_inf = make_subgraph(g = g_b, v = inf_res$influence_table$factor)
    # node size
    V(g_inf)$size = inf_res$influence_table$factor_fc[match(V(g_inf)$name, inf_res$influence_table$factor)]
    # node color
    V(g_inf)$col[V(g_inf)$col == ""] = "gray"
    V(g_inf)$color <- V(g_inf)$col
    # weight/width
    E(g_inf)$rel_w = rel(E(g_inf)$weight)
  } else {
    inf_res = "no Influence provided."
    g_inf = "no Influence provided."
  }
  
  #' Various metrics together in a DF
  message("#### Dataframe of comparative metrics ####")
  compare_df <- data.frame(
    metric = c(
      "Num Genes in Network",
      "Num Active TFs",
      "Mean Connections per TF",
      "Median Connections per TF",
      "Num Self-Regulated TFs",
      "Num Target Genes",
      "Mean Connections per TG",
      "Median Connections per TG",
      "Mean Emitted Connections",
      "Mean Received Connections",
      "Number of common TFs",
      "Number of exclusive TFs",
      "Number of common targets",
      "Number of exclusive targets"
    ),
    a = c(
      vcount(g_a),
      length(which(V(g_a)$is_TF)),
      mean(V(g_a)$outdegree[V(g_a)$is_TF]),
      median(V(g_a)$outdegree[V(g_a)$is_TF]),
      length(which(  tail_of(g_a,E(g_a)) == head_of(g_a,E(g_a))  )),
      length(which(V(g_a)$is_tg)),
      mean(V(g_a)$indegree[V(g_a)$is_tg]),
      median(V(g_a)$indegree[V(g_a)$is_tg]),
      mean(V(g_a)$outdegree),
      mean(V(g_a)$indegree),
      length(tf_ab$TFs_common_ab),
      length(tf_ab$TFs_exclusive_a),
      length(tg_ab$tgs_common_ab),
      length(tg_ab$tgs_exclusive_a)
    ),
    b = c(
      vcount(g_a),
      length(which(V(g_b)$is_TF)),
      mean(V(g_b)$outdegree[V(g_b)$is_TF]),
      median(V(g_b)$outdegree[V(g_b)$is_TF]),
      length(which(  tail_of(g_b,E(g_b)) == head_of(g_b,E(g_b))  )),
      length(which(V(g_b)$is_tg)),
      mean(V(g_b)$indegree[V(g_b)$is_tg]),
      median(V(g_b)$indegree[V(g_b)$is_tg]),
      mean(V(g_b)$outdegree),
      mean(V(g_b)$indegree),
      length(tf_ab$TFs_common_ab),
      length(tf_ab$TFs_exclusive_b),
      length(tg_ab$tgs_common_ab),
      length(tg_ab$tgs_exclusive_b)
    )
  )
  # Gene set of interest
  if (is.null(geneset_interest) == FALSE){	
    message("#### Genes of Interest ####")
    compare_df <-
      rbind(
        compare_df,
        data.frame(
          metric = "No. genes of interest in TGs",
          a = length(which(tg_ab$tgs_exclusive_a %in% geneset_interest)),
          b = length(which(tg_ab$tgs_exclusive_b %in% geneset_interest))
        )
      )
  }
  
  #' Centrality per TF class
  message("#### Centrality per TF class ####")
  message("Retrieving centrality per TF class; graph A")
  f = V(g_a)$TFclass != ""
  class_cent_a <-
    data.frame(
      id = V(g_a)$name[f],
      class = V(g_a)$TFclass[f],
      centr = V(g_a)$centr[f],
      betw = V(g_a)$betw[f]
    )
  message("Retrieving centrality per TF class; graph B")
  f = V(g_b)$TFclass != ""
  class_cent_b <-
    data.frame(
      id = V(g_b)$name[f],
      class = V(g_b)$TFclass[f],
      centr = V(g_b)$centr[f],
      betw = V(g_b)$betw[f]
    )
  
  #' GO analyses
  if( all(sapply(list(gene_universe,id2go),function(x) !(is.null(x)))) ) {
    message("#### Gene Ontologies ####")
    message("#### GOs of Transcription Factors ####")
    GOs_TFs <- getGOs(
      tf_ab,
      gene_universe = gene_universe,
      gene2GO = id2go,
      alg = "elim"
    )
    
    message("#### GOs of Target Genes (using only top 10% of common) ####")
    GOs_targets <- getGOs(
      tg_ab_top,
      gene_universe = gene_universe,
      gene2GO = id2go,
      alg = "elim"
    )
  } else {
    GOs_TFs <- "no GO available or requested."
    GOs_targets <- "no GO available or requested."
  }
  
  message("#### Retrieving results ####")
  res <-
    list(
      comparative = compare_df,
      gene_overlap = OL_ab,
      TF_overlap = tf_ab,
      tg_overlap = tg_ab,
      tg_overlap_top = tg_ab_top,
      transdev_overlap = td_ab,
      centrality_changes = centr_ab,
      centrality_changes_TFs = centr_ab_tf,
      centrality_changes_transdev = centr_ab_td,
      betweenness_changes = betw_ab,
      betweenness_changes_TFs = betw_ab_tf,
      betweenness_changes_transdev = betw_ab_td,
      rel_outdegree_changes = reloutdg_ab,
      rel_outdegree_changes_TFs = reloutdg_ab_tf,
      rel_outdegree_changes_transdev = reloutdg_ab_td,
      subgraph_TFs_exclusive_a = g_tfs_a,
      subgraph_common_TFs_a = g_tfs_ab_a,
      subgraph_common_TFs = g_tfs_ab_b,
      subgraph_TFs_exclusive_b = g_tfs_b,
      influence_results = 
        list(
          influence_table = inf_res[[1]],
          influence_plot = inf_res[[2]],
          influence_subgraph = g_inf
        ),
      centrality_TFclass_a = class_cent_a,
      centrality_TFclass_b = class_cent_b,
      GOs_TFs = GOs_TFs,
      GOs_targets	= GOs_targets
    )
  
  return(res)
}


# 3. funcion para plotear la centrality de los TFs over time (fuzz/broadlines depicting quantile or median or whatever)


#' 1.10 concordancia motivos atac homer y motivos de los tfs de la network (los + centrales, o los q tienen mas targets, etc.)



