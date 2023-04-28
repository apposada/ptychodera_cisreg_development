#' Functions ANANSE
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
  colnames (nw) <- c("tf","tg","prob")
  # output of function: a table/dataframe in the structure of igraph
  return(nw)
}

#' Filter network
FilterNetwork <- function(nw,q=0.85){
  q_filt <- nw$prob > quantile(nw$prob, q)
  ngenes <- length(
    unique(
      c(nw$tf,nw$tg)
    )
  )
  nw_filt <- nw[q_filt,]
  ngenes_filt <- length(
    unique(
      c(nw_filt$tf,nw_filt$tg)
    )
  )
  print(
    paste0(
      "Original network is ",
      ngenes,
      " genes. ",
      "Filtered network is ",
      ngenes_filt,
      " genes."
    )
  )
  # output of function: a filtered network in the shape of a data frame
  return(nw_filt)
}

#' Generate network
GenerateNetwork <- function(nw,q=0.85){
  q_filt <- nw$prob > quantile(nw$prob, q)
  ngenes <- length(
    unique(
      c(nw$tf,nw$target)
    )
  )
  nw_G <- graph.data.frame( # the initial igraph object
    d = nw[q_filt,],
    directed = T
  )
  E(nw_G)$width <- nw[q_filt,]$prob
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

#' Network stats and metrics
NetworkStats <- function(nw,graph,att, N = 10, C = 1,gene_universe,id2go){
  
  print("Basic Stats")
  numgenes_network <- length(V(graph))
  # top N connected TFs measured as # of stemming edges
  connectionxtf <- c(table(nw$tf))
  numtfs <- length(connectionxtf)
  meanconnectionxtf <- mean(connectionxtf)
  # topEmmitters <- names(rev(sort(table(nw$tf)))[1:N])
  
  # top N connected TGs measured as # of incoming edges
  connectionxtg <- c(table(nw$tg))
  numtgs <- length(connectionxtg)
  meanconnectionxtg <- mean(connectionxtg)
  # topReceivers <- names(rev(sort(table(nw$tg)))[1:N])
  
  # Num outgoing vs incoming interactions
  nw_InOut <- computeInVsOut(nw)
  topEmitters <- nw_InOut$id[rev(order(nw_InOut$emit))][1:N]
  topReceivers <- nw_InOut$id[rev(order(nw_InOut$recv))][1:N]
  
  # Num TFs that self-regulate in the network
  selfreg <- length(which(nw$tf == nw$tg))
  
  # Divide strongest and weakest components
  print("Components")
  graph_components <- data.frame(
    id=names(components(graph,mode=c("strong"))$membership),
    member=components(graph,mode=c("strong"))$membership
  )
  
  # new subgraph with only genes from largest CC
  mainCC <- names(
    rev(sort(table(
      graph_components$member
    )))[C])
  print(mainCC)
  
  graph_strong =
    induced_subgraph(
      graph,
      graph_components$id[graph_components$member == mainCC],
      impl = "auto"
    )
  
  ccsizes <- table(graph_components$member)[1:10]
  mainccsize <- max(table(graph_components$member))
  
  mainCC_genes <- graph_components$id[graph_components$member == mainCC]
  
  print("Centrality and category-based metrics")
  #' number TF genes
  numTFs <- length(V(graph_strong)[V(graph_strong)$TFEG == "TF"])
  #' number effector genes
  numEGs <- length(V(graph_strong)[V(graph_strong)$TFEG == "EG"])
  #' ratio TF/effectors in the network
  ratioTFEG <- ifelse(numEGs == 0, 0, numTFs/numEGs) #workaround to remove the error of having infinite TFEG ratio that messes with the plotting
  
  print("Centrality")
  #' centrality of all TFs --> top central TFs
  print("Centrality of TFs")
  tf_central <- getcentral(graph, att, "TFclass")
  centralTFs <-  rev(sort(unlist(tf_central)))[1:N]
  
  #' numero y centrality de transdevs en las networks,
  print("Centrality of Trans-Devs")
  td_central <- getcentral(graph, att, "TDHK")
  centralTDs <- rev(sort(unlist(td_central$td)))[1:N]
  
  #' number incoming edges per funcat family
  print("Edges of Functional Categories")
  funcat_edges <- merge(
    nw_InOut,
    att[,c("id","funcat")],
    by.x = 1,
    by.y = 1,
    all.x = T
  )
  #' centrality per funcat family
  print("Centrality of Functional Categories")
  funcat_central <- getcentral(graph,att,"funcat")
  
  #' GOs
  print("Getting Gene Ontology")
  GO_outputs <- getGOs( # re-define how these categories are defined
    list(
      tfs = V(graph_strong)$name[V(graph_strong)$name %in% att$id[att$TFclass != ""] ],
      active_tfs = unique(nw$tf),
      egs = V(graph_strong)$name[!(V(graph_strong)$name %in% att$id[att$TFclass != ""])],
      target_genes = unique(nw$tg),
      whole_graph = V(graph)$name
      ),
    gene_universe = gene_universe,
    gene2GO = id2go,
    alg = "elim"
    )
  
  print("Generating output")
  report <- list(
    numgenes_network = numgenes_network,
    connection_per_TF = connectionxtf,
    num_active_TFs = numtfs,
    mean_connection_per_TF = meanconnectionxtf,
    connection_per_TG = connectionxtg,
    num_TGs = numtgs,
    mean_connection_per_TG = meanconnectionxtg,
    top_Emitters = topEmitters,
    In_Out_per_Gene = nw_InOut,
    top_Receivers = topReceivers,
    num_selfregulated_TFs = selfreg,
    genes_per_component = graph_components,
    Connected_component_sizes = ccsizes,
    main_Connected_Component = mainCC,
    main_Connected_Component_graph = graph_strong,
    main_Connected_Component_size = mainccsize,
    main_Connected_Component_genes = mainCC_genes,
    num_total_TFs = numTFs,
    num_Effector_Genes = numEGs,
    ratio_TF_per_EG = ratioTFEG,
    Centrality_per_TFclass = tf_central,
    most_central_TFs = centralTFs,
    Centrality_per_TransDev_HK = td_central,
    most_central_transdevs = centralTDs,
    edges_per_funcat_family = funcat_edges,
    Centrality_per_funcat_family = funcat_central,
    GO_tables = GO_outputs$GOtable
  )
  return(report)
  print("Done.")
}

#' Generate tables with the stat reports
digestStats <- function(stats){
  # 1.9 genera tablas con reports de varios tipos
  #
}

#' Plot the basic metrics of the networks
NetworkPlots <- function(stats, nw, tfcol, layout = TRUE, colpal = default_colpal, pdf = FALSE, pdfname = NULL, ...){
  if (layout == TRUE) {
    layout(
      matrix(
        c(1,2,3,4,5,6,7,8,9,10,10,11,12,12,13),
        nrow = 5,
        byrow = TRUE
      )
    )
  }
  
  density_network(nw)
  # plot_numgenes_network(stats)
  plot_typegenes_network(stats)
  plot_connections(stats)
  # plot_connections_per_gene(stats)
  # plot_genes_ratioconnections(stats)
  plot_ratio_connections(stats)
  # plot_selfreg(stats)
  # plot_size_CC(stats)
  plot_centrality_TDHK(stats)
  plot_centrality_TFclass(stats$Centrality_per_TFclass,  f = tfcol, ylim = c(0.00005,0.00012) )
  plot_top_central_tfs(stats)
  plot_behavior_per_category2(stats)
  plot_top_central_transdev(stats)
  # plot_GOs(stats$GO_plots)
  par(mfrow=c(1,1))
}


#' Network Comparison
compareNetworks <- function(
  nw_a,nw_b,graph_a,graph_b,stats_a,stats_b,
  influence = NULL, name_network_a = "a", name_network_b = "b", col_a = "darkorange",
  col_b = "purple", geneset_interest = NULL, 
  top = 0.9, tfs, gene_universe, id2go
){  # requires nw, graph, graph strong, InVsOut, stats
  
  require(VennDiagram)
  
  #common_genes
  print("Defining exclusive and common targets")
  
  a_top_tgt <- names(table(nw_a$tg)[table(nw_a$tg) >= quantile(table(nw_a$tg),top)])
  b_top_tgt <- names(table(nw_b$tg)[table(nw_b$tg) >= quantile(table(nw_b$tg),top)])
  
  list_targets <- list(
    targets_exclusive_a = a_top_tgt[!(a_top_tgt %in% b_top_tgt)],
    targets_exclusive_b = b_top_tgt[!(b_top_tgt %in% a_top_tgt)],
    targets_common_ab = a_top_tgt[a_top_tgt %in% b_top_tgt]
  )
  
  numtgts_exclusive_a <- length(list_targets[[1]])
  numtgts_exclusive_b <- length(list_targets[[2]])
  numtgts_common <- length(list_targets[[3]])
  
  venn <- venn.diagram(
    x = list(a_top_tgt, b_top_tgt),
    category.names = c(name_network_a, name_network_b),
    filename = NULL,
    cex = 1,
    cat.cex = 1,
    lwd = 0.5,
    fill=c(col_a,col_b)
  )
  
  # GOs
  print("getting GO terms")
  GOs_targets <- getGOs(
    list_targets,
    gene_universe = gene_universe,
    gene2GO = id2go,
    alg = "elim"
  )
  
  # Putting data together
  print("Dataframe of various metrics")
  compare_df <- data.frame(
    metric = c(
      "Num Genes in Network",
      "Num Active TFs",
      "Mean Connections per TF",
      "Num Self-Regulated TFs",
      "Mean Connections per TG",
      "Median Connections per TG",
      "Mean Emitted Connections",
      "Mean Received Connections",
      "Number of common targets",
      "Number of exclusive targets"
    ),
    a = c(
      stats_a$numgenes_network,
      stats_a$num_active_TFs,
      stats_a$mean_connection_per_TF,
      stats_a$num_selfregulated_TFs,
      stats_a$mean_connection_per_TG,
      median(stats_a$connection_per_TG),
      mean(stats_a$In_Out_per_Gene$emit),
      mean(stats_a$In_Out_per_Gene$recv),
      numtgts_common,
      numtgts_exclusive_a
    ),
    b = c(
      stats_b$numgenes_network,
      stats_b$num_active_TFs,
      stats_b$mean_connection_per_TF,
      stats_b$num_selfregulated_TFs,
      stats_b$mean_connection_per_TG,
      median(stats_b$connection_per_TG),
      mean(stats_b$In_Out_per_Gene$emit),
      mean(stats_b$In_Out_per_Gene$recv),
      numtgts_common,
      numtgts_exclusive_b
    )
  )
  
  # merge InVsOut
  print("Gene Behaviour across networks")
  inout_ab <- merge(
    stats_a$In_Out_per_Gene,
    stats_b$In_Out_per_Gene,
    by = 1,
    all = TRUE
  )
  colnames(inout_ab) <- c("id","emit.a","recv.a","ratio.a","emit.b","recv.b","ratio.b")
  # centrality
  print("Gene centrality of ALL common genes across networks")
  common_genes_in_nw <-
    V(graph_a)$name [
      V(graph_a)$name %in% V(graph_b)$name
    ]
  
  centrality_per_network <- data.frame(
    id = common_genes_in_nw,
    centrality_a = closeness(
      graph_a,
      vids = common_genes_in_nw,
      mode = c("all"),
      weights = NULL,
      normalized = F
    ),
    centrality_b = closeness(
      graph_b,
      vids = common_genes_in_nw,
      mode = c("all"),
      weights = NULL,
      normalized = F
    )
  )
  
  # merge tf central
  print("Gene centrality of TFs across networks")
  tfcntr_ab <- merge(
    data.frame(
      id = names(unlist(stats_a$Centrality_per_TFclass)),
      a = unlist(stats_a$Centrality_per_TFclass)
    ),
    data.frame(
      id = names(unlist(stats_b$Centrality_per_TFclass)),
      b = unlist(stats_b$Centrality_per_TFclass)
    ),
    by = 1
  )
  tfcntr_ab$class <- sub("\\.TCONS..*","",tfcntr_ab$id)
  tfcntr_ab$id <- sub("..*TCONS","TCONS",tfcntr_ab$id)
  
  # merge transdev central
  print("Gene centrality of trans-dev genes across networks")
  td_ab <- merge(
    data.frame(
      id = names(unlist(stats_a$Centrality_per_TransDev_HK$td)),
      a = unlist(stats_a$Centrality_per_TransDev_HK$td)
    ),
    data.frame(
      id = names(unlist(stats_b$Centrality_per_TransDev_HK$td)),
      b = unlist(stats_b$Centrality_per_TransDev_HK$td)
    ),
    by = 1
  )
  
  # Generate graph of TFs
  print("Graph generation - TF genes")
  graph_tf_a <- induced.subgraph(
    graph_a,
    vids = tfcntr_ab$id[
      tfcntr_ab$b/tfcntr_ab$a < 
        0.95
    ],
    impl = "create_from_scratch"
  )
  V(graph_tf_a)$size <- 
    1 /
    (tfcntr_ab$b/tfcntr_ab$a)[
      match(
        V(graph_tf_a)$name,
        tfcntr_ab$id
      )
    ]
  
  graph_tf_b <- induced.subgraph(
    graph_b,
    vids = tfcntr_ab$id[
      tfcntr_ab$b/tfcntr_ab$a > 
        1
    ],
    impl = "create_from_scratch"
  )
  V(graph_tf_b)$size <- 
    (tfcntr_ab$b/tfcntr_ab$a)[
      match(
        V(graph_tf_b)$name,
        tfcntr_ab$id
      )
    ]
  
  print("Graph generation - Trans-dev genes")
  graph_td_a <- induced.subgraph(
    graph_a,
    vids = td_ab$id[td_ab$b/td_ab$a < 0.95 ],
    impl = "create_from_scratch"
  )
  V(graph_td_a)$size <- 
    1 /
    (td_ab$b/td_ab$a)[
      match(
        V(graph_td_a)$name,
        td_ab$id
      )
    ]
  
  graph_td_b <- induced.subgraph(
    graph_b,
    vids = td_ab$id[td_ab$b/td_ab$a > 1],
    impl = "create_from_scratch"
  )
  V(graph_td_b)$size <- 
    (td_ab$b/td_ab$a)[
      match(
        V(graph_td_b)$name,
        td_ab$id
      )
    ]
  
  # plots
  print("Plots")
  
  print("Mean connections per TF gene")
  barplot(
    main = "Mean connections per TF gene",
    height = c(compare_df[3,2],compare_df[3,2]),
    names = c(name_network_a, name_network_b),
    col = c(col_a, col_b)
  )
  
  print("Number of self-regulating TFs")
  barplot(
    main = "Number of self-regulating TFs",
    height = c(
      stats_a$num_selfregulated_TFs,
      stats_b$num_selfregulated_TFs),
    names.arg = c(name_network_a, name_network_b),
    col = c(col_a,col_b),
    sub = paste0("Chi-Sq pval ", chisq.test( #check this
      x = matrix(
        c(
          stats_a$num_selfregulated_TFs, 
          stats_a$num_active_TFs - stats_a$num_selfregulated_TFs,
          stats_b$num_selfregulated_TFs,
          stats_b$num_active_TFs - stats_b$num_selfregulated_TFs
        ),
        nrow=2)
    )$p.value)
  )
  
  print("connections per target gene")
  boxplot(
    main = "connections per target gene",
    list(
      c(table(nw_a$tg)),
      c(table(nw_b$tg))
    ),
    names = c(name_network_a, name_network_b),
    col = c(col_a,col_b),
    sub = paste0(
      "Wilcox p.value ",
      wilcox.test(
        x = c(table(nw_a$tg)),
        y = c(table(nw_b$tg))
      )$p.value
    )
  )
  
  # venn
  print("Venn Diagram Plot")
  plot(0,type="n", 
       xlab = "",
       ylab = "",
       main = "Targets shared between networks",
       axes = F
  )
  grid::grid.draw(venn)
  
  # Genes of interest
  
  if ( is.null(geneset_interest) == FALSE ) { 
    print("Number of genes of interest in targets")
    genes_interest_innetwork = c(
      name_network_a = length(which(list_targets$targets_exclusive_a %in% geneset_interest)),
      name_network_b = length(which(list_targets$targets_exclusive_b %in% geneset_interest))
    )
    barplot(
      main="Genes of interest\nin network",
      height=c(
        genes_interest_innetwork[1]/length(list_targets$targets_exclusive_a)*100,
        genes_interest_innetwork[2]/length(list_targets$targets_exclusive_b)*100
      ),
      col = c(
        col_a,
        col_b
      ),
      names=c(name_network_a, name_network_b),
      ylab="gene percent",
      las=1
    )
  } else {
    genes_interest_innetwork = "No genes of interest provided."
  }
  
  #Gene behavior
  print("Gene Behaviour plot")
  plot(
    main = "% emitted connections/gene across networks",
    inout_ab$ratio.a,
    inout_ab$ratio.b,
    pch = 1,
    col = alpha(col_a,0.4),
    cex = 0.75,
    xlab = name_network_a,
    ylab = name_network_b
  )
  
  print("Gene Behaviour Box plot")
  boxplot(
    x = list(
      a = stats_a$In_Out_per_Gene$ratio[
        stats_a$In_Out_per_Gene$ratio > 0
      ],
      b = stats_b$In_Out_per_Gene$ratio[
        stats_b$In_Out_per_Gene$ratio > 0
      ]
    ),
    col = c(col_a,col_b),
    main = "gene behavior in networks",
    sub = paste0(
      "Wilcox p.value ",
      wilcox.test(
        x = stats_a$In_Out_per_Gene$ratio[
          stats_a$In_Out_per_Gene$ratio > 0
        ],
        y = stats_b$In_Out_per_Gene$ratio[
          stats_b$In_Out_per_Gene$ratio > 0
        ]
      )$p.value
    )
  )
  
  # Centrality changes
  print("Changes in centrality for ALL genes")
  plot(
    main = "Changes in centrality across networks",
    centrality_per_network$centrality_a,
    centrality_per_network$centrality_b,
    xlab="network a",
    ylab = "network b",
    pch = 19,
    cex = 0.75,
    col = alpha(col_a,0.25)
  )
  # text(names of genes in top 10-15 of either one or the other network)
  
  # foldchange centrality of all tfs
  print("Changes in centrality for TF genes")
  plot(
    main = "Changes in TF centrality across networks",
    x = tfcntr_ab$a,
    y = tfcntr_ab$b,
    xlab="network a",
    ylab = "network b",
    pch = 19,
    cex = 0.75,
    col = alpha(col_a,0.25)
  )
  print("Foldchange in centrality of TFs")
  plot(
    main = "Foldchange in centrality of TFs",
    sort(tfcntr_ab$b/tfcntr_ab$a) # highlight here where are the topcentrals of each network using dynamic labels or something similar
  )
  
  print("graph plots of TFs that change")
  par(mfrow = c(1,2))
  plot(
    main = paste0("TFs ",name_network_a),
    graph_tf_a,
    vertex.size = 5*V(graph_tf_a)$size,
    vertex.color = col_a,
    vertex.label.cex = 0.6,
    edge.color = "gray",
    edge.arrow.size = 0.2,
    layout = layout_nicely(graph_tf_a)
  )
  plot(
    main = paste0("TFs ",name_network_b),
    graph_tf_b,
    vertex.color = col_b,
    vertex.size = 5*V(graph_tf_b)$size,
    vertex.label.cex = 0.6,
    edge.color = "gray",
    edge.arrow.size = 0.2,
    layout = layout_nicely(graph_tf_b)
  )
  par(mfrow=c(1,1))
  
  # foldchange centrality of all transdevs
  print("Changes in trans-dev centrality across networks")
  plot(
    main = "Changes in trans-dev centrality across networks",
    x = td_ab$a,
    y = td_ab$b,
    xlab = name_network_a,
    ylab = name_network_b,
    pch = 19,
    cex = 0.75,
    col = alpha(col_a,0.25)
  )
  print("Foldchange % emitted connections per gene")
  plot(
    main = "Foldchange % emitted connections per gene",
    sort(td_ab$b/td_ab$a), # highlight here where are the topcentrals of each network using dynamic labels or something similar
    ylab = "Foldchange"
  )
  
  print("graph plots of trans-dev that change")
  par(mfrow = c(1,2))
  plot(
    main = paste0("Trans-Devs ",name_network_a),
    graph_td_a,
    vertex.size = 5*V(graph_td_a)$size,
    vertex.color = col_a,
    vertex.label.cex = 0.6,
    edge.color = "gray",
    edge.arrow.size = 0.2,
    layout = layout_nicely(graph_td_a)
  )
  plot(
    main = paste0("Trans-Devs ",name_network_b),
    graph_td_b,
    vertex.color = col_b,
    vertex.size = 5*V(graph_td_b)$size,
    vertex.label.cex = 0.6,
    edge.color = "gray",
    edge.arrow.size = 0.2,
    layout = layout_nicely(graph_td_b)
  )
  par(mfrow=c(1,1))
  
  
  #Influence
  if ( is.null(influence) == FALSE ){
    print("Influence")
    
    influence_results <- influ_table(influence,tfs)
    
    print("Influence graph")
    influence_graph <- induced.subgraph(
      graph_b,
      vids = influence_results$factor,
      impl = "create_from_scratch"
    )
    
    V(influence_graph)$size <- 
      influence_results$factor_fc[
        match(V(influence_graph)$name, influence_results$factor)
      ]
    V(influence_graph)$color[
      V(influence_graph)$color == ""
    ] <- "gray"
    
    E(influence_graph)$weight <- E(influence_graph)$width
    E(influence_graph)$weight <-
      (E(influence_graph)$weight - min(E(influence_graph)$weight)) / 
      (max(E(influence_graph)$weight) - min(E(influence_graph)$weight))
    
    E(influence_graph)$width <- 2*E(influence_graph)$weight
    E(influence_graph)$color <- alpha("#2c2c2c",0.5*E(influence_graph)$weight)
    
    print("Plot influence graph")
    plot(
      main = "Influence network",
      influence_graph,
      vertex.label = NA,
      edge.arrow.size = 0.2,
      layout = layout_nicely(influence_graph)
    )
  } else {
    influence_results = "no Influence provided."
    influence_graph = "no Influence provided."
  }
  
  res = list(
    comparison_table = compare_df,
    target_genes = list_targets,
    connections_per_tgt_gene = list(
      connects_per_tgt_gene_a = table(nw_a$tg),
      connects_per_tgt_gene_b = table(nw_b$tg)
    ),
    emit_recv_connections = inout_ab,
    tf_centrality_across_networks = tfcntr_ab,
    transdev_centrality_across_networks = td_ab,
    graph_tfs_a = graph_tf_a,
    graph_tfs_b = graph_tf_b,
    graph_transdevs_a = graph_td_a,
    graph_transdevs_b = graph_td_b,
    genes_interest_innetwork,
    GOs_targets = GOs_targets,
    centrality_per_network = centrality_per_network,
    influence_results = influence_results,
    influence_graph = influence_graph
  )
  
  return(res)
}

# 3. funcion para plotear la centrality de los TFs over time (fuzz/broadlines depicting quantile or median or whatever)


#' 1.10 concordancia motivos atac homer y motivos de los tfs de la network (los + centrales, o los q tienen mas targets, etc.)

