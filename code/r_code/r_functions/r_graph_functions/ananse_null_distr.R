ananse_null_distribution_degree <- 
  function(g, tfs, n_genes, n_tfs, N = 100000, seed = 1234, verbose_subsampling = FALSE, message_every_N = 100, plots = FALSE, append_list_graphs = FALSE){
    require(igraph)
    
    is_tf = V(g)$name %in% tfs
    pool1 <- V(g)$name[is_tf]
    pool2 <- V(g)$name[!(is_tf)]
    
    n_tgs <- n_genes - n_tfs
    
    lg <- list()
    j = 0
    
    message("Starting subsampling, ",N," iterations, ", n_genes," genes, of which forcing ", n_tfs, " to be TFs")
    set.seed(seed)
    for(i in 1:N){
      vids <- c(sample(pool1,n_tfs),sample(pool2,n_tgs))
      lg[[i]] <- induced_subgraph(
        g,
        vids = V(g)[V(g)$name %in% vids],
        impl = "auto"
      )
      
      if(verbose_subsampling){
        if(i - j == message_every_N){
          message(i)
          j = j + message_every_N
        }
      }
    
    }
    
    message("calculating the relative frequency of disconnected genes...")
    null_dist <-
      sapply(
        lg,
        function(x){
          length(which(igraph::degree(x)==0)) / vcount(x) # how many genes are disconnected in relation to the number of genes in each graph
        }
      )
    null_dist <- null_dist[complete.cases(null_dist)]
    
    p_exp = length(null_dist[null_dist==0])/length(null_dist)
    message("Expected probability of finding ",n_genes," genes (", n_tfs, " of which are TFs) all connected:\n","p(0) = ", p_exp)
    
    if(plots == TRUE){
      plot(density(null_dist), main = "null distribution", ylab = paste0("Density ; p(0) = ", p_exp))
    }
    
    res = list(
      null_distribution = null_dist,
      expected_probability = p_exp
    )
    
    if(append_list_graphs) res = append(res,list(lg))
    
    return(res)
  }
