comparABle <- function(
  a_name,
  b_name,
  a,
  b,
  o,
  f,
  ma,
  mb,
  ga,
  gb,
  cog_a,
  cog_b,
  a_samples,
  b_samples,
  cooc_n = 1000,
  cooc_h = c(0.75,0.95),
  cooc_clustering_algorithm = "hclust",
  cooc_clustering_method = "average",
  cooc_cor_method = "pearson",
  cooc_p = 0.1,
  cooc_vargenes = rownames(merge_ab$ab_o),
  highlyvariable = TRUE,
  common_evo_nodes,
  a_universe,
  a_id2go,
  sep = ",\ ",
  ...
) {
  
  # TidyUp
  print("Tidy up data")
  samples_a = a_samples
  if (any(apply(a,2,is.character)) == TRUE) {
    a = a[,!(sapply(a, is.character))]
  }
  a = qnorm(a)
  a = tidyup(a, highlyvariable = highlyvariable) # remove genes with 0 tpms
  a = rep2means(samples_a,a)
  
  if (any(apply(b,2,is.character)) == TRUE) {
    b = b[,!(sapply(b, is.character))]
  }
  samples_b = b_samples
  b = qnorm(b)
  b = tidyup(b, highlyvariable = highlyvariable)
  b = rep2means(samples_b,b) # remove genes with 0 tpms
  
  o = pair2id(o)
  
  colnames(ma) <- c("id","module")
  colnames(mb) <- c("id","module")
  
  # MERGE
  print ("Merge data")
  merge_ab <- mergedata(a,b,o)
  
  # CORRELATIONS
  print("Correlations")
  cors <- rawcorsp(merge_ab$a_o,merge_ab$b_o) # FIX JSD
  
  print("PCA")
  pi <- prcomp(t(merge_ab$ab_o))
  
  # CO-OCCURENCE MATRIX
  print("Co-Occurrence")
  set.seed(4343)
  cooc <- treeFromEnsembleClustering(
    x=merge_ab$ab_o, p=cooc_p, h=cooc_h,  n = cooc_n, vargenes = rownames(merge_ab$ab_o), bootstrap=FALSE,
    clustering_algorithm=cooc_clustering_algorithm, clustering_method=cooc_clustering_method, 
    cor_method=cooc_cor_method
  )
  
  # Co-occurrence heatmap
  cooc_hm <- Heatmap(
    cooc$cooccurrence,
    col = colorRamp2(
      c(seq(min(cooc$cooccurrence),
            max(cooc$cooccurrence),
            length=9
      )
      ),
      colors=c(
        c("white",'#ffffe5','#fff7bc','#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#990000')
      )
    ), 
    name="co-occurence"
  )
  
  # COMMON GENES IN CORRELATIONS
  print("Common genes in Correlations")
  ab_common_genes_cor <- get_high_cor_genes(
    mat = cors$js,
    a_o = merge_ab$a_o,
    b_o = merge_ab$b_o,
    o = o
  )
  
  print("Common genes in Correlations (GO)")
  ab_common_genes_cor_GOs <- 
    getGOs(
      genelist = 
        lapply(
          ab_common_genes_cor$hicor_topgenes, 
          function(sub_list) {
            setNames(sub_list$top_genes$a, names(sub_list))
            }),
      gene_universe = a_universe,
      gene2GO = a_id2go
    )
  
  print("Common genes in Correlations (age)")
  ab_common_genes_cor_age <- 
    lapply(
      lapply(
        ab_common_genes_cor$hicor_topgenes,
        function(sub_list) {
          data.frame(
            id = a_universe,
            module = 
              ifelse(
                a_universe %in% sub_list$top_genes$a, "common","not_common"
              )
          )
        }
      ),
      function(x){
        gene_age_enrichment(
          x_modules = x,
          x_age = ga[ga$age %in% common_evo_nodes,]
        )
      }
    )
  
  # COMPARE MODULES
  print("Pairwise Orthology Overlap Strategy across modules -- hypergeometric and binonmial tests")
  modulecomp_ab <- comparemodules(ma,mb,f)
  
  # Genes in key fams across modules
  print("Getting info on genes from shared families across modules")
  ab_common_genes_details <- genes_in_key_fams(
    stats = modulecomp_ab$stats,
    top_comparisons = 20,
    f = f,
    ma = ma,
    mb = mb,
    age_a = ga[ga$age %in% common_evo_nodes,], # custom vector should be a variable of common evol.nodes
    cog_a = cog_a,
    gene2go_a = a_id2go,
    a_universe = a_universe,
    universe = a_universe,
    sep = sep,
    module_a = FALSE,
    module_b = FALSE,
    common = TRUE,
    exclusive = FALSE , 
    same_species = FALSE
  )
  
  
  # PLOTS
  print("Plots")
  # STORE ALL PLOTS IN FUNCTIONS AS DONE BY @cartwheel ON https://stackoverflow.com/questions/29583849/save-a-plot-in-an-object
  # PCA
  ab_pca <- function(){plot_pca_ab(pca = pi,ab_o = merge_ab)}
  
  # Correlation, Coocurrence
  ab_spearman <- Heatmap(
    column_title = paste0("correlation ",a_name," vs ",b_name),
    cors$sp,
    name = "spearman",
    cluster_rows = F,
    cluster_columns = F,
    col = rev(sequential_hcl(10,"YlOrRd"))
  )
  
  ab_pearson <- Heatmap(
    cors$pe,
    name = "correlation",
    cluster_rows = F,
    cluster_columns = F,
    col = rev(sequential_hcl(10,"YlOrRd"))
  )
  
  ab_jsd <- Heatmap(
    cors$js,
    name = "correlation",
    cluster_rows = F,
    cluster_columns = F,
    col = rev(sequential_hcl(10,"YlOrRd"))
  )
  
  # Comparison of modules, binomial test
  # turn this into complexheatmap
  modulecomp_ab_logbinom_hm <- Heatmap(
    modulecomp_ab$logbinom,
    col = rev(sequential_hcl(10,"Blues 3")),
    cluster_columns = F,
    cluster_rows = F
  )
  
  # Comparison of modules, hypergeometric test
  modulecomp_ab_loghypg_hm <- Heatmap(
    modulecomp_ab$loghypg,
    col = rev(sequential_hcl(10,"YlOrBr")),
    cluster_columns = F,
    cluster_rows = F
  )
  
  print("Generating results")
  
  res <- list(
    input = list(
      a = a,
      b = b,
      o = o,
      f = f,
      ma = ma,
      mb = mb,
      ga = ga,
      cog_a = cog_a
    ),
    merged_data = merge_ab,
    pairwise_correlations = cors,
    pca_analysis = pi,
    coocurrence_analysis = cooc,
    high_corr_genes = list(
      pairwise_data = ab_common_genes_cor,
      GOs = ab_common_genes_cor_GOs,
      age = ab_common_genes_cor_age
    ),
    orthology_overlap_modules = list(
      pairwise_module_comparison = modulecomp_ab,
      genes_in_common_fams = ab_common_genes_details
    ),
    plots = list( # these functions will NOT work outside the function call I think... possible fix: make them store the data they plot
      pca = ab_pca(),
      spearman_cor = ab_spearman,
      pearson_cor = ab_pearson,
      jensen_shannon = ab_jsd,
      coocurrence_Heatmap = cooc_hm,
      orthology_overlap_binomial_hm = modulecomp_ab_logbinom_hm,
      orthology_overlap_hypgeom_hm = modulecomp_ab_loghypg_hm
    )
  )
  
  return(res)
}