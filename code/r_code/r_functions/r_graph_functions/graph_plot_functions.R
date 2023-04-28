#' Plotting functions
density_network <- function(nw, q = 0.85) {
  den <- density(nw$prob)
  plot(den, col = "#613f86", main = "Distribution of network scores")
  polygon(den, border = "#613f86", col= "#ccabf4")
  polygon(
    c(den$x[den$x >= quantile(nw$prob,q)], quantile(nw$prob,q)),
    c(den$y[den$x >= quantile(nw$prob,q)], 0),
    col = "#881eae",
    border = "#613f86"
  )
  abline(v = quantile(nw$prob, q), lwd = 1.5, col = "red")
}

plot_numgenes_network <- function(s){
  barplot(
    s$numgenes_network, 
    ylim = c(0,s$numgenes_network + s$numgenes_network * 0.4),
    main = "Network\nsize"
  )
}

plot_typegenes_network <- function(s){
  barplot(
    log(c(
      s$num_active_TFs,
      s$num_total_TFs,
      s$num_Effector_Genes,
      s$num_TGs
      )+1,10),
    col = colorRampPalette(c("#5e4589","mediumpurple2"))(4),
    border=NA,
    names.arg = c(
      "Active TFs",
      "Total TFs",
      "Effector Genes",
      "Target genes"
    ),
    las = 2,
    cex.lab = 0.75,
    main = "Types of Genes in Network"
  )
}

plot_connections <- function(s){
  barplot(
    log(c(
      s$mean_connection_per_TF,
      s$mean_connection_per_TG
    )),
    col = c("red","blue"),
    main = "Mean connections per gene",
    names.arg = c("stemming\nconnections", "received\nconnections")
  )
}

plot_connections_per_gene <- function(s){
  boxplot(
    list(
      "stemming\nconnections" = log(s$connection_per_TF),
      "received\nconnections" = log(s$connection_per_TG)
    ),
    col = c("red","blue"),
    main = "Connections per gene"
  )
}

plot_genes_ratioconnections <- function(s){
  boxplot(
    log(s$In_Out_per_Gene$ratio),
    main = "Genes by out/in\nratio of connections",
    col = "#ccabf4"
  )
}

plot_ratio_connections <- function(s, remove_zeroes = TRUE){
  if (remove_zeroes == TRUE){
    data <- s$In_Out_per_Gene$ratio[
      s$In_Out_per_Gene$ratio > 0
    ]
  } else {
    data <- s$In_Out_per_Gene$ratio
  }
  
  plot(
    seq(1:length(data)),
    data,
    main = "Gene behavior\n(types & amt of connections)",
    col = colorRampPalette(
      c("blue", "red")
    )(length(data)),
    xlab = "",
    ylab = "% emitting connections"
  ) # + abline(h = mean(s$In_Out_per_Gene$ratio))
}

plot_selfreg <- function(s){
  barplot(
    s$num_selfregulated_TFs,
    ylim = c(0,s$num_selfregulated_TFs+s$num_selfregulated_TFs*0.4),
    main = "Number of\nself-regulating TFs"
  )
}

plot_size_CC <- function(s){
  barplot(
    c(
      s$Connected_component_sizes
    ),
    main = "Sizes of ten largest\nconnected components"
  )
}

plot_centrality_TFclass <- function(s,f, sub_s = NULL, main = "Centrality per TF class", ...){
  
  if (!missing(sub_s)) {
    s <- s[names(s) %in% sub_s]
    f <- f[names(f) %in% sub_s]
  }
  f <- f[names(f) %in% names(s)]
  s <- s[!(names(s) == "")]
  
  s <- s[order(names(s))]
  f <- f[order(names(f))]
  
  par(cex.axis = 0.7)
  
  boxplot(
    s,
    col = f[order(names(f))],
    main = main,
    ylab = "Centrality",
    las = 2,
    ...
  )
  
}

plot_centrality_TDHK <- function(s){
  boxplot(
    s$Centrality_per_TransDev_HK[c(2,3)],
    col = c("purple","orange"),
    main = "Centrality of\ntrans-dev genes"
  )
}

plot_top_central_tfs <- function(s){
  barplot(
    rev(s$most_central_TFs[1:15]),
    las=2,
    horiz=T,
    col = "#ccabf4",
    border = NA,
    main = "Top Central TFs"
  )
}

plot_top_central_transdev <- function(s){
  barplot(
    rev(s$most_central_transdevs)[1:15],
    las=2,
    horiz=T,
    col = "#ccabf4",
    border = NA,
    main = " Top Central Trans-Dev genes"
  )
}

plot_inedges_per_category <- function(s){
  boxplot(
    s$edges_per_funcat_family$recv ~
      s$edges_per_funcat_family$funcat,
    main = "Most regulated functional categories",
    col = rainbow(length(unique(s$edges_per_funcat_family$funcat)))
  )
}

plot_outedges_per_category <- function(s){
  boxplot(
    s$edges_per_funcat_family$emit ~
      s$edges_per_funcat_family$funcat,
    main = "Emitting edges by functional category",
    col = rainbow(length(unique(s$edges_per_funcat_family$funcat)))
  )
}

plot_behavior_per_category2 <- function(s){
  a <- s$edges_per_funcat_family[,c(4,5)]
  a <- a[a$funcat != "",]
  a$funcat <- factor(a$funcat, levels = LETTERS)
  a$ratio[a$ratio > min(a$ratio)] <- 1
  a <- a[a$ratio == 1,]
  x <- as.vector(table(a$funcat))
  names(x) <- names(table(a$funcat))
  x <- log(x,10)
  x[is.infinite(x)] <- 0
  barplot(
    x,
    col = rainbow(length(levels(a$funcat)))
  )
}


plot_behavior_per_category <- function(s){
  boxplot(
    s$edges_per_funcat_family$ratio ~
      s$edges_per_funcat_family$funcat,
    main = "Regulatory behavior by\nfunctional category",
    col = rainbow(length(unique(s$edges_per_funcat_family$funcat))) # select based on available COG category data
  )
}

plot_GOs <- function(listplots) {
  for (i in listplots){
    plot(i)
  }
}


default_colpal <- c(
  base = "",
  outline = "",
  filling = "",
  highlight_outline = "",
  highlight_filling = "",
  shade = "",
  outgoing = "red",
  incoming = "blue",
  setA = "purple",
  setB = "orange"
)

