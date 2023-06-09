# Basic functions
computeInVsOut <- function(n){
  a <- merge(
    as.data.frame(table(n$tf)),
    as.data.frame(table(n$tg)),
    by.x = 1,
    by.y = 1,
    all = T
  )
  colnames(a) <- c("id", "emit", "recv")
  a[is.na(a)] <- 0
  a$ratio <- a$emit / (a$emit + a$recv)
  a <- a[order(a$ratio),]
  return(a)
}

nw_to_hm <- function(nw) {
  hm <- matrix(
    0,
    nrow = length(unique(nw$tf)),
    ncol = length(unique(nw$tg))
  )
  rownames(hm) <- unique(nw$tf)
  colnames(hm) <- unique(nw$tg)
  for (i in unique(nw$tf)){
    a <- nw[nw$tf == i,]
    for (j in a$tg){
      val = a$prob[a$tg == j]
      hm[
        rownames(hm) == i,
        colnames(hm) == j
      ] <- val
    }
  }
  hm <- hm[order(rownames(hm)),]
  return(hm)
  # test <- nw_to_hm(pfla_EG_nw[
  #   pfla_EG_nw$tf %in% V(pfla_EG_stats$main_Connected_Component_graph)$name &
  #     pfla_EG_nw$tg %in% V(pfla_EG_stats$main_Connected_Component_graph)$name,
  # ])
  # heatmap(1-test)
  
}

getcentral <- function(g, d, category, ...){ # a much, much quicker option would be to calculate centrality for ALL genes, and then define 'operators' to retrieve them and put them in lists
  l <- list()
  for (i in unique(d[[category]]) ){
    v <- d$index[d[[category]] == i]
    l[[i]] <- closeness(
      g,
      vids = v,
      mode = c("all"),
      weights = NULL,
      normalized = F
    )
  }
  return(l)
}


#' Params
#' @genelist must be a object of type list() containing an undetermined number of character vectors. Each component of the list is a char vector of gene names. This function will perform gene ontology enrichment using the char vector as test sample and the totality of your organism's genes as population (Universe). 
#' @gene_universe is a vector of characters containing the genes you want to use as universe. This is normally the totality of genes in your species, or the totality of genes with recollected expression in your dataset.
#' @alg is the algorithm of choice, default is 'Classic', but 'elim' is advised
#' @cols is a char vector of the colors used in the barplots
#' 
require(topGO)
# getGOs(list(tfs,egs,active_tfs,target_genes,V(graph)$name),toy_universe)
getGOs <- function(
  genelist,
  gene_universe,
  gene2GO,
  alg = "Classic",
  stat = "fisher",
  cols=rainbow(length(genelist))
) {
  res <- list(
    GOtable = list(),
    GOplot = list()
  )
  geneID2GO <- gene2GO
  for (j in 1:length(genelist)){
    x <- genelist[[j]]
    numgenes <- length(x)
    allgenes <- data.frame(
      id = gene_universe)
    allgenes$diffreg <- 0
    allgenes$diffreg[allgenes$id %in% x] <- 1
    genelist_x <- allgenes$diffreg
    names(genelist_x) <- allgenes$id
    signif <- function (allScore) {return(allScore == 1)} # function for TopGO which retrieves the significant genes as those of interest
    
    #load data
    GOdata_x <- new(
      "topGOdata",
      ontology = "BP",
      allGenes = genelist_x,
      geneSelectionFun = signif,
      annot = annFUN.gene2GO,
      gene2GO = geneID2GO
    )
    
    #test
    res_Fisher.x <- runTest(GOdata_x, algorithm = alg, statistic = "fisher")
    
    #table
    allres_x <- GenTable(
      GOdata_x,
      classicFisher = res_Fisher.x,
      orderBy = "classicFisher",
      ranksOf = "classicFisher",
      topNodes = 30,
      numChar = 5000
    )
    allres_x <- allres_x[allres_x$Annotated > 5,]
    allres_x$classicFisher[allres_x$classicFisher=="< 1e-30"] <- 1e-30
    
    maxgenesplot <- ifelse(nrow(allres_x) > 10, 10, nrow(allres_x))
    res$GOtable[[j]] <- allres_x
    
    #plot
    res$GOplot[[j]] <- barplot(
      height = as.vector(
        rev(
          -log(
            as.numeric(allres_x$classicFisher[1:maxgenesplot]) * 
              length(x)
          )
        )
      ),
      horiz=T,
      border=F,
      names.arg=rev(allres_x$Term[1:maxgenesplot]),
      sub="-log10 p-value",
      main="GOs (elim,\n> 5 annot genes)",
      las=1,
      col=cols[j],
      space=c(0.5),
      cex.names=.7,
      plot=FALSE
    )
  }
  names(res$GOtable) <- names(genelist)
  return(res)
}

# ANANSE influence analysis
influ_table <- function(x,tftable) {
  influ_tbl <- x 
  influ_tbl <- merge(influ_tbl,tftable,by.x=1,by.y=1,all.x=T)
  influ_tbl$TFclass[is.na(influ_tbl$TFclass)] <- " "
  influ_tbl$col[is.na(influ_tbl$col)] <- "gray"
  influ_tbl <- influ_tbl [rev(order(influ_tbl$sumScaled)),]
  rownames(influ_tbl) <- NULL
  plot(
    influ_tbl$factor_fc,
    influ_tbl$sumScaled,
    pch=19,
    col=influ_tbl$col,
    bg="black",
    xlab="log2fold change of TF",
    ylab="ANANSE influence score",
    main="Main factors",
    bty="n",
    xlim=c(0,max(influ_tbl$factor_fc)+1)
  )
  text(
    influ_tbl$factor_fc[1:20],
    influ_tbl$sumScaled[1:20]+0.01,
    paste(influ_tbl$factor[1:20],influ_tbl$TFclass[1:20], sep = " "),
    cex=0.7
  )
  return(influ_tbl)
}
