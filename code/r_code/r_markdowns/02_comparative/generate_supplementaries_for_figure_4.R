## SUPP Y.A
h1 <- 
  Heatmap(
    PFLA_SPUR_COMPARISON$pairwise_correlations$pe,
    name = "Pearson",
    cluster_rows = F, cluster_columns = F,
    show_row_names = TRUE,
    col = sequential_hcl(10,"BluYl", rev = TRUE),
    left_annotation = devstages_ha_rows(),
    top_annotation = quick_ha(colnames(PFLA_SPUR_COMPARISON$pairwise_correlations$pe),"Purple-Orange", rev = TRUE)
  )

h2 <- 
  Heatmap(
    PFLA_SPUR_COMPARISON$pairwise_correlations$sp,
    name = "Spearman",
    cluster_rows = F, cluster_columns = F,
    show_row_names = TRUE,
    col = sequential_hcl(10,"YlOrRd", rev = TRUE),
    left_annotation = devstages_ha_rows(),
    top_annotation = quick_ha(colnames(PFLA_SPUR_COMPARISON$pairwise_correlations$sp),"Purple-Orange", rev = TRUE)
  )

h3 <- 
  Heatmap(
    name = "co-occurrence",
    PFLA_SPUR_COMPARISON$coocurrence_analysis$cooccurrence[1:16,17:33],
    cluster_rows = F,
    cluster_columns = F,
    col = sequential_hcl(10,"Heat", rev = TRUE),
    left_annotation = devstages_ha_rows(),
    top_annotation = quick_ha(colnames(PFLA_SPUR_COMPARISON$coocurrence_analysis$cooccurrence[1:16,17:33]),"Purple-Orange", rev = TRUE)
  )

h_list <- h1+h2+h3

pdf(
  file = "graphics/Supp_Y_A.pdf",
  width = 9,
  height = 4
)
draw(h_list, auto_adjust = FALSE)
dev.off()

## SUPP Y.B
go_plots <- 
  plot_grid(
    pfla_spur_hicor_genes_GOs$GOplot[[1]],
    pfla_spur_hicor_genes_GOs$GOplot[[2]],
    PFLA_SPUR_COMPARISON$high_corr_genes$GOs$GOplot[[1]],
    PFLA_SPUR_COMPARISON$high_corr_genes$GOs$GOplot[[2]],
    PFLA_SPUR_COMPARISON$high_corr_genes$GOs$GOplot[[3]],
    PFLA_SPUR_COMPARISON$high_corr_genes$GOs$GOplot[[4]],
    PFLA_SPUR_COMPARISON$high_corr_genes$GOs$GOplot[[5]],
    ncol = 1
  )

pdf(
  file = "graphics/Supp_Y_B.pdf",
  width = 6,
  height = 30
)
go_plots
dev.off()

## SUPP Y.C
pdf(
  file = "graphics/Supp_Y_C.pdf",
  width = 4.5,
  height = 4.5
)
PFLA_SPUR_COMPARISON$plots$orthology_overlap_hypgeom_hm
dev.off()

## SUPP Z.A
h1 <- 
  Heatmap(
    PFLA_BLAN_COMPARISON$pairwise_correlations$pe,
    name = "Pearson",
    cluster_rows = F, cluster_columns = F,
    show_row_names = TRUE,
    col = sequential_hcl(10,"BluYl", rev = TRUE),
    left_annotation = devstages_ha_rows(),
    top_annotation = quick_ha(colnames(PFLA_BLAN_COMPARISON$pairwise_correlations$pe),"Purple-Yellow", rev = TRUE)
  )

h2 <- 
  Heatmap(
    PFLA_BLAN_COMPARISON$pairwise_correlations$sp,
    name = "Spearman",
    cluster_rows = F, cluster_columns = F,
    show_row_names = TRUE,
    col = sequential_hcl(10,"YlOrRd", rev = TRUE),
    left_annotation = devstages_ha_rows(),
    top_annotation = quick_ha(colnames(PFLA_BLAN_COMPARISON$pairwise_correlations$sp),"Purple-Yellow", rev = TRUE)
  )

h3 <- 
  Heatmap(
    name = "co-occurrence",
    PFLA_BLAN_COMPARISON$coocurrence_analysis$cooccurrence[1:16,17:28],
    cluster_rows = F,
    cluster_columns = F,
    col = sequential_hcl(10,"Heat", rev = TRUE),
    left_annotation = devstages_ha_rows(),
    top_annotation = quick_ha(colnames(PFLA_BLAN_COMPARISON$coocurrence_analysis$cooccurrence[1:16,17:28]),"Purple-Yellow", rev = TRUE)
  )

h_list <- h1+h2+h3

pdf(
  file = "graphics/Supp_Z_A.pdf",
  width = 9,
  height = 4
)
draw(h_list, auto_adjust = FALSE)
dev.off()

## SUPP Z.B
go_plots <- 
  plot_grid(
    PFLA_BLAN_COMPARISON$high_corr_genes$GOs$GOplot[[1]],
    PFLA_BLAN_COMPARISON$high_corr_genes$GOs$GOplot[[2]],
    PFLA_BLAN_COMPARISON$high_corr_genes$GOs$GOplot[[3]],
    PFLA_BLAN_COMPARISON$high_corr_genes$GOs$GOplot[[4]],
    PFLA_BLAN_COMPARISON$high_corr_genes$GOs$GOplot[[5]],
    ncol = 1
  )

pdf(
  file = "graphics/Supp_Z_B.pdf",
  width = 6,
  height = 30
)
go_plots
dev.off()

## SUPP Z.C
pdf(
  file = "graphics/Supp_Z_C.pdf",
  width = 4.5,
  height = 4.5
)
PFLA_BLAN_COMPARISON$plots$orthology_overlap_hypgeom_hm
dev.off()

## SUPP AA.A
pf_aj_js <- 
  Heatmap(
    PFLA_AJAP_COMPARISON$pairwise_correlations$js,
    cluster_rows = F,
    cluster_columns = F, 
    show_row_names = TRUE, 
    name = "JSD", 
    col = brewer.pal(10,"RdBu"),
    left_annotation = devstages_ha_rows(),
    top_annotation = quick_ha(colnames(PFLA_AJAP_COMPARISON$pairwise_correlations$js),"Heat", rev = TRUE)
  )
pdf(
  file = "graphics/Supp_AA_A.pdf",
  height = 4.5,
  width = 4.5
)
draw(pf_aj_js)
dev.off()

## SUPP AA.B
scatters <-
  plot_grid(
    PFLA_AJAP_COMPARISON$high_corr_genes$pairwise_data$hicor_topgenes[[5]]$plot_topgenes,
    PFLA_AJAP_COMPARISON$high_corr_genes$pairwise_data$hicor_topgenes[[4]]$plot_topgenes
  )
pdf(
  file = "graphics/Supp_AA_B.pdf",
  height = 4.5,
  width = 9
)
  scatters
dev.off()

## SUPP AA.C
go_plots <-
  plot_grid(
    PFLA_AJAP_COMPARISON$high_corr_genes$GOs$GOplot[[5]],
    PFLA_AJAP_COMPARISON$high_corr_genes$GOs$GOplot[[4]],
    ncol = 2
  )
pdf(
  file = "graphics/Supp_AA_C.pdf",
  height = 4.5,
  width = 20
)
go_plots
dev.off()

## SUPP AA.D
pdf(
  file = "graphics/Supp_AA_D.pdf",
  width = 4.5,
  height = 4.5
)
PFLA_AJAP_COMPARISON$plots$orthology_overlap_hypgeom_hm
dev.off()

## SUPP AA.E
pdf(
  file = "graphics/Supp_AA_E.pdf",
  width = 6,
  height = 12
)
PFLA_AJAP_COMPARISON$orthology_overlap_modules$
  genes_in_common_fams$commonfams$
  cog_a_comon$heatmap
dev.off()

## SUPP AB.A
pf_bf_js <- 
  Heatmap(
    cors$js,
    cluster_rows = F,
    cluster_columns = F, 
    show_row_names = TRUE, 
    name = "JSD", 
    col = brewer.pal(10,"RdBu"),
    left_annotation = devstages_ha_rows(),
    top_annotation = quick_ha(colnames(cors$js), rev = TRUE)
  )
pdf(
  file = "graphics/Supp_AB_A.pdf",
  height = 4.5,
  width = 4.5
)
draw(pf_bf_js)
dev.off()

## SUPP AB.B
scatters <-
  plot_grid(
    # ab_common_genes_cor$hicor_topgenes[[1]]$plot_topgenes,
    # ab_common_genes_cor$hicor_topgenes[[2]]$plot_topgenes,
    # ab_common_genes_cor$hicor_topgenes[[3]]$plot_topgenes,
    ab_common_genes_cor$hicor_topgenes[[4]]$plot_topgenes,
    ab_common_genes_cor$hicor_topgenes[[5]]$plot_topgenes,
    ncol = 2
  )
pdf(
  file = "graphics/Supp_AB_B.pdf",
  width = 9,
  height = 4.5
)
scatters
dev.off()

## SUPP AB.C
go_plots <-
  plot_grid(
    # ab_common_genes_cor_GOs$GOplot[[1]],
    # ab_common_genes_cor_GOs$GOplot[[2]],
    # ab_common_genes_cor_GOs$GOplot[[3]],
    ab_common_genes_cor_GOs$GOplot[[4]],
    ab_common_genes_cor_GOs$GOplot[[5]]
  )
pdf(
  file = "graphics/Supp_AB_C.pdf",
  width = 12,
  height = 4
)
go_plots
dev.off()

## SUPP ORTHOF.A
h1 <- 
  Heatmap(
    PFLA_SPUR_COMPARISON_OF$pairwise_correlations$js,
    name = "JSD",
    cluster_rows = F, cluster_columns = F,
    show_row_names = TRUE,
    col = brewer.pal(10,"RdBu"),
    left_annotation = devstages_ha_rows(),
    top_annotation = quick_ha(colnames(PFLA_SPUR_COMPARISON_OF$pairwise_correlations$js),"Purple-Orange", rev = TRUE)
  )

h2 <- 
  Heatmap(
    PFLA_SPUR_COMPARISON_OF$pairwise_correlations$pe,
    name = "Pearson",
    cluster_rows = F, cluster_columns = F,
    show_row_names = TRUE,
    col = sequential_hcl(10,"BluYl", rev = TRUE),
    left_annotation = devstages_ha_rows(),
    top_annotation = quick_ha(colnames(PFLA_SPUR_COMPARISON_OF$pairwise_correlations$pe),"Purple-Orange", rev = TRUE)
  )

h3 <- 
  Heatmap(
    PFLA_SPUR_COMPARISON_OF$pairwise_correlations$sp,
    name = "Spearman",
    cluster_rows = F, cluster_columns = F,
    show_row_names = TRUE,
    col = sequential_hcl(10,"YlOrRd", rev = TRUE),
    left_annotation = devstages_ha_rows(),
    top_annotation = quick_ha(colnames(PFLA_SPUR_COMPARISON_OF$pairwise_correlations$sp),"Purple-Orange", rev = TRUE)
  )

h4 <- 
  Heatmap(
    name = "co-occurrence",
    PFLA_SPUR_COMPARISON_OF$coocurrence_analysis$cooccurrence[1:16,17:33],
    cluster_rows = F,
    cluster_columns = F,
    col = sequential_hcl(10,"Heat", rev = TRUE),
    left_annotation = devstages_ha_rows(),
    top_annotation = quick_ha(colnames(PFLA_SPUR_COMPARISON_OF$coocurrence_analysis$cooccurrence[1:16,17:33]),"Purple-Orange", rev = TRUE)
  )

h_list <- h1+h2+h3+h4

pdf(
  file = "graphics/Supp_Orthof_A.pdf",
  width = 12,
  height = 3.5
)
draw(h_list, auto_adjust = FALSE)
dev.off()

## SUPP ORTHOF.B
pdf(
  file = "graphics/Supp_Orthof_B.pdf",
  width = 4.5,
  height = 4.5
)
draw(PFLA_SPUR_COMPARISON_OF$plots$orthology_overlap_hypgeom_hm)
dev.off()

## SUPP ORTHOF.C
pdf(
  file = "graphics/Supp_Orthof_C.pdf",
  width = 6,
  height = 7
)
draw(PFLA_SPUR_COMPARISON_OF$orthology_overlap_modules$
  genes_in_common_fams$commonfams$
  cog_a_comon$heatmap)
dev.off()

## SUPP ORTHOF.D
h1 <- 
  Heatmap(
    PFLA_BLAN_COMPARISON_OF$pairwise_correlations$js,
    name = "JSD",
    cluster_rows = F, cluster_columns = F,
    show_row_names = TRUE,
    col = brewer.pal(10,"RdBu"),
    left_annotation = devstages_ha_rows(),
    top_annotation = quick_ha(colnames(PFLA_BLAN_COMPARISON_OF$pairwise_correlations$js),"Purple-Yellow", rev = TRUE)
  )

h2 <- 
  Heatmap(
    PFLA_BLAN_COMPARISON_OF$pairwise_correlations$pe,
    name = "Pearson",
    cluster_rows = F, cluster_columns = F,
    show_row_names = TRUE,
    col = sequential_hcl(10,"BluYl", rev = TRUE),
    left_annotation = devstages_ha_rows(),
    top_annotation = quick_ha(colnames(PFLA_BLAN_COMPARISON_OF$pairwise_correlations$pe),"Purple-Yellow", rev = TRUE)
  )

h3 <- 
  Heatmap(
    PFLA_BLAN_COMPARISON_OF$pairwise_correlations$sp,
    name = "Spearman",
    cluster_rows = F, cluster_columns = F,
    show_row_names = TRUE,
    col = sequential_hcl(10,"YlOrRd", rev = TRUE),
    left_annotation = devstages_ha_rows(),
    top_annotation = quick_ha(colnames(PFLA_BLAN_COMPARISON_OF$pairwise_correlations$sp),"Purple-Yellow", rev = TRUE)
  )

h4 <- 
  Heatmap(
    name = "co-occurrence",
    PFLA_BLAN_COMPARISON_OF$coocurrence_analysis$cooccurrence[1:16,17:28],
    cluster_rows = F,
    cluster_columns = F,
    col = sequential_hcl(10,"Heat", rev = TRUE),
    left_annotation = devstages_ha_rows(),
    top_annotation = quick_ha(colnames(PFLA_BLAN_COMPARISON$coocurrence_analysis$cooccurrence[1:16,17:28]),"Purple-Yellow", rev = TRUE)
  )

h_list <- h1+h2+h3+h4

pdf(
  file = "graphics/Supp_Orthof_D.pdf",
  width = 12,
  height = 3.5
)
draw(h_list, auto_adjust = FALSE)
dev.off()

## SUPP ORTHOF.E
pdf(
  file = "graphics/Supp_Orthof_E.pdf",
  width = 4.5,
  height = 4.5
)
draw(PFLA_BLAN_COMPARISON_OF$plots$orthology_overlap_hypgeom_hm)
dev.off()

## SUPP ORTHOF.F
pdf(
  file = "graphics/Supp_Orthof_F.pdf",
  width = 6,
  height = 7
)
draw(PFLA_BLAN_COMPARISON_OF$orthology_overlap_modules$
  genes_in_common_fams$commonfams$
  cog_a_comon$heatmap)
dev.off()

## SUPP ORTHOF.G
a <- pfla_rna_counts
b <- spur_vsd
a_samples = levels(condition_x)
b_samples = unique(sub("_.$", "", colnames(b)))
o = unique(rbh_pfla_spur)
a = qnorm(a)
a = tidyup(a, highlyvariable = FALSE) # remove genes with 0 tpms
a = rep2means(a_samples,a)
b = qnorm(b)
b = tidyup(b, highlyvariable = FALSE)
b = rep2means(b_samples,b) # remove genes with 0 tpms
o = pair2id(o)
pf_sp_merge_ab <- mergedata(a,b,o)
pf_sp_rbh_cors <- rawcorsp(pf_sp_merge_ab$a_o,pf_sp_merge_ab$b_o) # FIX JSD

h1 <- 
  Heatmap(
    pf_sp_rbh_cors$js,
    name = "JSD",
    cluster_rows = F, cluster_columns = F,
    show_row_names = TRUE,
    col = brewer.pal(10,"RdBu"),
    left_annotation = devstages_ha_rows(),
    top_annotation = quick_ha(colnames(pf_sp_rbh_cors$js),"Purple-Orange", rev = TRUE)
  )

h2 <- 
  Heatmap(
    pf_sp_rbh_cors$pe,
    name = "Pearson",
    cluster_rows = F, cluster_columns = F,
    show_row_names = TRUE,
    col = sequential_hcl(10,"BluYl", rev = TRUE),
    left_annotation = devstages_ha_rows(),
    top_annotation = quick_ha(colnames(pf_sp_rbh_cors$pe),"Purple-Orange", rev = TRUE)
  )

h3 <- 
  Heatmap(
    pf_sp_rbh_cors$sp,
    name = "Spearman",
    cluster_rows = F, cluster_columns = F,
    show_row_names = TRUE,
    col = sequential_hcl(10,"YlOrRd", rev = TRUE),
    left_annotation = devstages_ha_rows(),
    top_annotation = quick_ha(colnames(pf_sp_rbh_cors$sp),"Purple-Orange", rev = TRUE)
  )

h_list <- h1+h2+h3

pdf(
  file = "graphics/Supp_Orthof_G.pdf",
  width = 9,
  height = 3.5
)
draw(h_list,auto_adjust = FALSE)
dev.off()

## SUPP ORTHOF.H
a <- pfla_rna_counts
b <- blan_counts
a_samples = levels(condition_x)
b_samples = unique(sub("_.$", "", colnames(b)))
o = unique(rbh_pfla_blan)
a = qnorm(a)
a = tidyup(a, highlyvariable = FALSE) # remove genes with 0 tpms
a = rep2means(a_samples,a)
b = qnorm(b)
b = tidyup(b, highlyvariable = FALSE)
b = rep2means(b_samples,b) # remove genes with 0 tpms
o = pair2id(o)
pf_bl_merge_ab <- mergedata(a,b,o)
pf_bl_rbh_cors <- rawcorsp(pf_bl_merge_ab$a_o,pf_bl_merge_ab$b_o) # FIX JSD

h1 <- 
  Heatmap(
    pf_bl_rbh_cors$js,
    name = "JSD",
    cluster_rows = F, cluster_columns = F,
    show_row_names = TRUE,
    col = brewer.pal(10,"RdBu"),
    left_annotation = devstages_ha_rows(),
    top_annotation = quick_ha(colnames(pf_bl_rbh_cors$js),"Purple-Yellow", rev = TRUE)
  )

h2 <- 
  Heatmap(
    pf_bl_rbh_cors$pe,
    name = "Pearson",
    cluster_rows = F, cluster_columns = F,
    show_row_names = TRUE,
    col = sequential_hcl(10,"BluYl", rev = TRUE),
    left_annotation = devstages_ha_rows(),
    top_annotation = quick_ha(colnames(pf_bl_rbh_cors$pe),"Purple-Yellow", rev = TRUE)
  )

h3 <- 
  Heatmap(
    pf_bl_rbh_cors$sp,
    name = "Spearman",
    cluster_rows = F, cluster_columns = F,
    show_row_names = TRUE,
    col = sequential_hcl(10,"YlOrRd", rev = TRUE),
    left_annotation = devstages_ha_rows(),
    top_annotation = quick_ha(colnames(pf_bl_rbh_cors$sp),"Purple-Yellow", rev = TRUE)
  )

h_list <- h1+h2+h3

pdf(
  file = "graphics/Supp_Orthof_H.pdf",
  width = 9,
  height = 3.5
)
draw(h_list,auto_adjust = FALSE)
dev.off()

## SUPP 9BARPLOTS
pdf(
  file = "graphics/Supp_9BARPLOTS.pdf",
  width = 8,
  height = 8
)
ggarrange(
  pfla_tfs_ngenes_plot,
  spur_tfs_ngenes_plot,
  blan_tfs_ngenes_plot,
  pfla_tfs_expgenes_plot,
  spur_tfs_expgenes_plot,
  blan_tfs_expgenes_plot,
  pfla_tf_EXPNGEN_plot,
  spur_tf_EXPNGEN_plot,
  blan_tf_EXPNGEN_plot,
  ncol=3,
  nrow=3,
  common.legend=T
)
dev.off()

## SUPP FOXA-A
h1 <- 
  Heatmap(
    pf_sp_rbh_cors_tfs$js,
    name = "JSD",
    cluster_rows = F, cluster_columns = F,
    show_row_names = TRUE,
    col = brewer.pal(10,"RdBu"),
    left_annotation = devstages_ha_rows(),
    top_annotation = quick_ha(colnames(pf_sp_rbh_cors_tfs$js),"Purple-Orange", rev = TRUE)
  )
h2 <- 
  Heatmap(
    pf_bl_rbh_cors_tfs$js,
    name = "JSD",
    cluster_rows = F, cluster_columns = F,
    show_row_names = TRUE,
    col = brewer.pal(10,"RdBu"),
    left_annotation = devstages_ha_rows(),
    top_annotation = quick_ha(colnames(pf_bl_rbh_cors_tfs$js),"Purple-Yellow", rev = TRUE)
  )
h_list <- h1+h2

pdf(
  file = "graphics/Supp_FoxA_A.pdf",
  width = 6,
  height = 3.5
)
draw(h_list, auto_adjust = FALSE)
dev.off()

## SUPP FOXA-B
# plots of correlation changes in centrality tfs and transdev
pdf(
  file = "graphics/Supp_FoxA_B.pdf",
  width = 4,
  height = 4
)
plot(
  main = "Changes in TF centrality\nacross networks",
  x = relativise(EB_vs_LG$tf_centrality_across_networks$a),
  y = relativise(EB_vs_LG$tf_centrality_across_networks$b),
  xlab = "Early Blastula",
  ylab = "Late Gastrula",
  pch = 19,
  cex = 0.75,
  col = alpha("#E58745",0.25)
)
points(
  x = relativise(EB_vs_LG$tf_centrality_across_networks$a)[EB_vs_LG$tf_centrality_across_networks$id == "TCONS_00009611"],
  y = relativise(EB_vs_LG$tf_centrality_across_networks$b)[EB_vs_LG$tf_centrality_across_networks$id == "TCONS_00009611"],
  pch = 21,
  col = "black",
  bg = "purple"
)
text(
  x = relativise(EB_vs_LG$tf_centrality_across_networks$a)[EB_vs_LG$tf_centrality_across_networks$id == "TCONS_00009611"],
  y = relativise(EB_vs_LG$tf_centrality_across_networks$b)[EB_vs_LG$tf_centrality_across_networks$id == "TCONS_00009611"],
  "FoxA",
  pos = 3
)
dev.off()

## SUPP FOXA-C
pdf(
  file = "graphics/Supp_FoxA_C.pdf",
  width = 4,
  height = 4
)
plot(sort(relativise(pfla_LG_stats$Centrality_per_TFclass$Forkhead)), main = "centrality of\nForkhead TFs",type = "p", ylab = "relative centrality", xlab = "")
points(21,1,pch = 21, col = "black", bg = "purple")+text(21,1,"FoxA",pos = 2)
dev.off()

## SUPP FOXA-D
ananse_EB_to_LG <- 
  read.table("outputs/ananse/influence_LG_EB.tsv",header=T) #change this path

ananse_EB_to_LG$factor <-
  sub("TCONS","TCONS_",ananse_EB_to_LG$factor)
influ_tbl <- ananse_EB_to_LG 
influ_tbl <- merge(influ_tbl,pfla_tfs_graph_analysis,by.x=1,by.y=1,all.x=T)
influ_tbl$TFclass[is.na(influ_tbl$TFclass)] <- " "
influ_tbl$col[is.na(influ_tbl$col)] <- "gray"
influ_tbl <- influ_tbl [rev(order(influ_tbl$sumScaled)),]
rownames(influ_tbl) <- NULL

# 'translate' the gene ids to known gene names using a small function
influ_tbl$genename <- 
  translate_ids(x = influ_tbl$factor,dict = pfla_genenames)

pdf(
  file = "graphics/Supp_FoxA_D.pdf",
  width = 4,
  height = 4
)
plot(
  influ_tbl$factor_fc,
  influ_tbl$sumScaled,
  pch=19,
  col=alpha("black",0.15),
  xlab="log2fold change of TF",
  ylab="ANANSE influence score",
  main="Main factors",
  bty="n",
  xlim=c(0,max(influ_tbl$factor_fc)+1)
)
points(
  influ_tbl$factor_fc[influ_tbl$factor == foxa],
  influ_tbl$sumScaled[influ_tbl$factor == foxa],
  pch = 21,
  col = "black",
  bg = "purple"
)
text(
  influ_tbl$factor_fc[influ_tbl$factor == foxa],
  influ_tbl$sumScaled[influ_tbl$factor == foxa],
  "FoxA",
  pos = 2
)
dev.off()

pf_annot <- devstages_ha_rows()
sp_annot <- quick_ha(colnames(pf_sp_rbh_cors_tfs$js),"Purple-Orange", rev = TRUE)
bl_annot <- quick_ha(colnames(pf_bl_rbh_cors_tfs$js),"Purple-Yellow", rev = TRUE)
aj_annot <- quick_ha(colnames(PFLA_AJAP_COMPARISON$pairwise_correlations$js),"Heat", rev = TRUE)
bf_annot <-  quick_ha(unique(sub("_.$", "", colnames(bflo_counts))), rev = TRUE)
draw(pf_annot)
draw(sp_annot)
draw(bl_annot)
draw(aj_annot)
draw(bf_annot)
