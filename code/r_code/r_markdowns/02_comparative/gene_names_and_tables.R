load("outputs/rda/species_comparison.rda")

pfla_genenames <-
  read.delim2(
    file = "outputs/functional_annotation/eggnog/emapper.annotations",
    skip = 3,
    header = TRUE
  )[,c(1,5,13)]
pfla_genenames <- pfla_genenames[pfla_genenames$predicted_gene_name != "",]


Pfla_Blan_Spur_comparison_genes <-
  rbind(
    cbind(
      pfla_genenames[pfla_genenames$X.query_name %in% PFLA_BLAN_COMPARISON$high_corr_genes$pairwise_data$hicor_topgenes$cor_05_MG__08_18h$top_genes$a,],
      comparison = "Pfla MG/Blan 18h"
    ),
    cbind(
      pfla_genenames[pfla_genenames$X.query_name %in% PFLA_BLAN_COMPARISON$high_corr_genes$pairwise_data$hicor_topgenes$cor_10_Me__14_Premet$top_genes$a,],
      comparison = "Pfla Metchk/Blan Premet"
    ),
    cbind(
      pfla_genenames[pfla_genenames$X.query_name %in% pfla_spur_hicor_genes$hicor_topgenes$cor_04_EG__03_MeBl_24h1$top_genes$a,],
      comparison = "Pfla EG/Spur MeBl.1"
    ),
    cbind(
      pfla_genenames[pfla_genenames$X.query_name %in% pfla_spur_hicor_genes$hicor_topgenes$cor_04_EG__03_MeBl_24h2$top_genes$a,],
      comparison = "Pfla EG/Spur MeBl.2"
    )
  )

colnames(Pfla_Blan_Spur_comparison_genes) <-
  c(
    "Pfla_geneID", "predicted_gene_name",
    "eggNOG description", "comparison"
  )

write.table(
  Pfla_Blan_Spur_comparison_genes,
  file = "outputs/Pfla_Blan_Spur_comparison_genes.tsv",
  sep = "\t",
  dec = ".",
  quote = F,
  row.names = F
)


pf_bl_ma__mb <- PFLA_BLAN_COMPARISON$orthology_overlap_modules$genes_in_common_fams$commonfams$table_a_common

pf_bl_ma__mb_genenames <-
  merge(
    pfla_genenames,
    pf_bl_ma__mb,
    by.x = 1,
    by.y = 1
  )

pf_bl_ma__mb_genenames <-
  pf_bl_ma__mb_genenames[
    pf_bl_ma__mb_genenames$module %in% 
      c(
        "14__18","15__21","17__21","17__22","17__24",
        "19__21","19__22","19__24","20__21","20__24",
        "21__21","22__21","22__24"
      ),
  ]

colnames(pf_bl_ma__mb_genenames) <-
  c(
    "Pfla_geneID", "predicted_gene_name",
    "eggNOG description", "pair of stage-specific clusters"
  )

write.table(
  pf_bl_ma__mb_genenames,
  file = "outputs/Pfla_Blan_commongenes_latedevelopmentClusters.tsv",
  sep = "\t",
  dec = ".",
  quote = F,
  row.names = F
)
write.table(
  rbind(
    data.frame(
      stage = colnames(PFLA_BLAN_COMPARISON$input$a),
      species = "P.flava"
    ),
    data.frame(
      stage = colnames(PFLA_BLAN_COMPARISON$input$b),
      species = "B.lanceolatum"
    ),
    data.frame(
      stage = colnames(PFLA_SPUR_COMPARISON$input$b),
      species = "S.purpuratus"
    )
  ),
  file = "graphics/stages.tsv",
  sep = "\t",
  quote = F,
  row.names = F
)
