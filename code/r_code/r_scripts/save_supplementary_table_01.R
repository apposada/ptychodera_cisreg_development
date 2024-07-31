# to run this, one needs to load before the following: pfla_deseq2 rda, pfla_graphs rda (for attributes list), ALGs rda, markdown 02l rda for orthologues

pfla_supp01 <-
  data.frame(
    transcript_id = rownames(pfla_rna_dev),
    pfla_rna_dev
  )

colnames(pfla_supp01) <- gsub("X","",colnames(pfla_supp01))

pfla_supp01 <-
  merge(pfla_supp01, xloc_tcons, by.x = 1, by.y = 1, all.x = TRUE)

pfla_supp01 <- 
  merge(
    pfla_supp01,
    pfla_attributes_list[[2]][,c(1,2)],
    by.x = 1, by.y = 1, all.x = TRUE
  )

pfla_supp01 <- 
  merge(
    pfla_supp01,
    pfla_attributes_list[[4]][,c(1,2)],
    by.x = 1, by.y = 1, all.x = TRUE
  )

pfla_supp01 <- 
  merge(
    pfla_supp01,
    pfla_age[,c(1,3)],
    by.x = 1, by.y = 1, all.x = TRUE
  )

pfla_supp01 <- 
  merge(
    pfla_supp01,
    pfla_attributes_list[[6]][,c(1,2)],
    by.x = 1, by.y = 1, all.x = TRUE
  )

pfla_supp01 <-
  merge(
    pfla_supp01,
    deut_or,
    by.x = 1, by.y = 2,
    all.x = TRUE
  )

pfla_supp01[is.na(pfla_supp01)] <- "-"

colnames(pfla_supp01) <- c(
  "transcript_id","00_Uf", "01_16", "02_EB", "03_LB", "04_EG","05_MG", "06_MG", "07_LG", "08_To", "09_He", "10_Me", "11_Kr", "12_Sp", "13_Ag", "14_TL", "15_Ad", "temporal_cluster", "gene_id",  "TF_class", "transdev_housekeeping",  "predicted_OMA_gene_age", "predicted_eggNOG_gene_name", "ALG_LinEtAl2024", "bflo_deut_orthologue", "spur_id_orthologue", "deuterostome_orthologue_triplet_id"
)

pfla_supp01 <- pfla_supp01[,c(19,1:18,20:22,24:27,23)]

write.table(
  pfla_supp01,
  file = "outputs/20240520_supp_material_01.tsv",
  sep = "\t",
  dec = ".",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
