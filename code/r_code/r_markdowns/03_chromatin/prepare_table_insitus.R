# 20211029_list_endomesoecto_genes

pfla_ISH <- 
  read.delim(
    "outputs/functional_annotation/germlayers_ISH/20211029_ecto_meso_endo_insitu.tsv",
    header = F, sep = "\t",
    col.names=c("gene","TSO","geneID","Accession","BlastDescription","Endoderm","Mesoderm","Ectoderm","Reference","Provider","Notes"))

xloc_tcons <- 
  read.table(
    "/home/ska/aperpos/Def_Pfla/data/dropbox_collaborators/annotation/v20210202/20210202_20211029_xloc_to_tcons.txt",
    header = F,
    col.names = c("id","geneID")
    )

pfla_ISH$id <-
  translate_ids(x = pfla_ISH$geneID, dict = xloc_tcons[,c(2,1)])


germlayers_colors <- 
  data.frame(
    germlayer = c(
      "Ectoderm",
      "EndoEcto",
      "All",
      "EctoMeso",
      "Mesoderm",
      "EndoMeso",
      "Endoderm",
      "None"
      ), 
    germlayer_col = c(
      '#a8d1f1', # ectoderm, blue
      '#d4f7e2' , # endo + ectoderm
      '#bcbcbc', # all three ,gray
      '#d9d2e9', # ecto + mesoderm
      '#f4cccc',	#mesoderm, red
      '#fce5cd', #endo + mesoderm
      '#fff2cc',	#endoderm, yellow
      'lightgray' # NONE
      )
    )

pfla_ISH$Endoderm <- ifelse(is.na(pfla_ISH$Endoderm),0,ifelse(pfla_ISH$Endoderm == "",0,1))
pfla_ISH$Mesoderm <- ifelse(is.na(pfla_ISH$Mesoderm),0,ifelse(pfla_ISH$Mesoderm == "",0,1))
pfla_ISH$Ectoderm <- ifelse(is.na(pfla_ISH$Ectoderm),0,ifelse(pfla_ISH$Ectoderm == "",0,1))


pfla_ISH$germlayer <- 
  ifelse(
  pfla_ISH$Endoderm==1 & pfla_ISH$Mesoderm==1 & pfla_ISH$Ectoderm==1, "All",
  ifelse(
    pfla_ISH$Endoderm==1 & pfla_ISH$Mesoderm==1 & pfla_ISH$Ectoderm==0, "EndoMeso",
  ifelse(
    pfla_ISH$Endoderm==1 & pfla_ISH$Mesoderm==0 & pfla_ISH$Ectoderm==1, "EndoEcto",
  ifelse(
    pfla_ISH$Endoderm==0 & pfla_ISH$Mesoderm==1 & pfla_ISH$Ectoderm==1, "EctoMeso",
  ifelse(
    pfla_ISH$Endoderm==1 & pfla_ISH$Mesoderm==0 & pfla_ISH$Ectoderm==0, "Endoderm",
  ifelse(
    pfla_ISH$Endoderm==0 & pfla_ISH$Mesoderm==1 & pfla_ISH$Ectoderm==0, "Mesoderm",
  ifelse(
    pfla_ISH$Endoderm==0 & pfla_ISH$Mesoderm==0 & pfla_ISH$Ectoderm==1, "Ectoderm",
    "None"
  )))))))

pfla_ISH <- 
  merge(
    pfla_ISH,
    germlayers_colors,
    by.x="germlayer",by.y=1,all.x=T
    )[,c(13,1,14,2,3:12)]

pfla_ISH <-
  pfla_ISH[pfla_ISH$id != "",]

write.table(
  pfla_ISH,
  file = "outputs/functional_annotation/germlayers_ISH/in_situ_data.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
