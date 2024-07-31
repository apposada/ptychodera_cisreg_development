df <- read.delim2("outputs/comparative/20240404_orthofinder/proteomes/OrthoFinder/Results_Apr05/Orthogroups/Orthogroups.tsv")
df_ <- data.frame(
  Orthogroup = character(),
  id = character()
)

for(i in 2:ncol(df)){
  d <- df[,c(1,i)]
  colnames(d) = c("Orthogroup","id")
  d <- d[d$id != "",]
  d <- as.data.table(d)
  
  # melt
  d <- d[, c(id=strsplit(id, ", ")), by=Orthogroup]
  
  # back to df
  d <- as.data.frame(d)
  
  df_ <- rbind(df_,d)
  
  message(i-1)
}
colnames(df_) <- c("gfam","id")
write.table(df_[,c(2,1)], "outputs/comparative/20240404_orthofinder/Orthogroups_melt.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
