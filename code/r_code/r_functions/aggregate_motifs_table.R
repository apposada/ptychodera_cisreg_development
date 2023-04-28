aggregate_motifs_table <- function(tsv, fun = max){
  # took the idea from this R API for Homer : https://robertamezquita.github.io/marge/articles/marge-workflow.html
  aggreg <- 
    aggregate(
      tsv[,c(8,11,15)],
      by = list(
        class = tsv$class,
        category = tsv[[1]]),
      FUN = fun
    ) 
  
  colnames(aggreg) <- 
    c("class","category","pct","logqval","size_pct")
  
  return(aggreg)
}
