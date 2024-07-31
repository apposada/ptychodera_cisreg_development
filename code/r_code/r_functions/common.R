# R common functions
fcha <- function(){ gsub("-","",Sys.Date()) }

#' Params:
#' a = a vector of sample names/condition names/etc ;
#' tbl = a feature x condition dataset, column names 
#' harmonised and parseable using vector 'a'
rowMeans_by_repl <- function(a,tbl){
  numsple <- length(grep(a,colnames(tbl)))
  if(numsple > 1 ){
    c <- rowMeans(tbl[,grep(a, colnames(tbl))])
  } else c <- tbl[,grep(a, colnames(tbl))]
  return(c)
}

devstages_ha_columns_simple <- function(x = NULL){
  library(ComplexHeatmap)
  library(circlize)
  
  if(is.null(x)){
    x = levels(condition_x)
  }
  
  ha <- HeatmapAnnotation(
    stage = anno_simple(
      x,
      col  = setNames(dev_palette,x),
      pch = as.character(0:(length(x)-1))
      ),
    show_annotation_name = FALSE
  )
  
  return(ha)
}

devstages_ha_columns <- function(x = NULL){
  library(ComplexHeatmap)
  library(circlize)
  
  if(is.null(x)){
    x = levels(condition_x)
  }
  
  ha <- HeatmapAnnotation(
    stage = x,
    col  = list(stage = setNames(dev_palette,x)),
    which = "column",
    show_legend = FALSE,
    show_annotation_name = FALSE
  )
  
  return(ha)
}

devstages_ha_rows_simple <- function(x = NULL){
  library(ComplexHeatmap)
  library(circlize)
  
  if(is.null(x)){
    x = levels(condition_x)
  }
  
  ha <- HeatmapAnnotation(
    stage = anno_simple(
      x,
      col  = setNames(dev_palette,x),
      pch = as.character(0:(length(x)-1))
    ),
    which = "row",
    show_legend = FALSE,
    show_annotation_name = FALSE
  )
  
  return(ha)
}


devstages_ha_rows <- function(x = NULL){
  library(ComplexHeatmap)
  library(circlize)
  
  if(is.null(x)){
    x = levels(condition_x)
  }
  
  ha <- HeatmapAnnotation(
    stage = x,
    col  = list(stage = setNames(dev_palette,x)),
    which = "row",
    show_legend = FALSE,
    show_annotation_name = FALSE
  )
  
  return(ha)
}

# quantile normalisation
# from: https://davetang.org/muse/2014/07/07/quantile-normalisation-in-r/
quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

tfs_ha_columns <- function(tfs){
  library(ComplexHeatmap)
  library(circlize)
  
  ha <- HeatmapAnnotation(
    class = tfs,
    col = 
      list(
        class = 
          c(topclasses_col,otherclasses_col)[
            names(c(topclasses_col,otherclasses_col)) %in%
            tfs
            ]
        ),
    show_annotation_name = FALSE
    )
  
  return(ha)
}

quick_ha <- function(x, col = "SunsetDark", rev = FALSE, side = "column"){
  require(ComplexHeatmap)
  require(circlize)
  
  l <- length(x)
  
  cols <- sequential_hcl(l,col,rev = rev)
  
  ha <-
    HeatmapAnnotation(
      stage = x,
      col = list(stage = setNames(cols,x) ),
      which = side,
      show_legend = FALSE,
      show_annotation_name = FALSE
    )
  
  return(ha)
}  