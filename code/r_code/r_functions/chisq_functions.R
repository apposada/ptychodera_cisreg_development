chisq_and_posthoc <- function(tbl,method_ph="BH"){
  require(rstatix)
  require(chisq.posthoc.test)
  require(dplyr)
  require(tidyverse)
  require(tidyr)
  
  # chi squared
  chsq <- chisq_test(tbl)
  
  sample_cols <- ncol(tbl)
  
  # posthoc test
  ph <- chisq.posthoc.test(tbl, method = method_ph) %>% 
    pivot_longer(2+1:sample_cols,names_to = "category", values_to = "estimate") %>%
    dplyr::rename("Statistic"=Value, "Value"=estimate, "Cluster" = Dimension) %>%
    dplyr::select(Cluster,category, Statistic,Value)
  
  ## Aggregate at the cluster, category level(s) so pvalues and residuals are not in separate rows
  ph_tidy <- 
    left_join(
      ph[ph$Statistic=="p values",], 
      ph[ph$Statistic=="Residuals",],
      by = c("Cluster","category")
    ) %>%
    dplyr::rename("pval"=Value.x, "res"= Value.y) %>%
    dplyr::select(Cluster, category, pval, res) %>%
    arrange(Cluster, desc(category)) %>% 
    mutate(
      significance = cut(
        pval,
        breaks = c(0,0.001,0.01,0.05,1),
        labels = c("***","**","*"," "), include.lowest = T
      )
    )
  
  res <- 
    list(
      chi_squared_test = chsq,
      posthoc_test = ph,
      result = ph_tidy
    )
  
  return(res)
}

make_chisq_heatmap <- function(x, reorder_values = NULL){
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  require(ComplexHeatmap)
  require(circlize)
  
  # Handle the matrices
  ph_res <- x$result[,c(1,2,4)] %>% pivot_wider(names_from = Cluster, values_from = res) %>% column_to_rownames("category") %>% as.matrix()
  ph_pva <- x$result[,c(1,2,3)] %>% pivot_wider(names_from = Cluster, values_from = pval) %>% column_to_rownames("category") %>% as.matrix()
  
  if (!is.null(reorder_values)){
    ph_res <- ph_res[,match(levels(reorder_values),colnames(ph_res))]
    ph_pva <- ph_pva[,match(levels(reorder_values),colnames(ph_pva))]
  }
  
  # stuff for heatmap
  col_enr <- c("#90529f","white","#45aaa7")
  col_enr_ramp <- 
    c(
      colorRampPalette(c(col_enr[1],col_enr[2]))(20)[1:19],
      colorRampPalette(c(col_enr[2],col_enr[3]))(21)
    )
  enr_pal_broad_1 <- circlize::colorRamp2(
    breaks = c(seq(min(ph_res),-0.1,length = 19),0,seq(0.1,max(ph_res),length = 20)), # this is NOT perfect. should be centered around zero no matter what. How to do this??
    colors = col_enr_ramp
  )
  pval_fun <- function(j,i,x,y,width,height,fill){
    if(ph_pva[i,j] < 0.01){
      if(ph_pva[i,j] < 0.001){
        grid.text("**", x, y, gp = gpar(fontsize = 10))
      } else{
        grid.text("*", x, y, gp = gpar(fontsize = 10))
      }
    }
  }
  
  # Make the heatmap
  ph_hm <-
    ComplexHeatmap::Heatmap(
      name = "residuals",
      ph_res,
      col = enr_pal_broad_1,
      cluster_columns = FALSE,
      cluster_rows = FALSE,
      column_names_side = "top",
      column_names_gp = gpar(fontsize = 6),
      row_names_side = "left",
      cell_fun = pval_fun
    )
  
  # Gather output
  res <- list(
    residuals_matrix = ph_res,
    pvalues_matrix = ph_pva,
    heatmap = ph_hm
  )
  
  # return
  return(res)
}

make_enr_pal <- function(values, pal = c("#90529f","white","#45aaa7"), q = .99){
  require(circlize)
  
  ramp = 
    c(
      colorRampPalette(c(pal[1],pal[2]))(20)[1:19],
      colorRampPalette(c(pal[2],pal[3]))(21)
    )
  
  lower = quantile(values, 1-q) # min(values) # should be `quantile(values, .01)` ?
  upper = quantile(values, q)
  
  breaks = c(
    seq(lower, -.1, length = 19),
    0,
    seq( .1, upper, length = 20)
  )
  
  enr_pal = 
    circlize::colorRamp2(
      breaks = breaks,
      colors = ramp
    )
  
  return(enr_pal)
}





