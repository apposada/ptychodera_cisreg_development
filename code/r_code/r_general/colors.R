library(RColorBrewer)
library(viridis)
library("colorspace")
library(ggpubfigs)
library(circlize)
library(ComplexHeatmap)

dev_palette <- 
  c(
    "#ec7c84",
    "#ec7c84",
    "#efaa90",
    "#efaa90",
    "#fbcf99",
    "#fbcf99",
    "#fbcf99",
    "#fbcf99",
    "#f2f7a5",
    "#ceeaa7",
    "#a7deaa",
    "#69cbae",
    "#4ac1b0",
    "#4c959c",
    "#2B6B8C",
    "#2B6B8C"
  )

pal <- friendly_pal("contrast_three", 10, type = "continuous")

app_pal1=colorRampPalette(
  rev(c(
    "#FFEA46",
    "#E5B43F",
    "#A07D2C",
    "#213351", # try this 142c4b
    "#458193",
    "#54958D",
    "#7FB7C2"
  ))) (10)



app_pal2 =colorRampPalette(
  rev(c(
    "#FFEA46",
    "#E5B43F",
    "#A07D2C",
    "#404e65",
    "#465D7E",
    "#587B98",
    "#8AA2B6"
  ))) (10)


#Color palette for heatmap of TF-TF Spearman correlation
tfs_sp_hm_col <- colorRamp2(
  c( # breaks, clipped
    seq(-1,0,len=10),
    seq(0,1,len=10)
  ),
  colorRampPalette( # colors
    rev(brewer.pal(7,"RdYlBu"))
  )(20)
)

# Color palette for heatmap of expression
tfs_expr_hm_col <- colorRamp2(
  c( # breaks, clipped
    seq(0,4,len=10)
  ),
  rev(sequential_hcl(10,"YlGnBu")) # colors
)


viridis_pastel <-
  c(
    "#ffee61", "#96e88e", "#5dc9ac",
    "#4da2ba", "#6b6eab", "#552761"
  )
  
tf_clusters_hm_col <- 
  colorRamp2(
    seq(0,0.5,len=20),
    colorRampPalette(
      rev(viridis_pastel)
      )(20)
    )

gb_col <- 
  c(
    yellow1 = "#f5f8aa",
    yellow2 = "#d0be7b",
    orange1 = "#efad84",
    orange2 = "#d78c91",
    purple1 = "#9d95e2",
    purple2 = "#7b6b94",
    blue1 = "#92baef",
    blue2 = "#4a42b5",
    teal1 = "#6bad94",
    teal2 = "#2f738a",
    green1 = "#9cce7b",
    green2 = "#4a9c5a",
    white1 = "#eff0f9",
    gray1 = "#9696a4",
    gray2 = "#5a5a63",
    black1 = "#1a1a1a"
  )