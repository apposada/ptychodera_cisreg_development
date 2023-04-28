lookup_table <- list(
  "[0-9]\\.$" = "",
  "[0-9]\\(" = "\\(",
  "Hox..*" = "Hox",
  "AP-2..*" = "AP-2",
  "Fox..*" = "Fox",
  "CEBP..*" = "CEBP",
  "GATA" = "GATA",
  "DLX" = "Dlx",
  "PAX" = "Pax",
  "CUX" = "Cux",
  "LHX" = "Lhx",
  "Oct4..*" = "Oct"
)

tidyup_motifnames <- function(x,lookup_table){
  x <-
    ifelse(
      grepl("[Pp][-]*[567]3", x, perl = TRUE) == TRUE,
      "p53/63/73",
      gsub("[0-9]*$","",gsub("\\(..*","",x))
    ) # remove all the numerics at the end of names indicating actual gene (e.g. Sox21 becomes Sox ), but preserve p53/63/73
  x <- gsub("[0-9]\\.$","", x)
  x <- gsub("[0-9]\\(","\\(", x)
  x <- gsub("Hox..*","Hox", x)
  x <- gsub("AP-2..*","AP-2", x)
  x <- gsub("Fox..*","Fox", x,ignore.case=T)
  x <- gsub("CEBP..*","CEBP", x,ignore.case=T)
  x <- gsub("GATA","GATA", x,ignore.case=T)
  x <- gsub("DLX","Dlx", x,ignore.case=T)
  x <- gsub("PAX","Pax", x,ignore.case=T)
  x <- gsub("CUX","Cux", x,ignore.case=T)
  x <- gsub("LHX","Lhx", x,ignore.case=T)
  x <- gsub("Oct4..*","Oct", x,ignore.case=T)
}