library(XML)
library(plyr)
library(methods)
options(scipen=999)

source("~/projects/ptychodera_cisreg_development/code/r_code/r_functions/sourcefolder.r")
sourceFolder("~/projects/ptychodera_cisreg_development/code/r_code/r_functions",recursive = TRUE)

# Load orthoxml HierarchicalGroups

hogs <- 
  xmlParse(
    file = "/home/ska/aperpos/projects/ptychodera_cisreg_development/outputs/comparative/20240404_oma/precompiled/Output/HierarchicalGroups.orthoxml",
    useInternalNodes = TRUE
  )

# Convert xml doc into List

xL <- xmlToList(hogs) 

xL_speciesgeneids <- 
  unlist(xL[1:(length(xL)-2)]) # if this works then it can be extended to any oma result

xL_speciesgeneids <- 
  sub(" .*","",xL_speciesgeneids)

# Generate the dictionary 'integer id --> geneid'

speciesgeneids <- 
  data.frame(
    id = xL_speciesgeneids[
      seq(1,length(xL_speciesgeneids),by=2)
    ],
    geneid = xL_speciesgeneids[
      seq(2,length(xL_speciesgeneids),by=2)
    ]
  )

speciesgeneids <- 
  speciesgeneids[speciesgeneids$id != "", ]

speciesgeneids <- 
  speciesgeneids[speciesgeneids$geneid != "0", ]

write.table(
  speciesgeneids,
  "/home/ska/aperpos/projects/ptychodera_cisreg_development/outputs/comparative/20240404_oma/dictionary_oma_integer_to_omaid.tsv",
  sep="\t",
  dec=".",
  quote=F,
  row.names=F,
  col.names=F
)

xL_groups <- xL[[(length(xL)-1)]] # xL[[length(xL-1)]]


library(stringr)

simpleogrouplist <- list()

for( ogroup in xL_groups ){
  
  a <- ogroup
  
  newname <- 
    unlist(a)[
      which(
        str_detect( names(unlist(a)) , "attrs" ) # What does this do?
      )
    ]
  
  newogroup <- 
    list(
      c( unlist(a)[
        which(
          str_detect( names(unlist(a)) , "geneRef")
        )
      ]
      )
    )
  
  names(newogroup) <- newname
  
  simpleogrouplist <- append(simpleogrouplist,newogroup)
}

og_L <- simpleogrouplist
num_digits <- magn_order(length(og_L))+1
names(og_L) <- paste0("HOG",formatC(1:length(og_L),width=num_digits,flag="0"))
data <- ldply(og_L, data.frame)

colnames(data) <- c("gfam","oma_integer")

write.table(
  data[,c(2,1)], "/home/ska/aperpos/projects/ptychodera_cisreg_development/outputs/comparative/20240404_oma/oma_integer_gfam.tsv",
  sep="\t",
  dec=".",
  quote=F,
  row.names=F
)

data$omaid <- 
  translate_ids(
    data$oma_integer,
    dict = speciesgeneids
  )

write.table(
  data[,c(3,1)],
  "/home/ska/aperpos/projects/ptychodera_cisreg_development/outputs/comparative/20240404_oma/oma_omaid_gfam.tsv",
  sep="\t",
  dec=".",
  quote=F,
  row.names=F
)
