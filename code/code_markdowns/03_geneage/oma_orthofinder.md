---
title: "Ptychodera Cisreg Development: Orthology Annotation"
author: "Alberto Perez-Posada"
date: "3/30/2023"
output:
  html_document: default
  pdf_document: default
editor_options: 
  markdown: 
    wrap: sentence
---

## About

In this notebook we provide the instructions we followed to generate orthology data for Ptychodera flava.

## Running OMA

We downloaded a bundle of the OMA standalone software together with pre-computed data for a number of species (see article).

(species tree here: `"(Yeast,(Capsaspora,(Nematostella,((Schmidtea,(C. elegans,Drosophila)),(((((Mouse,Human,Xenopus),(Medaka,(Zebrafish,Cavefish))),Amphioxus),((Ptychodera,Saccoglossus),(Anneissia japonica,(Amphiura filiformis,((Apostychopus californicum,Apostychopus japonicus),(Lytechinus variegatus,(Heliocidaris sp,Strongylocentrotus purpuratus))))))))))));"` )

After this, we added several of our species of interest, including aneissia japonica, strongylocentrotus purpuratus, and ptychodera flava.

Since Amphioxus was already present in our data albeit from a different version, we would have to convert between the sequence identifiers down the line.

```sh
#!/bin/bash
#SBATCH -J OMA_aperpos
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
export NR_PROCESSES=1
~/projects/ptychodera_cisreg_development/outputs/comparative/oma/OMA_precompiled/bin/oma
```

### Extracting one to one orthologs

We can do that with:

```sh
# Ptychodera and Amphioxus
grep "1:1" OMA_precompiled/Output/PairwiseOrthologs/STRPU_ok-PFLAV.txt | \
awk '{OFS = "\t"} {print $4,$3}' | \
perl -pe "s/\.p[0-9]+//" > 1to1_pfla_spur.tsv 

# Ptychodera and Sea Urchin
grep "1:1" OMA_precompiled/Output/PairwiseOrthologs/PFLAV-BRALA.txt | \
awk '{OFS = "\t"} {print $4,$3}' | \
perl -pe "s/\.p[0-9]+//" > 1to1_pfla_blan.tsv 
```

We will also isolate the ptychodera genes and their associated oma group (a proxy for "gene family") that will come in handy to infer the transcriptomic age of the embryonic development of this species.

```sh
grep PFLAV OrthologousGroups.txt | \
perl -pe "s/\t/\n/g" | \
grep -P "PFLAV|OMA[0-9]{2,}" | \
paste - - | \
perl -pe "s/PFLAV\://" | \
perl -pe "s/\.p[0-9]+//"  > pfla_omanames.tsv
```

Lastly, we will adapt OrthologousMatrix's format to account for the extra column with the OMA gfam identifiers:

```sh
awk 'BEGIN {print "OGroup"} !/^#/ {print $1}' ../oma/OrthologousGroups.txt | \
paste - <(awk '!/^#/ {print substr($0, index($0,$2))}' ../oma/OrthologousMatrix.txt ) \
> OrthologousMatrix.txt.2
```

## Extracting gene family information

We aim at obtaining a 1gene - 1genefamily long format with every gene in the comparative dataset and the Orthogroup to which they are annotated.

For this we loaded in R the Hierachical Orthogroups OrthoXML file that is returned by OMA and we run the following:

```r
library(XML)
library(plyr)
library(methods)

hogs <- 
xmlParse(
  file = "/home/ska/aperpos/Def_Pfla/outputs/oma/20210821_OMA.2.4.1_ECHI_HEMI/Output/HierarchicalGroups.orthoxml", # change to new path
  useInternalNodes = TRUE
  )

# Convert xml doc into List
xL <- 
  xmlToList(hogs)

xL_speciesgeneids <- unlist(xL[1:lengt(xL)-2]) # this way it can be extended to any oma run
xL_speciesgeneids <- sub(" .*","",xL_speciesgeneids)

# Generate the dictionary 'integer id --> geneid'
speciesgeneids <- 
  data.frame(
    id = xL_speciesgeneids[seq(1,length(xL_speciesgeneids),by=2)],
    geneid = xL_speciesgeneids[
      seq(2,length(xL_speciesgeneids),by=2)
    ]
  ) 
speciesgeneids <- speciesgeneids[speciesgeneids$id != "", ]
speciesgeneids <- speciesgeneids[speciesgeneids$geneid != "0", ]

write.table(
  speciesgeneids,
  file = "/home/ska/aperpos/Def_Pfla/outputs/oma/20210821_OMA.2.4.1_ECHI_HEMI/20210823_results/20210824_dictionary_oma_integer_to_omaid.tsv", # change to new path
  sep = "\t",
  dec = ".",
  quote = F,
  row.names = F,
  col.names = F
  )

xL_groups <- xL[[length(xL-1)]]

# Create a list of integer id and gfam
library(stringr)
simpleogrouplist <- list()
for( ogroup in xL_groups ){
  a <- ogroup
  newname <- 
    unlist(a)[
      which(str_detect(names(unlist(a)),"attrs"))
      ]

  newogroup <- 
    list(
      c(
        unlist(a)[
          which(str_detect(names(unlist(a)),"geneRef"))
          ]
      )
    )
  names(newogroup) <- newname

  simpleogrouplist <- append(simpleogrouplist,newogroup)
}

data <- ldply(simpleogrouplist, data.frame)
colnames(data) <- c("gfam","id")

data$gfam <- as.integer(data$gfam) # actually... these are the longformat 1gene-1family. I would just need to retrieve those from pfla and seaurchin, and there it is.

write.table(
  data,
  "/home/ska/aperpos/Def_Pfla/outputs/oma/20210821_OMA.2.4.1_ECHI_HEMI/20210823_results/20210824_gfams_integer.tsv", # change to new path
  sep = "\t",
  dec = ".",
  quote = F,
  row.names = F
  )

```

## Running Count Software

From the OMA output, we can use the OrthologousMatrix.txt file to quantify the presence and composition of gene families in each species across our species tree, in turn figuring out the gene age of each gene in our species of interest. For this, we can use the Count software. Since it has a graphical interface, we will simply detail our steps here.

[...]

This output can enter directly into the markdown of Gene Age.

## Running OrthoFinder

For OrthoFinder, which we used also as an alternative way to estimate conservation of gene families throughout development, we ran orthofinder with standard parameters and the same set of proteomes we used for OMA.

```sh

```