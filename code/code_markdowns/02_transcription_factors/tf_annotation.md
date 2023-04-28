---
title: "Ptychodera Cisreg Development: Transcription Factor Annotation"
author: "Alberto Perez-Posada"
date: "3/25/2023"
output:
  html_document: default
  pdf_document: default
editor_options: 
  markdown: 
    wrap: sentence
---

## About

## Running Pfam Scan

```sh
PATH_TO_PFAM_DB=~/DATA/static/databases/pfam/pfam-A.hmm

PATH_TO_PROTEINS=./outputs/functional_annotation/transdecoder/pfla.transdecoder.pep

hmmscan --cpu 8 --noali -E 10e-3 --domE 10e-3 \
-o outputs/functional_annotation/TF_annotation/TFannot_hmmscan/pfla.pfam.outstd 
--tblout outputs/functional_annotation/TF_annotation/TFannot_hmmscan/pfla.pfam.tblout \
--domtblout outputs/functional_annotation/TF_annotation/TFannot_hmmscan/pfla.pfam.domtblout \
--pfamtblout outputs/functional_annotation/TF_annotation/TFannot_hmmscan/pfla.pfam.pfamtblout \
$PATH_TO_PFAM_DB $PATH_TO_PROTEINS
```

Subsetting for pfam TF domains:

```sh
grep -w -f assets/TF_domains_XG.txt outputs/functional_annotation/TF_annotation/TFannot_hmmscan/pfla.pfam.tblout | \
awk 'BEGIN{OFS="\t"} {print $3}' | \
sort | uniq \
> outputs/functional_annotation/TF_annotation/TFannot_hmmscan/pfla_pfam_tfs.txt
```


## Running OrthoFinder


```sh
orthofinder -t 12 -a 12 -X -f outputs/TF_annotation_TFannot_orthofinder/proteomes/
```

For simplicity, we will adapt the output to long format, and only keeping those orthogroups with genes from our species of interest.

```sh
awk '
BEGIN {FS="\t"; OFS="\t"}
NR==1 {for (i=1; i<=NF; i++) header[i]=$i}
NR>1 {
  for (i=2; i<=NF; i++) {
    split($i, values, ", ");
    for (j in values) {
      print $1, values[j]
    }
  }
}' outputs/functional_annotation/TF_annotation/TFannot_orthofinder/OrthoFinder/Results_*/Orthogroups/Orthogroups.tsv | \
awk 'BEGIN {OFS = "\t"} {if (NF == 2) print $0}' | awk '$2 ~ /^[A-Za-z0-9]+/' \
> outputs/functional_annotation/TF_annotation/TFannot_orthofinder/orthogroups_longformat.tsv
```

Since we work with pre-computed data for several of these species (amphioxus and sea urchin), we will translate their IDs to their equivalents in their gene expression reanalyses.

```python
# read the contents of dict file into a dictionary
with open("assets/species.dct", "r") as f:
    b_dict = dict(line.strip().split() for line in f)

# create an output file to write the translated values
with open("outputs/functional_annotation/TF_annotation/TFannot_orthofinder/orthogroups_longformat_translated.tsv", "w") as out:
    # iterate over the lines in gfams file
    with open("outputs/functional_annotation/TF_annotation/TFannot_orthofinder/orthogroups_longformat.tsv", "r") as f:
        for line in f:
            # split the line into two columns
            columns = line.strip().split("\t")
            # get the value from file "B" based on the key in column 2 of file "A"
            b_value = b_dict.get(columns[1], columns[1])
            # write the translated line to the output file
            out.write(columns[0] + "\t" + b_value + "\n")
```

We will later bring this file to R for parsing and annotating the TFs.

## Running ANANSE motif2factors

Following gimmemotifs documentation, we create a directory in genomepy's default directory with the name of our organism and a .pep file with the predicted protein sequences (from transdecoder).

After this, we can run motif2factors to transfer TF and motif annotation from well-known organisms to Ptychodera based on sequence homology (specifically, Orthofinder)
```sh
cd outputs/functional_annotation/TF_annotation/TFannot_motif2factors/

# Command for ananse motif2factors
gimme motif2factors --new-reference pfla --outdir mynewdatabase

# The version of ANANSE we used did not support (yet) underscores in name ids. Thus we remove them
perl -pi.bak -e "s/TCONS_/TCONS/" pfla.gimme.vertebrate.v5.0.motif2factors.txt
perl -pi.bak -e "s/TCONS_/TCONS/" pfla.gimme.vertebrate.v5.0.pfm

# Parse the association file and keep motif names as indication of TF class
cut -f1,2 pfla.gimme.vertebrate.v5.0.motif2factors.txt.bak | \
perl -pe "s/^([A-Za-z0-9]+\.){3}([A-Za-z0-9-_]+).*\t/\2\t/" | \
awk 'BEGIN{OFS="\t"} {print $2,$1}' > pfla_tfs_ananse.tsv
```

Since this is a parallel approach to our detection of TFs, we must take these into account when we assign a TF class to each Ptychodera TF gene.

## Download animalTFDB

We downloaded the identifiers for zebrafish and human transcription factors from the animalTFDB website.

We concatenated them into a single file that can work as a lookup table later.

```sh
cat Drer_tfs.tsv Hsap_tfs.tsv > assets/animaltfdb_drer_hsap.tsv
```

## Integrating the TF annotation

```r
# Load Orthogroups
gfams <-
  read.table(
    file = "outputs/functional_annotation/TF_annotation/TFannot_orthofinder/orthogroups_longformat_translated.tsv",
    header = FALSE,
    col.names = c("og","id")
  )

gfams$id <- sub("\\.p[1-9]+$","",gfams$id)

# Load AnimalTFDB database (Human and Zebrafish)
animaltfdb <-
  read.delim2(
    file = "assets/animaltfdb_drer_hsap.tsv",
    header = FALSE,
    col.names = c("id","gene","genename","class")
  ) [,c(1,4)]

# Create a dictionary with animaltfdb as key and motif2factor as values
dict_anim_m2f <- 
  split(
    read.table("assets/animaltfdb_ananse_equiv.tsv",header = T)[,2], # keys
    read.table("assets/animaltfdb_ananse_equiv.tsv",header = T)[,1] # values
    )

# Find the indices of the elements of my_list that match the values in class column
indices <- match(animaltfdb$class, names(dict_anim_m2f))

# Replace the values in class column with the corresponding values from animaltfdb/motif2factors dictionary
animaltfdb$class <- 
  sapply(
    animaltfdb$class,
    function(x) {
      ifelse(
        x %in% names(dict_anim_m2f),
        unlist(dict_anim_m2f[x]),
        x
        )
      }
    )


# Lookup table orthogroup-TFclass based on Hsap and Drer TFs
og_class <- 
  merge(
    gfams,
    animaltfdb,
    by.x = 2,
    by.y = 1,
    all.y = TRUE
    )[,c(2,3)]

og_class <-unique(og_class)

og_class <-
  og_class[complete.cases(og_class),]

# Identify OGs with more than 1 TFclass associated
dup_tfs <-
  names(
    rev(
      sort(
        table(
          og_class$og
        )
      )
    )
  )[1:(length(og_class$og)-length(unique(og_class$og)))]

# These ambiguous cases are treated as Others
og_class$class <-
  ifelse(
    og_class$og %in% dup_tfs,
    "Unknown",
    og_class$class
    )

og_class <- unique(og_class)


## Parse TF classes for all species

# ptychodera
pfla_tfs <-
  merge(
    gfams[grep("TCONS",gfams$id),],
    og_class,
    by.x = 1,
    by.y = 1,
  )[,c(2,3)]

# ptychodera is a special case. We must add motif2factors TFs

pfla_tfs_motif2factors <- read.table(
  file = "outputs/functional_annotation/TF_annotation/TFannot_motif2factors/pfla_tfs_motif2factors.tsv",
  header = TRUE,
  col.names = c("id","class")
)

pfla_tfs_motif2factors <- unique(pfla_tfs_motif2factors)

# We will add of course only the classes that are not already present
pfla_tfs <-
  rbind(
    pfla_tfs,
    pfla_tfs_motif2factors[
      which(!(pfla_tfs_motif2factors$id %in% pfla_tfs$id)),
      ]
  )

# amphioxus
blan_tfs <- 
  merge(
    gfams[grep("BL",gfams$id),],
    og_class,
    by.x = 1,
    by.y = 1
  )[,c(2,3)]

# sea urchin
spur_tfs <-
  merge(
    gfams[grep("STR",gfams$id),],
    og_class,
    by.x = 1,
    by.y = 1
  )[,c(2,3)]


# final output: longformat of the 4 species; gene <--> TF category
tfs_all_spp <- 
  cbind(
    rbind(pfla_tfs,blan_tfs,spur_tfs),
    mean = 0,
    counts = 0,
    cv = 0
  )

## Adding top classes of interest and color coding

topclasses <-
  read.table(
    file = "assets/top_tf_classes.txt"
  )[,1]

otherclasses <- 
  unique(pfla_tfs$class)[
    !(unique(pfla_tfs$class) %in% topclasses)
    ]

topclasses_col <- 
  setNames(
  c(
  "#e4eeb9", "#96e88e", "#50cc83", "#9ce5df", "#6dabd4",
  "#878cbc", "#5a5a82", "#2f5685", "#5a799a", "#8d909b",
  "#a4a3a2", "#a55297", "#f6afba", "#b92f7d", "#f38d97",
  "#e16265", "#ff9b7e", "#ffc694", "#ffebb5", "#ffd966"
  ),
  topclasses
  )

set.seed(121)
otherclasses_col <- 
  setNames(
  sample(colors(),length(otherclasses)),
  otherclasses
)

# Savint the data
save(
  pfla_tfs,
  blan_tfs,
  spur_tfs,
  topclasses,
  topclasses_col,
  otherclasses,
  otherclasses_col,
  file = "outputs/rda/TF_annotation.rda"
)

```
