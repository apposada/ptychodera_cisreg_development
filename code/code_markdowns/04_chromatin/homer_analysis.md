---
title: "Ptychodera Cisreg Development: HOMER analysis"
author: "Alberto Perez-Posada"
date: "3/31/2023"
output:
  html_document: default
  pdf_document: default
editor_options: 
  markdown: 
    wrap: sentence
---


## About

This markdowns gathers the main commands and instructions followed to perform Motif Enrichment Analyses using the HOMER suite.

## HOMER de novo motif finding

We ran De novo Motif Finding as follows:

```sh
#!/bin/bash

GENOME="~/projects/ptychodera_cisreg_development/data/refs/genome/ptyFlav3.fa"
ALL_PEAKS="~/projects/ptychodera_cisreg_development/data/ATAC/mapped/idr/pfla_all_peaks.bed"
OUTDIR="~/projects/ptychodera_cisreg_development/outputs/homer/denovo"

findMotifsGenome.pl $ALL_PEAKS $GENOME $OUTDIR -S 200 -p 12
```

## HOMER motif enrichment

For every set of peaks, we ran the `findMotifsGenome.pl` command. We always used the whole set of peaks of Ptychodera as background.

```sh
#!/bin/bash

GENOME="~/projects/ptychodera_cisreg_development/genome/ptyFlav3.fa"
PEAKS_DIRECTORY="~/projects/ptychodera_cisreg_development/outputs/atac/clusters_development/"
BG_PEAKS="~/projects/ptychodera_cisreg_development/data/ATAC/mapped/idr/pfla_all_peaks.bed"
OUTDIR="~/projects/ptychodera_cisreg_development/outputs/homer/atac_clusters/motif_enrichments/"

for i in ${PEAKS_DIRECTORY}/*.bed ; do
	
	x=${i##*/}
	z=${x%.sorted.bed}
	
	echo "treating sample $z"
	
	srun -c 12 -N 1 -n 1 -J homer_appos \
	findMotifsGenome.pl $i $GENOME \
	${OUTDIR}/homer_output_${z} -bg $BG_PEAKS \
	-p 12 -mset vertebrates

done

echo "Done."
```

Once it was over, we concatenated the results in a single table using:

```sh
#!/bin/bash

HOMERDIR="~/projects/ptychodera_cisreg_development/outputs/homer/atac_clusters/motif_enrichments/"

for i in ${HOMERDIR}/homer_output_pfla_ATAC_cluster_* ; do \
    x=${i##*/}
    awk -v var1=$x 'BEGIN {sple=var1} {print sple"\t"$0}' $i/knownResults.txt 
    done \
> ${HOMERDIR}/all_motifs_per_atac_cluster.tsv

```
