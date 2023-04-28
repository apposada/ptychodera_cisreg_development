---
title: "Ptychodera Cisreg Development: Mapping of RNA-Seq data"
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

This markdown gathers the instructions followed to map the RNA-seq data of Ptychodera using a transcriptome database of the longest isoforms per gene.

We first did a per-sample concatenation of libraries of different lanes using `cat` on the R1 gzipped files of ech sample, and hte same on the R2 gzipped files of each sample.

We then indexed the kallisto transcriptome.

```sh
kallisto index -i ~/projects/ptychodera_cisreg_development/genome/ptyFlav3.kallisto.index ~/projects/ptychodera_cisreg_development/genome/PO1410_Ptychodera_flava.transcript.fasta
```

Finally we mapped all samples using kallisto with standard parameters:

```sh
KALLISTO_INDEX="~/projects/ptychodera_cisreg_development/genome/ptyFlav3.kallisto.index"
KALLISTO_OUTPUT_FOLDER="~/projects/ptychodera_cisreg_development/data/RNA/kallisto/kallisto_raw_output"

for f in ~/projects/ptychodera_cisreg_development/data/RNA/raw/Sample* ; do
x=${f##*/}
echo Starting with sample ${x} ...
ls $f
echo "mkdir ${KALLISTO_OUTPUT_FOLDER}/kallisto_out_${x}"
mkdir ${KALLISTO_OUTPUT_FOLDER}/kallisto_out_${x}
echo "kallisto quant -t 10 -i $KALLISTO_INDEX $f/*R1*.f*.gz $f/*R2*.f*.gz -o ${KALLISTO_OUTPUT_FOLDER}/kallisto_out_${x}"
kallisto quant -t 10 -i $KALLISTO_INDEX $f/*R1*.f*.gz $f/*R2*.f*.gz -o ${KALLISTO_OUTPUT_FOLDER}/kallisto_out_${x}
echo "done with sample ${x} ..."
done
echo
echo "Done."
```

The resulting .h5ad files were also moved to a set of folders found at `data/RNA/kallisto/abundances/`. We also created a `sample_table.tsv` file in this path with the information on sample and batch for all the mapped libraries.