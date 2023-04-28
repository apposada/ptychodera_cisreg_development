---
title: "Ptychodera Cisreg Development: ANANSE analysis"
author: "Alberto Perez-Posada"
date: "4/24/2023"
output:
  html_document: default
  pdf_document: default
editor_options: 
  markdown: 
    wrap: sentence
---

## About

This markdowns gathers the main commands and instructions followed to run ANANSE using the RNA-seq and ATAC-seq data of Ptychodera.

More info on: https://gimmemotifs.readthedocs.io/en/master/reference.html#command-gimme-motif2factors

## Installation

These runs were performed using ANANSE developing branch in August/September 2021.

```sh
#install the virtual env
python3 -m venv ./ananse_pip_venv/
source ananse_pip_venv/bin/activate
pip install --upgrade pip
pip install cython
pip install pyarrow
pip install git+https://github.com/vanheeringen-lab/ANANSE.git@develop
python -m pip install 'fsspec>=0.3.3'
python -m pip install "dask[distributed]" --upgrade
```

## Running gimme motif2factors

Because we are working with a non-model organism, we had to adapt ANANSE's motif database to Ptychodera using motif2factors. At the time of running, this tool used Orthofinder to infer TF class and associated motifs on your organism of interest using a curated database of species. These were downloaded using genomepy.

```sh
#!/bin/bash
# downloaded Hsap, Mmus, ...  using genomepy
# added P.flava genome to the genomepy directory at ~/.local/share/genomes/PFLAV/ using a compliant format

#Finally ran motif2factors in the form of
gimme motif2factors --new-reference PFLAV --lenient --outdir 20210830_motif2factors
```

We also used the output TFs from this run to annotate more TFs (see the markdown of TF annotation.)

## Running ANANSE binding

To run ANANSE binding we did:

```sh
#!/bin/bash

M2F_DB="20210830_motif2factors/PFLAV.gimme.vertebrate.v5.0.pfm"
GENOME="~/projects/ptychodera_cisreg_development/data/refs/genome/ptyFlav3.fa"
ALL_PEAKS="~/projects/ptychodera_cisreg_development/data/ATAC/mapped/idr/pfla_all_peaks.bed"

ananse binding -A atac/EB-127-5.bam atac/EB-177-9.bam -g $GENOME -p $M2F_DB -o EB_binding -r $ALL_PEAKS -n 12
ananse binding -A atac/LG-124-6.bam atac/LG-130-7.bam -g $GENOME -p $M2F_DB -o EB_binding -r $ALL_PEAKS -n 12
```

## Running ANANSE network

To run ANANSE network we did:

```sh
#!/bin/bash

GENOME="~/projects/ptychodera_cisreg_development/data/refs/genome/ptyFlav3.fa"
ALL_PEAKS="~/projects/ptychodera_cisreg_development/data/ATAC/mapped/idr/pfla_all_peaks.bed"

ananse network -n 4 -e rna/02_EB_02.tsv -g $GENOME -a $ALL_PEAKS -o pfla_EB.network pfla_EB_binding/binding.h5
ananse network -n 4 -e rna/07_LG_01.tsv -g $GENOME -a $ALL_PEAKS -o pfla_LG.network pfla_LG_binding/binding.h5
```

The resulting networks were analysed in the markdown `ananse_graph_analysis.rmd`

## Running ANANSE influence

To run ANANSE influence we did:

```sh
# EB to LG
ananse influence -s pfla_EB.network -t pfla_LG.network -d 20220927_diff_RNA_EG_over_EB.tsv -o EG_EB_influence.txt -n 8 
```

