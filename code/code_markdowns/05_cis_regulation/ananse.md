---
title: "Ptychodera Cisreg Development: running ANANSE"
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

## Preparation

### Installing genomes

```sh
activate_miniconda3_venv
mamba create -n ananse_venv -c bioconda -c conda-forge ananse #agat too?
conda activate ananse_venv
genomepy install -a GRCh38.p13 -p Ensembl
genomepy install -a GRCm38.p6 -p Ensembl
genomepy install -a BraLan2 -p Ensembl
```

### Preparing the Pfla genome

```sh
cd ~/.local/share/genomes/
mkdir -p PFLAV
cd ./PFLAV/
cp ${PATH_TO_GENOME} ./PFLAV.fa
cp ${PATH_TO_PROTEOME} ./PFLAV.pep.fa
cp ${PATH_TO_GTF} ./PFLAV.annotation.gtf
cp ${PATH_TO_BED12} ./PFLAV.annotation.bed12.bed
```

We do the same with a number of genomes of different species: 
Ofus Hmia Aque Cowc Djap BraLan2 Bflo STRPU LYTVA SACKO DROME NEMVE


We will also create a folder to store everything related to ANANSE:

```sh
mkdir -p ~/projects/ptychodera_cisreg_development/outputs/ananse/
cd ~/projects/ptychodera_cisreg_development/outputs/ananse/
mkdir -p data/{atac,rna} prep outs/{binding,influence,network}
```
### Running Motif2Factors

see gimmemotifs.readthedocs.io

```sh
gimme motif2factors --new-reference PFLAV --lenient --tmpdir ./tmp_m2f/ --ortholog-reference Ofus Hmia Aque Cowc Djap BraLan2 Bflo STRPU LYTVA SACKO DROME NEMVE --threads 12 --outdir ./data/m2f/pfla_gimme_m2f
```

### get peak summits & max read depths

We will center the peak coordinates around the summit. For that we will do a pileup and a 

```python
import pandas as pd
import numpy as np
import subprocess as sp
import os	

peakfile = "~/projects/ptychodera_cisreg_development/data/ATAC/mapped/idr/pfla_all_peaks.bed"
df = pd.read_table(peakfile, comment="#", header=None)
df.columns = ["chrom", "start", "end", "name"]
df

# no comments/headers
df.to_csv("/mnt/sda/alberto/projects/smed_cisreg/outputs/ananse/prep/peaks.bed", sep="\t", index=False, header=False)
```

We run pileup using a 

```sh
bedtools slop -b 1000 -i ~/projects/ptychodera_cisreg_development/data/ATAC/mapped/idr/pfla_all_peaks.bed -g ~/projects/ptychodera_cisreg_development/data/refs/genome/sizes.genome > peaks_extended.bed

samtools merge -@12 all_reads.bam ../data/atac/*.bam

samtools merge ~/projects/ptychodera_cisreg_development/data/ATAC/mapped/*.bam all_reads.bam
samtools mpileup -A -B -l peaks_extended.bed -o peak_extended_pileup.tsv all_reads.bam
```

```python
import pandas as pd
import numpy as np
import subprocess as sp
import os	
import matplotlib

peakfile = "~/projects/ptychodera_cisreg_development/outputs/ananse/prep/peaks_extended.bed"
df = pd.read_table(peakfile, comment="#", header=None)
df.rename(index = {0 : "chrom", 1 : "start" , 2 : "end"})

pu = pd.read_table(
    "~/projects/ptychodera_cisreg_development/outputs/ananse/prep/peak_extended_pileup.tsv", 
    header=None,
    usecols=[0,1,3],
    names=["chrom", "pos", "count"],
)

summits = []
depths = []
last_chrom = ""
for idx, (chrom, start, end) in df.iterrows():
    if last_chrom != chrom:
        print("chrom:", chrom)  # print progress
        chrom_pu = pu[pu["chrom"] == chrom]
        last_chrom = chrom
    
    sub_pu = chrom_pu[chrom_pu["pos"].between(start, end)]
    maxima = sub_pu["count"] == sub_pu["count"].max()
    if sum(maxima) == 0:
        summit = None
    elif sum(maxima) == 1:
        summit = int(sub_pu[maxima]["pos"])
    else:
        idx = sub_pu[maxima].index
        middle_summit = idx[len(idx) // 2]
        summit = int(sub_pu.loc[middle_summit]["pos"])
    depth = sub_pu["count"].max()
    
    depths.append(depth)
    summits.append(summit)



df2 = df.copy()
df2["summit"] = summits
df2["depth"] = depths
df2 = df2.dropna()
df2["summit"] = df2["summit"].astype(np.int64)
df2["depth"] = df2["depth"].astype(np.int64)

depth = list(range(0, df2["depth"].max()))
n_regions = []
for n in depth:
    nr = sum(df2['depth'] == n)
    n_regions.append(nr)
    if n <= 10:
        print(
            f"regions with depth {n}: {nr}"
        )



density_readsperpeak = pd.DataFrame({"reads": depth, "number of reads under peaks": n_regions}).plot(x="reads", xlim=(0, 500))

density_readsperpeak.get_figure().savefig('density_readsperpeak.pdf', format='pdf')
```



```python
min_depth = list(range(0, df2["depth"].max()))
n_regions = []
for n in min_depth:
    nr = sum(df2['depth'] >= n)
    n_regions.append(nr)

density_mindepth = pd.DataFrame({"min_depth": min_depth, "remaining_regions": n_regions}).plot(x="min_depth", xlim=(0, 500))#, ylim=(0, 20_000))

density_mindepth.get_figure().savefig('density_mindepth.pdf', format='pdf')
```

Filter peaks by min read depth
```python
min_reads = 75

df3 = df2.copy()
df3.rename(columns = {0 : "chrom", 1 : "start" , 2 : "end"}, inplace = True)
df3 = df3[df3['depth'] >= min_reads]

# standardize the peak widths
df3["start"] = df3["summit"]-100
df3["end"] = df3["summit"]+100
df3 = df3[["chrom", "start", "end"]]
df3
```

```python
df3.to_csv(f"/home/ska/aperpos/projects/ptychodera_cisreg_development/outputs/ananse/prep/peaks_normalized_min{min_reads}.bed", index=False, header=False, sep="\t")
```

## Running ANANSE standalone (one by one)

### ANANSE binding

```sh
cd ~/projects/ptychodera_cisreg_development/outputs/ananse
activate_miniconda3_venv
conda activate ananse_venv

M2F_DB="data/m2f/pfla_gimme_m2f/PFLAV.gimme.vertebrate.v5.0.pfm"
GENOME="~/.local/share/genomes/PFLAV/PFLAV.fa"
ALL_PEAKS="~/projects/ptychodera_cisreg_development/outputs/ananse/prep/peaks_normalized_min75_corrected.bed"

mkdir -p outs/binding/

# EB
samtools merge -@12 data/atac/EB.bam data/atac/EB-*.bam
samtools index -@12 EB.bam
ananse binding -A data/atac/EB.bam -g $GENOME -p $M2F_DB -o outs/binding/EB -r $ALL_PEAKS -n 12
# LG
samtools merge -@12 data/atac/LG.bam data/atac/LG-*.bam
samtools index -@12 LG.bam
ananse binding -A data/atac/LG.bam -g $GENOME -p $M2F_DB -o outs/binding/LG -r $ALL_PEAKS -n 12
```

### ANANSE network

We retrieve the RNA data and export it for ANANSE

```r
load("~/projects/ptychodera_cisreg_development/outputs/rda/deseq2.rda")
write.table(
  data.frame(tpm=rowMeans(vsd_allgen[,c(4,5)])/sum(rowMeans(vsd_allgen[,c(4,5)]))*10e6),
  file = "~/projects/ptychodera_cisreg_development/outputs/ananse/data/rna/EB.tsv", sep = "\t",
  quote = FALSE
  )
write.table(
  data.frame(tpm=rowMeans(vsd_allgen[,c(16,17)])/sum(rowMeans(vsd_allgen[,c(16,17)]))*10e6),
  file = "~/projects/ptychodera_cisreg_development/outputs/ananse/data/rna/LG.tsv", sep = "\t",
  quote = FALSE
)
```

```sh
mkdir -p outs/network
#cat ../../data/refs/genome/ptyFlav3_CYi_longest.bed12.bed | perl -pe "s/__length[0-9_]+pilon//" | perl -pe "s/__/_/" | sortBed > ~/.local/share/genomes/PFLAV/PFLAV.annotation.bed
GENOME="~/.local/share/genomes/PFLAV/PFLAV.fa"
ANNOTATION="~/.local/share/genomes/PFLAV/PFLAV.annotation.bed12.bed"
ALL_PEAKS="~/projects/ptychodera_cisreg_development/outputs/ananse/prep/peaks_normalized_min75_corrected.bed"
# EB
ananse network -n 12 -e data/rna/EB.tsv -g $GENOME -a $ANNOTATION --full-output -o outs/network/EB.network outs/binding/EB/binding.h5
# LG
ananse network -n 24 -e data/rna/LG.tsv -g $GENOME -a $ANNOTATION --full-output -o outs/network/LG.network outs/binding/LG/binding.h5
```

### ANANSE influence

```sh
mkdir -p outs/influence
# EB to LG
ananse influence -s outs/network/EB.network -t outs/network/LG.network -d data/diff_RNA_LG_EB.tsv -o outs/influence/EG_EB_influence.txt -n 8 -i 100000
```