---
title: "Ptychodera Cisreg Development: Mapping of ATAC-Seq data"
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

This markdowns gathers the main commands and instructions followed to map the ATAC-seq data of Ptychodera flava using a custom perl wrapper for bowtie2, and another custom bash wrapper of IDR analysis.

## Mapping using the internal pipeline of the lab

The internal pipeline consists of a custom perl wrapper that uses standard tools.

In essence, for each pair of FQ reads, the wrapper runs bowtie2 to map against the Ptychodera genome:

```perl
# [...]
  system ("echo bowtie2 -t --no-unal --no-mixed -X 2000 -p $cores -3 $trim3 -x $index -1 $fq1 -2 $fq2 > $log_file");
  system("bowtie2 -t --no-unal --no-mixed -X 2000 -p $cores -3 $trim3 -x $index -1 $fq1 -2 $fq2 2>> $log_file |
         samtools view -b - | samtools sort -m 5G -@ $sort_cores -T $temp - | samtools rmdup - $output 2> $tmp_log");
# [...]
```

After this, it transforms the bam to bed file while keeping only the aligned fragments of quality threshold above 10 and size above 130, to keep only the nucleosome-free regions.

```perl
# [...]
my $qual_thres = 10;    # quality threshold
my $min_length = 130;   # minimum length of fragments (130 for keeping only nucleosome free)
# [...]
while (my $read = $iter->next_seq) {
  my $len = $read->isize;
  unless(($read->qual >= $qual_thres) && (abs($len) <= $min_length) && (grep { $_ eq $read->flag } @flags)) {next};
# [...]
return $bedfile;
# [...]
```

Finally, it transforms the bed file to bigWig for visualisation in the genome browser, if needed:

```perl
# [...]
  system("bedSort $infile stdout | genomeCoverageBed -bg -i - -g $sizes | wigToBigWig stdin $sizes $file_out");
# [...]
```

It also creates a track file for genome visualisation of the nucfree regions:

```perl
# [...]
open(TRACK, ">", $track_file) or die "Could not open $track_file: $!\n";
    print TRACK "track type=bigWig name='$track_name' description='$track_name' color=0,0,0 visibility=2 maxHeightPixels=128:60:11 bigDataUrl=http://193.147.188.155/tracks/ATAC/$file_out\n";
    close TRACK;
# [..]
```
These files are stored and the .bed and bigwig files serve as input for a custom bash wrapper of the IDR.

Overall, the pipeline was run like this:

```sh
#!/bin/bash

genome_dir="~/projects/ptychodera_cisreg_development/data/refs/genome/"

ATAC_dir="~/projects/ptychodera_cisreg_development/data/ATAC/"

for i in ${ATAC_dir}/raw/* ; do

  x=${i##*/}
  
  echo "Treating sample ${x} ..."
  
  srun -n 1 -N 1 -c 12 -J atac_app \
    perl ~/ATAC_pipe.pl \
    -f1 ${i}/*_1.fq.gz -f2 ${i}/*_2.fq.gz \
    -o ${ATAC_dir}/mapped/nucf/${x} \
    -s ${genome_dir}/sizes.genome \
    -i ptyFlav3.fa.bowtie2 \
    -p 12 -bp $genome_dir \
    -ov_th
  
  echo "done ${x}"

done

echo "Done."
```

## Consensus peaks using the MACS2, IDR, and bedtools

### Important note

The code used for these steps is in the form of a bash wrapper was standardly utilised in the laboratory at the time of this project. Collaborative effort is paramount in the host laboratory for routinarily and reliably generate and analyse data in a systematic way. At the time of writing this it has been not possible for me, the first author, to identify the person and have confirmation of our rights to replicate this code here. For this reason I have adapted the minimum necessary information to showcase how the consensus peaks were generated. This is very likely to change in the future and I will make sure proper credit is given.

Peak calling using MACS2

```sh
# Replicate 1
macs2 callpeak -f BED -t ../1.bed --nomodel --extsize 100 --shift -45 --buffer-size 50000  -p $5 -g ${gsize} -n ${name_folder}_rep1 --outdir ./CONS_PEAKS &

# (Same is done for replicate 2)

# pool
macs2 callpeak -f BED -t ../1.bed ../2.bed --nomodel --extsize 100 --shift -45  --buffer-size 50000 -p $5 -g ${gsize} -n ${name_folder}_pooled_rep --outdir ./CONS_PEAKS &

wait %1 %2 %3
```

Top 500K peak extraction

```sh
# Replicate 1
sort -k8,8nr ./CONS_PEAKS/${name_folder}_rep1_peaks.narrowPeak | head -n 500000 > ./CONS_PEAKS/${name_folder}_rep1_toppeaks.bed &

# (Same is done for replicate 2)

# pool
sort -k8,8nr ./CONS_PEAKS/${name_folder}_pooled_rep_peaks.narrowPeak | head -n 500000 > ./CONS_PEAKS/${name_folder}_pool_toppeaks.bed &

wait %1 %2 %3
```

High-confidence peaks represent reproducible events and account for read sampling noise. The optimal set is more sensitive, especially when one of the replicates has lower data quality than the other.

For each biological replicate we generate two pseudorreplicates by random shuffling and splitting in two files. If a given region is reliably open chromatin-wise, splitting the reads in two should allow to reconstruct that peak in each separate pseudorreplicate.

Pseudorreplicate generation

```sh
# Replicate 1
nlines1=$( cat ../1.bed | wc -l ) 
halfnlines1=$(( (nlines1 + 1) / 2 )) 
cat ../1.bed | shuf | split -d -l ${halfnlines1} - ./OPT_PEAKS/${name_folder}_rep1_ps  ## -d indica que se pondran sufijos numericos

# (Same is done for replicate 2)

# Pool
cat ../1.bed ../2.bed | shuf > ./pooled.bed
nlinesp=$( cat ./pooled.bed | wc -l ) 
halfnlinesp=$(( (nlinesp + 1) / 2 )) 
cat ./pooled.bed | shuf | split -d -l ${halfnlinesp} - ./OPT_PEAKS/${name_folder}_pool_ps 
```

Pseudorreplicate peak calling

```sh
# Replicate 1 pseudorreplicates
macs2 callpeak -f BED -t ./OPT_PEAKS/${name_folder}_rep1_ps00 --nomodel --extsize 100 --shift -45  -p $5 -g ${gsize} -n ${name_folder}_rep1_ps00 --outdir ./OPT_PEAKS &
macs2 callpeak -f BED -t ./OPT_PEAKS/${name_folder}_rep1_ps01 --nomodel --extsize 100 --shift -45  -p $5 -g ${gsize} -n ${name_folder}_rep1_ps01 --outdir ./OPT_PEAKS &

# (Same is done for replicate 2)

# Pool pseudorreplicates
macs2 callpeak -f BED -t ./OPT_PEAKS/${name_folder}_pool_ps00 --nomodel --extsize 100 --shift -45  -p $5 -g ${gsize} -n ${name_folder}_pool_ps00 --outdir ./OPT_PEAKS &
macs2 callpeak -f BED -t ./OPT_PEAKS/${name_folder}_pool_ps01 --nomodel --extsize 100 --shift -45  -p $5 -g ${gsize} -n ${name_folder}_pool_ps01 --outdir ./OPT_PEAKS &

wait %1 %2 %3 %4 %5 %6

```

Top 500K peaks from each:


```sh
# Replicate 1 pseudorreplicates
sort -k8,8nr ./OPT_PEAKS/${name_folder}_rep1_ps00_peaks.narrowPeak | head -n 500000 > ./OPT_PEAKS/${name_folder}_rep1_ps00_toppeaks.bed &
sort -k8,8nr ./OPT_PEAKS/${name_folder}_rep1_ps01_peaks.narrowPeak | head -n 500000 > ./OPT_PEAKS/${name_folder}_rep1_ps01_toppeaks.bed &

# (Same is done for replicate 2)

# (Same is done for pseudorreplicates)

wait %1 %2 %3 %4 %5 %6
```

IDR analysis:

```sh
# Replicates 1,2, and pooled
idr --input-file-type narrowPeak --rank p.value --soft-idr-threshold 0.1 -s ./CONS_PEAKS/${name_folder}_rep1_toppeaks.bed ./CONS_PEAKS/${name_folder}_rep2_toppeaks.bed -p ./CONS_PEAKS/${name_folder}_pool_toppeaks.bed --verbose -l ./CONS_PEAKS/idr_consrv.log -o ./CONS_PEAKS/${name_folder}_rep1_rep2.txt --plot &

# Replicate 1 pseudorreplicates + Replicate 1
idr --input-file-type narrowPeak --rank p.value --soft-idr-threshold 0.1 -s ./OPT_PEAKS/${name_folder}_rep1_ps00_toppeaks.bed ./OPT_PEAKS/${name_folder}_rep1_ps01_toppeaks.bed -p ./CONS_PEAKS/${name_folder}_rep1_toppeaks.bed --verbose -l ./OPT_PEAKS/idr_opt_rep1.log -o ./OPT_PEAKS/${name_folder}_rep1_ps.txt --plot &

# (Same is done for replicate 2 pseudoreps and replicate 2)

# (Same is done for pool pseudoreps and pool)
wait %1 %2 %3 %4
```

Finally, some stats were generated using awk, including:

 - Number of conservative peaks
 - Number of filtered conservative peaks (GlobalIDR >=1) (Nt)
 - Number of optimal peaks for rep1
 - Number of filtered optimal peaks for rep1 (GlobalIDR >=1) (N1)
 - Number of intersected filtered optimal peaks for rep1 (GlobalIDR >=1)
 - Number of optimal peaks for rep2
 - Number of filtered optimal peaks for rep2 (GlobalIDR >=1) (N2)
 - Number of intersected filtered optimal peaks for rep2 (GlobalIDR >=1)
 - Number of optimal peaks for pooled data
 - Number of filtered optimal peaks for pooled data (GlobalIDR >=1) (Np)
 - Number of intersected filtered optimal peaks for pooled data (GlobalIDR >=1)
 - Rescue ratio = max (Np, Nt) / min (Np, Nt)
 - Self consistency ratio = max (N1, N2) / min (N1, N2)

Overall this bash wrapper was run as follows:

```sh
atac_path="~/data/ATAC/mapped/nucf/"

# Parameters in order:
# - sample id
# - genome size in Gb
# - replicate 1
# - p-value macs2

sbatch ~/jobs/idr_mod_job_custom.sh \
       EB \
       1.16e9 \
       ${atac_path}/EB-127-5/EB-127-5_nucfree.bed \
       ${atac_path}/EB-177-9/EB-177-9_nucfree.bed \
       10e-3

# (same for EG, MG, and LG)

```

## Generating a unique consensus peak file

```sh
cat *ConsPeaks.bed | \
sort -k1,1 -k2,2n | \
cut -f1,2,3,4,5,6 | \
mergeBed | sortBed | \
awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"peak"NR}' \
> pfla_all_peaks.bed
```

## Generating the counts file

```sh
atac_nucf_path="~/data/ATAC/mapped/nucf/"

for i in ${atac_nucf_path}/*/*_nucfree.bed ; do
  
  x=${i##*/}
  
  z=${x%_nucfree.bed}
  
  echo Starting intersect ... sample $x
  
  mkdir ${z}
  
  intersectBed -c -a pfla_all_peaks.bed -b $i -nonamecheck > ${z}/${z}_counts.txt
  
  echo done. Starting parsing to bed... sample $x
  
  awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"peak"NR,$4,"+"}' ${z}/${z}_counts.txt > ${z}/${z}.bed
  
  echo done. Starting parsing to counts... sample $x
  
  cut -f 4,5 ${z}/${z}.bed > ${z}/${z}.counts
  
  rm ${z}/${z}_counts.txt
  
  echo done with sample $x .

done

echo Done.
```

These .counts files can be passed to DESeq2 to analyse differential chromatin accessibility during development.