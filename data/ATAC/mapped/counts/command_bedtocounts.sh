#!/bin/bash
for i in /home/ska/aperpos/Def_Pfla/outputs/ATAC_ptyFlav3/devstages/*/*_nucfree.bed ; do
x=${i##*/}
z=${x%_nucfree.bed}
echo Starting intersect ... sample $x
mkdir ${z}
intersectBed -c -a /home/ska/aperpos/Def_Pfla/outputs/ATAC_ptyFlav3/devstages/idr/20200807_ptyFlav3_devstages_all_peaks_nonorm.bed -b $i -nonamecheck > ${z}/${z}_counts.txt
echo done. Starting parsing to bed... sample $x
awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"peak"NR,$4,"+"}' ${z}/${z}_counts.txt > ${z}/${z}.bed
echo done. Starting parsing to counts... sample $x
cut -f 4,5 ${z}/${z}.bed > ${z}/${z}.counts
rm ${z}/${z}_counts.txt
echo done with sample $x . Starting with another sample...
done
echo Done.
