#!/bin/bash

genome_dir="/home/ska/aperpos/projects/ptychodera_cisreg_development/data/refs/genome"

ATAC_dir="/home/ska/aperpos/projects/ptychodera_cisreg_development/data/ATAC/"

for i in /home/ska/aperpos/projects/ptychodera_cisreg_development/data/ATAC/raw/* ; do

  x=${i##*/}
  
  echo "Treating sample ${x} ..."
  
  echo "srun -n 1 -N 1 -c 12 -J atac_app \
    perl ~/ATAC_pipe.pl \
    -f1 ${i}/*_1.fq.gz -f2 ${i}/*_2.fq.gz \
    -o ${ATAC_dir}/mapped/nucf/${x} \
    -s ${genome_dir}/sizes.genome \
    -i ptyFlav3.fa.bowtie2 \
    -p 12 -bp $genome_dir \
    -ov_th"
    
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
