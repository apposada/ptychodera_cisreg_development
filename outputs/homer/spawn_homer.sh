#!/bin/bash

GENOME="/home/ska/aperpos/projects/ptychodera_cisreg_development/data/refs/genome/ptyFlav3.fa"
PEAKS_DIRECTORY=$(realpath $1)
BG_PEAKS_ALL=$(realpath $2)
OUTDIR=$(realpath $3)

mkdir -p ${OUTDIR}/backgrounds/

for i in ${PEAKS_DIRECTORY}/*.bed ; do
	
	x=${i##*/}
	z=${x%.bed}
	
	echo "treating sample $z"
	
	findMotifsGenome.pl $i $GENOME \
	${OUTDIR}/homer_output_${z} -bg $BG_PEAKS_ALL \
	-p 12 -mset vertebrates -size 200

done
