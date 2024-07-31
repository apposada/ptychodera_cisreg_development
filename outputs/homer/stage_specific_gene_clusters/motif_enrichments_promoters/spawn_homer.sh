#!/bin/bash

GENOME="/home/ska/aperpos/projects/ptychodera_cisreg_development/data/refs/genome/ptyFlav3.fa"
PEAKS_DIRECTORY="/home/ska/aperpos/projects/ptychodera_cisreg_development/outputs/homer/stage_specific_gene_clusters/motif_enrichments_promoters/"
BG_PEAKS="/home/ska/aperpos/projects/ptychodera_cisreg_development/outputs/functional_annotation/promoters/ptyFlav3_CYi_longest.promoters.bed"
OUTDIR="/home/ska/aperpos/projects/ptychodera_cisreg_development/outputs/homer/stage_specific_gene_clusters/motif_enrichments_promoters/homer_outputs/"

mkdir -p $OUTDIR

for i in ${PEAKS_DIRECTORY}/*.bed ; do

        x=${i##*/}
        z=${x%.bed}

        echo "treating sample $z"

        srun -c 12 -N 1 -n 1 -J homer_appos \
        findMotifsGenome.pl $i $GENOME \
        ${OUTDIR}/homer_output_${z} -bg $BG_PEAKS \
        -p 12 -mset vertebrates

done

echo "Done."
