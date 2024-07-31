#!/bin/bash
# from https://www.baeldung.com/linux/variable-preserve-linebreaks
RES_DIR=$1
COLUMN=$2
RES_FILES=$(find ${RES_DIR} -name "knownResults.txt" | sort | uniq)
head -n1 $(echo "$RES_FILES" | head -n1) | \
     grep Benjamini | \
     sort | uniq | \
     head -n1 | perl -pe "s/Motif\([of 0-9]+\)/Motif/g" | \
     awk 'BEGIN {OFS="\t"} {print "module",$0}' > ${RES_DIR}/motifs_all.tsv

for i in ${RES_FILES} ; do
    x=`echo $i | cut -d "/" -f2 | cut -d "_" -f $COLUMN `
    awk -v var1=$x 'BEGIN {OFS="\t"} {print var1,$0}' $i | tail -n +2  >> ${RES_DIR}/motifs_all.tsv
done

