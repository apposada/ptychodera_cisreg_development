---
title: "Ptychodera Cisreg Development: Gene Functional Annotation"
author: "Alberto Perez-Posada"
date: "4/16/2023"
output:
  html_document: default
  pdf_document: default
editor_options: 
  markdown: 
    wrap: sentence
---


## About

The main GO annotation of the predicted set of protein sequences was done using BLAST2GO. For the annotation of trans-developmental genes and the functional categories, we used EggNOG against the EggNOG metazoan database(meNOG). After this we blasted and scanned against this same database and the pfam database for a more conservative threshold in sequence homology and domain architecture.

## GO annotation using EggNOG

First the EggNOG mapping. We kept anything with significance below 10E-5 and only one-to-one orthologs. All kinds of evidence sources.

```sh
OUTDIR="~/projects/ptychodera_cisreg_development/outputs/functional_annotation/eggnog/pfla_meNOG_one2one"

emapper.py -m diamond --database meNOG --target_orthologs one2one --hmm_evalue 1e-5 -i /home/ska/aperpos/Def_Pfla/outputs/transdecoder/ptyFlav3_CYi_longest.cds.fa.transdecoder.pep --output $OUTDIR

```

### Blast against the meNOG database

At minimum anything below 10E-3, with a top of 3 similar target sequences and a minimum word size of 3.

```sh
blastp -db meNOG_all.fa -query ~/Def_Pfla/genome/ptyFlav3_CYi_longest.pep.faa -num_threads 12 -evalue 10e-3 -max_target_seqs 3 -word_size 3 -outfmt "6 qseqid sseqid qlen slen length gaps bitscore evalue" > pfla_vs_meNOG.tsv
```

### Pfam search

This search was done using standard parameters.

```sh
hmmscan --cpu 8 --noali -E 10e-3 --domE 10e-3 -o ptyFlav3_CYi_longest.cds.fa.transdecoder.pep.pfam.outstd --tblout ptyFlav3_CYi_longest.cds.fa.transdecoder.pep.pfam.tblout --domtblout ptyFlav3_CYi_longest.cds.fa.transdecoder.pep.pfam.domtblout --pfamtblout ptyFlav3_CYi_longest.cds.fa.transdecoder.pep.pfam.pfamtblou Pfam-A.hmm

perl -pe "s[ ]{2,}/\t/g" ptyFlav3_CYi_longest.cds.fa.transdecoder.pep.pfam.tblout | perl -pe "s/\ TCONS/\tTCONS/" > ptyFlav3_CYi_longest.cds.fa.transdecoder.pep.pfam.tblout.tsv
```

### Curate with Blast and pfam

The main points are to keep:

 - hits with query coverage higher than 50%,
 - hits with bitscore above 70
 - Add any GO hit based on Pfam domains that was not caught by EggNOG

```r
#EGGNOG CURATION W BLAST

pfmenog <- read.table("/home/ska/aperpos/Def_Pfla/outputs/eggnog/menog_1to1/pfla_vs_meNOG.tsv",header=F,sep="\t",dec=".")
colnames(pfmenog) <- c("qseqid", "sseqid", "qlen", "slen", "length", "gaps", "bitscore", "evalue")

pfmenog$qcov <- (pfmenog$length-pfmenog$gaps)/pfmenog$qlen

pfmenog <- pfmenog[pfmenog$qcov > 0.5,]
pfmenog <- pfmenog[pfmenog$bitscore > 70,]

pfla_meNOG_besthits <- unique(pfmenog[,1:2])
pfla_meNOG_besthits$queryseed <- paste(pfla_meNOG_besthits$qseqid,pfla_meNOG_besthits$sseqid,sep="_")


pfla_menog_1to1 <- read.table("/home/ska/aperpos/Def_Pfla/outputs/eggnog/menog_1to1/ptyFlav3_CYi_longest.transdecoder_clean.faa_meNOG_1to1.emapper.annotations",sep="\t",header=T,comment.char="#", stringsAsFactors=FALSE, quote="", fill=FALSE)
pfla_menog_1to1$queryseed <- paste(pfla_menog_1to1$query_name,pfla_menog_1to1$seed_eggNOG_ortholog,sep="_")
pfla_menog_1to1 <- pfla_menog_1to1[pfla_menog_1to1$queryseed %in% pfla_meNOG_besthits$queryseed,]
rownames(pfla_menog_1to1) <- NULL

#EGGNOG CURATION W PFAM

pfla_pfam <- read.delim("~/Def_Pfla/outputs/pfam/ptyFlav3_CYi_longest.transdecoder_clean.faa.pep.pfam.tblout.tsv", header=FALSE, comment.char="#") [,c(2,3)]
pfam2go <- read.delim("~/Def_Pfla/outputs/pfam/pfam2go.tsv", header=FALSE, comment.char="!")[,-c(2,3)]
pfam2go$V1 <- sub("Pfam:","",pfam2go$V1)
pfla_pfam$V2 <- sub("\\..*","",pfla_pfam$V2)
pfla_pfam <- pfla_pfam[pfla_pfam$V2!="",]
pfla_pfam2go <- merge(pfam2go,pfla_pfam,by=1)

for (i in 1:nrow(pfla_pfam2go)){
  tcons <- pfla_pfam2go[i,3]
  go <- pfla_pfam2go[i,2]
  i_golist <- pfla_menog_1to1[grep(tcons,pfla_menog_1to1$query_name),6]
  if(go %in% i_golist == F) i_golist <- c(i_golist,go)
  pfla_menog_1to1[grep(tcons,pfla_menog_1to1$query_name),6] <- paste(i_golist,sep="",collapse=",")
}
pfla_menog_1to1$GO_terms <- sub(",","",pfla_menog_1to1$GO_terms)

write.table(pfla_menog_1to1[,-14],"/home/ska/aperpos/Def_Pfla/outputs/eggnog/menog_1to1/curated_qcov_hmm/ptyFlav3_CYi_longest.transdecoder_clean.faa_meNOG_1to1.emapper.annotations.qcov.hmm.curated.tsv",sep="\t",quote=F,row.names=F)
```

## Trans-developmental and House-keeping genes

We defined "trans-developmental" genes as any genes containing one or more of the following GO terms in their GO annotation, and NONE of the "house-keeping" GO terms.

```
GO:0009790 embryo development
GO:0030154 cell differentiation
GO:0043565 DNA binding
GO:0007267 cell-cell signaling 
GO:0008380 RNA splicing 
GO:0003700 DNA-binding transcription factor activity  
GO:0140297 DNA-binding transcription factor binding
GO:0010467 gene expression
GO:0010468 regulation of gene expression 
GO:0016070 RNA metabolic process
```

We defined "house-keeping" genes as any genes containing one or more of the following GO terms in their GO annotation, and NONE of the "trans-developmental" GO terms.

```
GO:0044430 cytoskeletal part
GO:0006412 protein biosynthesis
GO:0003735 structural constituent of ribosome
GO:0005840 ribosome
GO:0005842 cytosolic large ribos. subunit (sensu Euk.)
GO:0030529 ribonucleoprotein complex
GO:0006414 translational elongation
GO:0009119 ribonucleoside metabolic process
GO:0044237 cellular metabolic process
GO:0009062 fatty acid catabolic process
GO:0044262 cellular carbohydrate metabolic process
GO:0006725 cellular aromatic compound metabolic process
```

And the command to fish out trans-dev and housekeeping genes:

```sh
GO_FILE="~/projects/ptychodera_cisreg_development/outputs/functional_annotation/eggnog/menog_1to1_qcov_hmm_curated.goterms"
OUTDIR="~/projects/ptychodera_cisreg_development/outputs/functional_annotation/TD_annotation"

TD_GOS="GO:0009790\|GO:0030154\|GO:0043565\|GO:0007267\|GO:0008380\|GO:0003700\|GO:0140297\|GO:0010467\|GO:0010468\|GO:0016070"
HK_GOS="GO:0044430\|GO:0006412\|GO:0003735\|GO:0005840\|GO:0005842\|GO:0030529\|GO:0006414\|GO:0009119\|GO:0009062\|GO:0044262\|GO:0006725"

grep "${TD_GOS}" $GO_FILE | grep -v "${HK_GOS}" | cut -f1 | sort | uniq > ${OUTDIR}/transdev_expanded_nohk.tsv
grep "${HK_GOS}" $GO_FILE | grep -v "${TD_GOS}" | cut -f1 | sort | uniq > ${OUTDIR}/housekeep_expanded_notd.tsv
```