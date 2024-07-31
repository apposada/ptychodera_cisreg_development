# Ptychodera Cisreg Development : Retrieving orthologues for kernel subgraphs

## Preparing stuff

### Human/P.flava reciprocal best hits 
We will run a whole reciprocal best hits (rbh) BLAST of the Human proteome we used against the P. flava proteome. This is to use in addition to OMA and Orthofinder outputs to retrieve putative orthologues.

```sh
# pfla,HUMAN rbh
cd ~/projects/ptychodera_cisreg_development/outputs/comparative/rbh/
ln -s ~/projects/ptychodera_cisreg_development/outputs/comparative/20240404_orthofinder/proteomes/HUMAN.fa HUMAN.faa
makeblastdb -dbtype prot -in HUMAN.faa
# pfla against human
blastp -db HUMAN.faa -query pfla.faa -num_threads 12 -evalue 1e-3 -max_target_seqs 1 -outfmt "6 qseqid sseqid qlen slen length gaps qcovs qcovhsp bitscore evalue" | tee pfla_HUMAN.out | cut -f1,2 | sort | uniq > pfla_HUMAN.tsv
# human against pfla
blastp -db pfla.faa -query HUMAN.faa -num_threads 12 -evalue 1e-3 -max_target_seqs 1 -outfmt "6 qseqid sseqid qlen slen length gaps qcovs qcovhsp bitscore evalue" | tee HUMAN_pfla.out | awk '{OFS = "\t"} {print $2,$1}' | sort | uniq > HUMAN_pfla.tsv
# script to keep the common lines
python common_lines.py HUMAN_pfla.tsv pfla_HUMAN.tsv pfla_HUMAN_rbh.tsv

# link to the working directory
cd ~/projects/ptychodera_cisreg_development/outputs/functional_annotation/grn_homologs/orthologues/
ln -s ~/projects/ptychodera_cisreg_development/outputs/comparative/rbh/pfla_HUMAN_rbh.tsv .
```
Then we link the pairwise orthologs from OMA and Orthofinder, for pfla-spur and pfla-hsap:

```sh
ln -s ~/projects/ptychodera_cisreg_development/outputs/comparative/20240404_oma/precompiled/Output/PairwiseOrthologs/STRPU_ok-PFLAV.txt oma_spur_pfla_all_pairwise_orthologues.tsv
ln -s ~/projects/ptychodera_cisreg_development/outputs/comparative/20240404_oma/precompiled/Output/PairwiseOrthologs/HUMAN-PFLAV.txt oma_human_pfla_all_pairwise_orthologues.tsv
ln -s ~/projects/ptychodera_cisreg_development/outputs/comparative/20240404_orthofinder/proteomes/OrthoFinder/Results_Apr05/Orthologues/Orthologues_HUMAN/HUMAN__v__PFLAV.tsv of_human_pfla_all_orthologues.tsv
ln -s ~/projects/ptychodera_cisreg_development/outputs/comparative/rbh/pfla_HUMAN_rbh.tsv .
```

## The Endomesoderm Kernel

We retrieved the echinobase/RefSeq Spur5.0 sequences of a list of endomesoderm kernel genes taken from Davidson, 2005.

As we do not have the exactly same Spur proteome (ours is a pre-compiled one from OMA), to be entirely sure that we are using the right sequences, we will blast the ones we retrieved to our proteome of Spur.

```sh
# link the proteome, make a db
ln -s ~/projects/ptychodera_cisreg_development/outputs/comparative/20240404_orthofinder/proteomes/STRPU_ok.fa .
makeblastdb -dbtype prot -in STRPU_ok.fa

# blast the sequences I got against the sequences in the proteomes I used in Orthofinder / OMA
blastp -db STRPU_ok.fa -query endomesod_spur_echinobase_refseq.faa -num_threads 12 -max_target_seqs 3 -evalue 10e-10 -outfmt 6 -out endomesod_strpu.out 
```

After this, we grep them out of the OMA Pfla/Spur pairwise orthologs output file from OMA, and we store them separately.

```sh
# keep these genes and these only
cat endomesod_strpu.out | awk '$11=="0.0"' | cut -f2 | grep -f - oma_spur_pfla_all_pairwise_orthologues.tsv | awk '{OFS="\t"} {print($3,$4,$5)}' | perl -pe "s/\.p[0-9]+//" > endomesod_pfla_spur.tsv
```
We manually tidied up this table to add the gene symbol column.

## The axial mesoderm kernel

We retrieve the Uniprot Human Sequences of a list of notochord/axial mesoderm GRN gene lists from Zhang et al, 2020 (a single cell atlas of C. intestinalis) and Di Gregorio et al, 2021 (a review on the notochord GRN).

Just like before, we BLAST them to the Human proteome we used in our analyses to keep the exact matches.

```sh
# link the proteome, make a db
ln -s ~/projects/ptychodera_cisreg_development/outputs/comparative/20240404_orthofinder/proteomes/HUMAN.fa .
makeblastdb -dbtype prot -in HUMAN.fa
# blast the sequences I got against the sequences in the proteomes I used in Orthofinder / OMA
blastp -db HUMAN.fa -query axmesod_hsap_zhang_digregorio_uniprot.faa -num_threads 12 -max_target_seqs 3 -evalue 10e-10 -outfmt 6 -out axmesod_human.out 
```

This is a bit more involved. We used hits from OMA, Orthofinder, and Reciprocal Best hits. These outputs were manually collected and tidied up. We kept all the hits with ambiguous assignation (e.g. we got two many:many orthologues of Snail in Orthofinder but just one in OMA, and tagged them as such). This did not happen with many genes, the other case being Plod1/2/3 which is again a many:1 case, or FoxD1.

```sh
cat axmesod_human.out | awk '$11=="0.0"' | cut -f2 | grep -f - oma_human_pfla_all_pairwise_orthologues.tsv | awk '{OFS="\t"} {print($4,$3)}' | perl -pe "s/\.p[0-9]+//" > axmesod_pfla_human.tsv

cat axmesod_human.out | awk '$11=="0.0" || $3>90' | cut -f2 | grep -f - pfla_HUMAN_rbh.tsv
cat axmesod_human.out | awk '$11=="0.0" || $3>90' | cut -f2 | grep -f - of_human_pfla_all_orthologues.tsv
cat axmesod_human.out | awk '$11=="0.0" || $3>90' | cut -f2 | grep -f - oma_human_pfla_all_pairwise_orthologues.tsv | awk '{OFS="\t"} {print($3,$4,$5)}' | perl -pe "s/\.p[0-9]+//" > axmesod_pfla_human.tsv
```

In addition we added Nodal, Delta, Brachyury, and FoxA as we have information from these from in situ hybridisations.
