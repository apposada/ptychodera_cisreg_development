# pfla,spur rbh
ln -s /home/ska/aperpos/Def_Pfla/outputs/transdecoder/using_ptyFlav3_CYi_longest.mrna/ptyFlav3_CYi_longest.transdecoder_clean.faa pfla.faa
makeblastdb -dbtype prot -in pfla.faa 

ln -s ~/projects/ptychodera_cisreg_development/outputs/comparative/oma/OMA_precompiled/DB/STRPU_ok.fa spur.faa
makeblastdb -dbtype prot -in spur.faa 

blastp -db spur.faa -query pfla.faa -num_threads 12 -evalue 1e-3 -max_target_seqs 1 -outfmt "6 qseqid sseqid qlen slen length gaps qcovs qcovhsp bitscore evalue" | tee pfla_spur.out | cut -f1,2 | sort | uniq > pfla_spur.tsv
blastp -db pfla.faa -query spur.faa -num_threads 12 -evalue 1e-3 -max_target_seqs 1 -outfmt "6 qseqid sseqid qlen slen length gaps qcovs qcovhsp bitscore evalue" | tee spur_pfla.out | awk '{OFS = "\t"} {print $2,$1}' | sort | uniq > spur_pfla.tsv

python common_lines.py spur_pfla.tsv pfla_spur.tsv pfla_spur_rbh.tsv

# pfla,blan rbh
ln -s /home/ska/aperpos/Def_Pfla/outputs/oma/depr_previous_versions/20210601_OMA.2.4.1_Spur5.0_ptyFlav3_singlebestonly/DB/BRALA.fa blan.faa
makeblastdb -dbtype prot -in blan.faa

blastp -db blan.faa -query pfla.faa -num_threads 12 -evalue 1e-3 -max_target_seqs 1 -outfmt "6 qseqid sseqid qlen slen length gaps qcovs qcovhsp bitscore evalue" | tee pfla_blan.out | cut -f1,2 | sort | uniq > pfla_blan.tsv
blastp -db pfla.faa -query blan.faa -num_threads 12 -evalue 1e-3 -max_target_seqs 1 -outfmt "6 qseqid sseqid qlen slen length gaps qcovs qcovhsp bitscore evalue" | tee blan_pfla.out | awk '{OFS = "\t"} {print $2,$1}' | sort | uniq > blan_pfla.tsv

python common_lines.py blan_pfla.tsv pfla_blan.tsv pfla_blan_rbh.tsv

