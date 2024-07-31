# Pfla
cat STRPU_ok-PFLAV.txt | cut -f3,4,5 | grep "1:1" | awk 'BEGIN {OFS="\t"}{print($2,$1)}' | perl -pe "s/\.p[0-9]+//" > 1to1_pfla_spur.tsv
cat AENJA_ok-PFLAV.txt | cut -f3,4,5 | grep "1:1" | awk 'BEGIN {OFS="\t"}{print($2,$1)}' | perl -pe "s/\.p[0-9]+//" > 1to1_pfla_ajap.tsv
cat PFLAV-BRALA.txt | cut -f3,4,5 | grep "1:1" | awk 'BEGIN {OFS="\t"}{print($1,$2)}' | perl -pe "s/\.p[0-9]+//" > 1to1_pfla_blan.tsv
cat BRAFL-PFLAV.txt | cut -f3,4,5 | grep "1:1" | awk 'BEGIN {OFS="\t"}{print($2,$1)}' | perl -pe "s/\.p[0-9]+//" > 1to1_pfla_bflo.tsv
cat OWEFU-PFLAV.txt | cut -f3,4,5 | grep "1:1" | awk 'BEGIN {OFS="\t"}{print($2,$1)}' | perl -pe "s/\.p[0-9]+//" > 1to1_pfla_ofus.tsv

# Spur
cat AENJA_ok-STRPU_ok.txt | cut -f3,4,5 | grep "1:1" | awk 'BEGIN {OFS="\t"}{print($2,$1)}' | perl -pe "s/\.p[0-9]+//" > 1to1_spur_ajap.tsv
cat STRPU_ok-BRALA.txt | cut -f3,4,5 | grep "1:1" | awk 'BEGIN {OFS="\t"}{print($1,$2)}' | perl -pe "s/\.p[0-9]+//" > 1to1_spur_ajap.tsv
cat STRPU_ok-BRAFL.txt | cut -f3,4,5 | grep "1:1" | awk 'BEGIN {OFS="\t"}{print($1,$2)}' | perl -pe "s/\.p[0-9]+//" > 1to1_spur_bflo.tsv
cat STRPU_ok-OWEFU.txt | cut -f3,4,5 | grep "1:1" | awk 'BEGIN {OFS="\t"}{print($1,$2)}' | perl -pe "s/\.p[0-9]+//" > 1to1_spur_ofus.tsv

# Ajap
cat AENJA_ok-BRALA.txt | cut -f3,4,5 | grep "1:1" | awk 'BEGIN {OFS="\t"}{print($1,$2)}' | perl -pe "s/\.p[0-9]+//" > 1to1_ajap_blan.tsv
cat AENJA_ok-BRAFL.txt | cut -f3,4,5 | grep "1:1" | awk 'BEGIN {OFS="\t"}{print($1,$2)}' | perl -pe "s/\.p[0-9]+//" > 1to1_ajap_bflo.tsv
cat AENJA_ok-OWEFU.txt | cut -f3,4,5 | grep "1:1" | awk 'BEGIN {OFS="\t"}{print($1,$2)}' | perl -pe "s/\.p[0-9]+//" > 1to1_ajap_ofus.tsv

# Blan
cat BRAFL-BRALA.txt | cut -f3,4,5 | grep "1:1" | awk 'BEGIN {OFS="\t"}{print($2,$1)}' | perl -pe "s/\.p[0-9]+//" > 1to1_blan_bflo.tsv
cat OWEFU-BRALA.txt | cut -f3,4,5 | grep "1:1" | awk 'BEGIN {OFS="\t"}{print($2,$1)}' | perl -pe "s/\.p[0-9]+//" > 1to1_blan_ofus.tsv

# Bflo
cat BRAFL-OWEFU.txt | cut -f3,4,5 | grep "1:1" | awk 'BEGIN {OFS="\t"}{print($1,$2)}' | perl -pe "s/\.p[0-9]+//" > 1to1_bflo_ofus.tsv
