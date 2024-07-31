## Pfla
cat ../proteomes/OrthoFinder/Results_Apr05/Orthologues/Orthologues_PFLAV/PFLAV__v__STRPU_ok.tsv | grep -v "," | perl -pe "s/\.p[0-9]+//" | cut -f2,3 | tail -n +2 > 1to1_pfla_spur.tsv

cat ../proteomes/OrthoFinder/Results_Apr05/Orthologues/Orthologues_PFLAV/PFLAV__v__AENJA_ok.tsv | grep -v "," | perl -pe "s/\.p[0-9]+//" | cut -f2,3 | tail -n +2 > 1to1_pfla_ajap.tsv

cat ../proteomes/OrthoFinder/Results_Apr05/Orthologues/Orthologues_PFLAV/PFLAV__v__BRALA.tsv | grep -v "," | perl -pe "s/\.p[0-9]+//" | cut -f2,3 | tail -n +2 > 1to1_pfla_blan.tsv

cat ../proteomes/OrthoFinder/Results_Apr05/Orthologues/Orthologues_PFLAV/PFLAV__v__BRAFL.tsv | grep -v "," | perl -pe "s/\.p[0-9]+//" | cut -f2,3 | tail -n +2 > 1to1_pfla_bflo.tsv

cat ../proteomes/OrthoFinder/Results_Apr05/Orthologues/Orthologues_PFLAV/PFLAV__v__OWEFU.tsv | grep -v "," | perl -pe "s/\.p[0-9]+//" | cut -f2,3 | tail -n +2 > 1to1_pfla_ofus.tsv

## Spur
cat ../proteomes/OrthoFinder/Results_Apr05/Orthologues/Orthologues_STRPU_ok/STRPU_ok__v__AENJA_ok.tsv | grep -v "," | perl -pe "s/\.p[0-9]+//" | cut -f2,3 | tail -n +2 > 1to1_spur_ajap.tsv

cat ../proteomes/OrthoFinder/Results_Apr05/Orthologues/Orthologues_STRPU_ok/STRPU_ok__v__BRALA.tsv | grep -v "," | perl -pe "s/\.p[0-9]+//" | cut -f2,3 | tail -n +2 > 1to1_spur_blan.tsv

cat ../proteomes/OrthoFinder/Results_Apr05/Orthologues/Orthologues_STRPU_ok/STRPU_ok__v__BRAFL.tsv | grep -v "," | perl -pe "s/\.p[0-9]+//" | cut -f2,3 | tail -n +2 > 1to1_spur_bflo.tsv

cat ../proteomes/OrthoFinder/Results_Apr05/Orthologues/Orthologues_STRPU_ok/STRPU_ok__v__OWEFU.tsv | grep -v "," | perl -pe "s/\.p[0-9]+//" | cut -f2,3 | tail -n +2 > 1to1_spur_ofus.tsv

## Ajap
cat ../proteomes/OrthoFinder/Results_Apr05/Orthologues/Orthologues_AENJA_ok/AENJA_ok__v__BRALA.tsv | grep -v "," | perl -pe "s/\.p[0-9]+//" | cut -f2,3 | tail -n +2 > 1to1_ajap_blan.tsv

cat ../proteomes/OrthoFinder/Results_Apr05/Orthologues/Orthologues_AENJA_ok/AENJA_ok__v__BRAFL.tsv | grep -v "," | perl -pe "s/\.p[0-9]+//" | cut -f2,3 | tail -n +2 > 1to1_ajap_bflo.tsv

cat ../proteomes/OrthoFinder/Results_Apr05/Orthologues/Orthologues_AENJA_ok/AENJA_ok__v__OWEFU.tsv | grep -v "," | perl -pe "s/\.p[0-9]+//" | cut -f2,3 | tail -n +2 > 1to1_ajap_ofus.tsv

## Blan
cat ../proteomes/OrthoFinder/Results_Apr05/Orthologues/Orthologues_BRALA/BRALA__v__BRAFL.tsv | grep -v "," | perl -pe "s/\.p[0-9]+//" | cut -f2,3 | tail -n +2 > 1to1_blan_bflo.tsv

cat ../proteomes/OrthoFinder/Results_Apr05/Orthologues/Orthologues_BRALA/BRALA__v__OWEFU.tsv | grep -v "," | perl -pe "s/\.p[0-9]+//" | cut -f2,3 | tail -n +2 > 1to1_blan_ofus.tsv

## Bflo
cat ../proteomes/OrthoFinder/Results_Apr05/Orthologues/Orthologues_BRAFL/BRAFL__v__OWEFU.tsv | grep -v "," | perl -pe "s/\.p[0-9]+//" | cut -f2,3 | tail -n +2 > 1to1_bflo_ofus.tsv
