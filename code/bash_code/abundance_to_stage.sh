while read LINE; do
stage=`echo $LINE | cut -d ' ' -f1`
sample_query=`echo $LINE | cut -d ' ' -f2`
sample_dir=`ls /home/ska/aperpos/Def_Pfla/outputs/RNA_ptyFlav3/devstages/kallisto_20210303/ | grep $sample_query`
echo "cp /home/ska/aperpos/Def_Pfla/outputs/RNA_ptyFlav3/devstages/kallisto_20210303/$sample_dir/abundance.h5 $stage/ "
cp /home/ska/aperpos/Def_Pfla/outputs/RNA_ptyFlav3/devstages/kallisto_20210303/$sample_dir/abundance.h5 $stage/
done < stage_to_sample.tsv
