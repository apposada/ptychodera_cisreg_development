while read LINE; do
stage=`echo $LINE | cut -d ' ' -f1`
sample_query=`echo $LINE | cut -d ' ' -f2`
sample_dir=`ls /home/ska/aperpos/Def_Pfla/other_spp/data_spur/spur_Tu/RNA/kallisto_20210330 | grep $sample_query`
echo "cp /home/ska/aperpos/Def_Pfla/other_spp/data_spur/spur_Tu/RNA/kallisto_20210330/$sample_dir/abundance.h5 $stage/ "
cp /home/ska/aperpos/Def_Pfla/other_spp/data_spur/spur_Tu/RNA/kallisto_20210330/$sample_dir/abundance.h5 $stage/
done < stage_to_sample.tsv
