for id in `cat $1`
do
pixy --stats pi \
--vcf Gpen147.mis20.filtered1.${id}.vcf.gz \
--sites_file Gpen147.mis20.filtered1.${id}.sites.txt \
--populations Gpen147_pop0.txt \
--window_size 200000 \
--output_folder ./pixy_Gpen147_total \
--output_prefix ${id} \
--n_cores 4
done1
