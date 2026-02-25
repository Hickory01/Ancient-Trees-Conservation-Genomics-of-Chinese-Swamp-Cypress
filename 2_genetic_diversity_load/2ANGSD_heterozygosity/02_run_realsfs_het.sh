for n in `cat ind147`
do
sh run_realsfs_het.sh ./listdir/list${n} chr11.list /home/ZhangWP/water_pine/Gpen_ref_index/Gpen_refgenome_100k.fa 1 &>./logout/log_het_list${n}
done
