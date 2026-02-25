####run in 171.24.64.50####
##use perl script to produce the input file of est-sfs
species_est_sfs="/panfs3/home/zhangwp/water/4_mutation_road/new_20250611/new_20251123/Gpen147_out2_DBN/sp_ind_est_sfs1"
config_jc="/panfs3/home/zhangwp/software/est-sfs-release-2.04/config-JC.txt"
config_kimura="/panfs3/home/zhangwp/software/est-sfs-release-2.04/config-kimura.txt"
config_rate6="/panfs3/home/zhangwp/software/est-sfs-release-2.04/config-rate6.txt"
seed="/panfs3/home/zhangwp/software/est-sfs-release-2.04/seedfile.txt"

for chr in `cat $1`
do
perl est_sfs_input.dfe.Gpen147.pl $species_est_sfs/Gpen147_$chr.whole.est_sfs.frq $species_est_sfs/out1_$chr.whole.est_sfs.frq $species_est_sfs/out2_$chr.whole.est_sfs.frq
nohup est-sfs $config_rate6 $species_est_sfs/Gpen147_$chr.whole.est_sfs.input.txt $seed $species_est_sfs/Gpen147_$chr.whole.sfs.output.txt $species_est_sfs/Gpen147_$chr.whole.est_sfs.output.p_anc.txt &
done
