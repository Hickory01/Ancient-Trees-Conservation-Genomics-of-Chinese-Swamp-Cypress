source /opt/software_by_fc/gcc/5.5.0.env

indvcffile="/home/ZhangWP/water_pine/mutation_road/new_4Jun2025/Gpen147_out2_20_missing/sp_ind_vcf"
species_anno_vcf="/home/ZhangWP/water_pine/mutation_road/new_4Jun2025/Gpen147_out2_20_missing/new20251122/sp_ind_hom_het_vcf"
ancestor1_sites="/home/ZhangWP/water_pine/mutation_road/new_4Jun2025/Gpen147_out2_20_missing/new20251122/allele_REF_is_ancestor.txt"
ancestor2_sites="/home/ZhangWP/water_pine/mutation_road/new_4Jun2025/Gpen147_out2_20_missing/new20251122/allele_ALT_is_ancestor.txt"

for id in `cat $1`
do

####step4_1 extract cds & ancestor1 & hom_het sites
vcftools --vcf $indvcffile/$id.anno.vcf --positions $ancestor1_sites --recode --recode-INFO-all --out $species_anno_vcf/$id.anno.out2.homo2.cds.ancestor1

grep $'\t0/1:' $species_anno_vcf/$id.anno.out2.homo2.cds.ancestor1.recode.vcf >$species_anno_vcf/$id.anno.out2.homo2.cds.ancestor1.het.vcf
grep $'\t0|1:' $species_anno_vcf/$id.anno.out2.homo2.cds.ancestor1.recode.vcf >>$species_anno_vcf/$id.anno.out2.homo2.cds.ancestor1.het.vcf

sort $species_anno_vcf/$id.anno.out2.homo2.cds.ancestor1.het.vcf | uniq >$species_anno_vcf/$id.anno.out2.homo2.cds.ancestor1.het1.vcf

grep $'\t1/1:' $species_anno_vcf/$id.anno.out2.homo2.cds.ancestor1.recode.vcf >$species_anno_vcf/$id.anno.out2.homo2.cds.ancestor1.hom.vcf
grep $'\t1|1:' $species_anno_vcf/$id.anno.out2.homo2.cds.ancestor1.recode.vcf >>$species_anno_vcf/$id.anno.out2.homo2.cds.ancestor1.hom.vcf

sort $species_anno_vcf/$id.anno.out2.homo2.cds.ancestor1.hom.vcf | uniq >$species_anno_vcf/$id.anno.out2.homo2.cds.ancestor1.hom1.vcf

####step4_2 extract cds & ancestor2 & hom_het sites
vcftools --vcf $indvcffile/$id.anno.vcf --positions $ancestor2_sites --recode --recode-INFO-all --out $species_anno_vcf/$id.anno.out2.homo2.cds.ancestor2

grep $'\t0/1:' $species_anno_vcf/$id.anno.out2.homo2.cds.ancestor2.recode.vcf >$species_anno_vcf/$id.anno.out2.homo2.cds.ancestor2.het.vcf
grep $'\t0|1:' $species_anno_vcf/$id.anno.out2.homo2.cds.ancestor2.recode.vcf >>$species_anno_vcf/$id.anno.out2.homo2.cds.ancestor2.het.vcf

sort $species_anno_vcf/$id.anno.out2.homo2.cds.ancestor2.het.vcf | uniq >$species_anno_vcf/$id.anno.out2.homo2.cds.ancestor2.het1.vcf

grep $'\t0/0:' $species_anno_vcf/$id.anno.out2.homo2.cds.ancestor2.recode.vcf >$species_anno_vcf/$id.anno.out2.homo2.cds.ancestor2.hom.vcf
grep $'\t0|0:' $species_anno_vcf/$id.anno.out2.homo2.cds.ancestor2.recode.vcf >>$species_anno_vcf/$id.anno.out2.homo2.cds.ancestor2.hom.vcf

sort $species_anno_vcf/$id.anno.out2.homo2.cds.ancestor2.hom.vcf | uniq >$species_anno_vcf/$id.anno.out2.homo2.cds.ancestor2.hom1.vcf

####step5 extract synonymous(synonymous_variant), loss of function(start_lost, stop_gained, stop_lost), missense(missense_variant) 
grep "synonymous_variant" $species_anno_vcf/$id.anno.out2.homo2.cds.ancestor1.hom1.vcf >$species_anno_vcf/$id.anno.out2.homo2.cds.ancestor1.hom.syn.vcf
grep "synonymous_variant" $species_anno_vcf/$id.anno.out2.homo2.cds.ancestor1.het1.vcf >$species_anno_vcf/$id.anno.out2.homo2.cds.ancestor1.het.syn.vcf

grep "start_lost" $species_anno_vcf/$id.anno.out2.homo2.cds.ancestor1.hom1.vcf >$species_anno_vcf/$id.anno.out2.homo2.cds.ancestor1.hom.lof.vcf
grep "stop_gained" $species_anno_vcf/$id.anno.out2.homo2.cds.ancestor1.hom1.vcf >>$species_anno_vcf/$id.anno.out2.homo2.cds.ancestor1.hom.lof.vcf
grep "stop_lost" $species_anno_vcf/$id.anno.out2.homo2.cds.ancestor1.hom1.vcf >>$species_anno_vcf/$id.anno.out2.homo2.cds.ancestor1.hom.lof.vcf

grep "start_lost" $species_anno_vcf/$id.anno.out2.homo2.cds.ancestor1.het1.vcf >$species_anno_vcf/$id.anno.out2.homo2.cds.ancestor1.het.lof.vcf
grep "stop_gained" $species_anno_vcf/$id.anno.out2.homo2.cds.ancestor1.het1.vcf >>$species_anno_vcf/$id.anno.out2.homo2.cds.ancestor1.het.lof.vcf
grep "stop_lost" $species_anno_vcf/$id.anno.out2.homo2.cds.ancestor1.het1.vcf >>$species_anno_vcf/$id.anno.out2.homo2.cds.ancestor1.het.lof.vcf

grep "missense_variant" $species_anno_vcf/$id.anno.out2.homo2.cds.ancestor1.hom1.vcf >$species_anno_vcf/$id.anno.out2.homo2.cds.ancestor1.hom.mis.vcf
grep "missense_variant" $species_anno_vcf/$id.anno.out2.homo2.cds.ancestor1.het1.vcf >$species_anno_vcf/$id.anno.out2.homo2.cds.ancestor1.het.mis.vcf


grep "synonymous_variant" $species_anno_vcf/$id.anno.out2.homo2.cds.ancestor2.hom1.vcf >$species_anno_vcf/$id.anno.out2.homo2.cds.ancestor2.hom.syn.vcf
grep "synonymous_variant" $species_anno_vcf/$id.anno.out2.homo2.cds.ancestor2.het1.vcf >$species_anno_vcf/$id.anno.out2.homo2.cds.ancestor2.het.syn.vcf

grep "start_lost" $species_anno_vcf/$id.anno.out2.homo2.cds.ancestor2.hom1.vcf >$species_anno_vcf/$id.anno.out2.homo2.cds.ancestor2.hom.lof.vcf
grep "stop_gained" $species_anno_vcf/$id.anno.out2.homo2.cds.ancestor2.hom1.vcf >>$species_anno_vcf/$id.anno.out2.homo2.cds.ancestor2.hom.lof.vcf
grep "stop_lost" $species_anno_vcf/$id.anno.out2.homo2.cds.ancestor2.hom1.vcf >>$species_anno_vcf/$id.anno.out2.homo2.cds.ancestor2.hom.lof.vcf

grep "start_lost" $species_anno_vcf/$id.anno.out2.homo2.cds.ancestor2.het1.vcf >$species_anno_vcf/$id.anno.out2.homo2.cds.ancestor2.het.lof.vcf
grep "stop_gained" $species_anno_vcf/$id.anno.out2.homo2.cds.ancestor2.het1.vcf >>$species_anno_vcf/$id.anno.out2.homo2.cds.ancestor2.het.lof.vcf
grep "stop_lost" $species_anno_vcf/$id.anno.out2.homo2.cds.ancestor2.het1.vcf >>$species_anno_vcf/$id.anno.out2.homo2.cds.ancestor2.het.lof.vcf

grep "missense_variant" $species_anno_vcf/$id.anno.out2.homo2.cds.ancestor2.hom1.vcf >$species_anno_vcf/$id.anno.out2.homo2.cds.ancestor2.hom.mis.vcf
grep "missense_variant" $species_anno_vcf/$id.anno.out2.homo2.cds.ancestor2.het1.vcf >$species_anno_vcf/$id.anno.out2.homo2.cds.ancestor2.het.mis.vcf

done

