# 逐位点导出GT（双等位SNP过滤建议事先完成）
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' /home/ZhangWP/water_pine/snp_vcf/Gpen147.DBN20.recode.vcf > all147.gt.tsv
