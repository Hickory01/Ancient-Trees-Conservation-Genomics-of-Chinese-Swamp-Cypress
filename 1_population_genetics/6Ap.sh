vcffile1="/home/ZhangWP/water_pine/snp_vcf/Gpen147.DBN20.recode.vcf"

nohup vcftools --vcf $vcffile1 --singletons --out 1_Gpen147_DBN20_Ap >1_Gpen147_DBN20_Ap_log 2>&1 &
