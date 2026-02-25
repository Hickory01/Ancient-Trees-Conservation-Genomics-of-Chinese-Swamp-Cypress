vcf00="/home/ZhangWP/water_pine/snp_vcf/Gpen66.DBN20.CDS200.MAF01.LD200.jc.vcf"
out00="/home/ZhangWP/water_pine/snp_vcf/1_pca/4th_pca"

vcftools --vcf $vcf00 --plink --out $out00/Gpen66_mis20
plink --file $out00/Gpen66_mis20 --make-bed --out $out00/Gpen66_mis20
plink --file $out00/Gpen66_mis20 --pca 3 --out $out00/Gpen66_mis20 --noweb
