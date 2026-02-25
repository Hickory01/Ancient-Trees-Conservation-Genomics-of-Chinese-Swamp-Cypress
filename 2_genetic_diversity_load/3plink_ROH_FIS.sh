vcf00="/home/ZhangWP/water_pine/snp_vcf/Gpen147.DBN20.recode.vcf"

#ROH
plink --vcf $vcf00 -allow-extra-chr --make-bed --out Gpen147_DBN20

plink \
        --bfile Gpen147_DBN20 --allow-extra-chr \
        --homozyg \
        --homozyg-density 50 \
        --homozyg-gap 500 \
        --homozyg-kb 100 \
        --homozyg-snp 50 \
        --homozyg-window-het 1 \
        --homozyg-window-snp 50 \
        --homozyg-window-threshold 0.05 \
        --homozyg-window-missing 5 \
        --homozyg-het 1 \
        --out Gpen147_DBN20_ROH

#FIS
vcftools --vcf $vcf00 --het --out Gpen147_DBN20_het
