#convert
plink --vcf 1m_Gpen147_chr_jc.vcf --make-bed --out 1m_Gpen147

#kinship
king -b 1m_Gpen147.bed --kinship --prefix new1m_Gpen147_kinship

#unrelated
king -b 1m_Gpen147.bed --unrelated --degree 3 --prefix new1m_Gpen147_unrelated3
