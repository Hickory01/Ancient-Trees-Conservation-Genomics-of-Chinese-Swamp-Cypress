## 样本顺序
#bcftools query -l /home/ZhangWP/water_pine/snp_vcf/Gpen147.DBN20.recode.vcf > samples.order.txt

# 计算（按需加上/去掉可选子群）
python 02_allele_count.py \
  --gt-tsv all147.gt.tsv \
  --samples-order samples.order.txt \
  --ancients ancients64.list \
  --cultivated cult50.list \
  --wild wild33.list \
  --anc-nat anc_nat26.list \
  --anc-cult anc_cult38.list \
  --anc-admix anc_admix4.list \
  --anc-min anc_Min43.list \
  --anc-zhu anc_Zhu17.list \
  --out allele_table.with_flags.tsv \
  --chunksize 200000

