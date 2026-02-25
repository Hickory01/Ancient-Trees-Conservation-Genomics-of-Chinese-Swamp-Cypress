python 03_ac_wac_irr.py \
  --allele-table allele_table.with_flags.tsv \
  --gt-tsv all147.gt.tsv \
  --samples-order samples.order.txt \
  --ancients ancients64.list \
  --anc-nat anc_nat26.list \
  --anc-cult anc_cult38.list \
  --anc-admix anc_admix4.list \
  --anc-min anc_Min43.list \
  --anc-zhu anc_Zhu17.list \
  --cultivated cult50.list \
  --wild wild33.list \
  --epsilon 1e-3 \
  --max-occ 2 \
  --irr-coverage cult \
  --out-prefix gpen_acwac_full

