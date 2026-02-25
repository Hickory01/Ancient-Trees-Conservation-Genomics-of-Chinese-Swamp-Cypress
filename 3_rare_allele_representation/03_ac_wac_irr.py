#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
03_ac_wac_irr.py  (fullset-147 + extra targets/covers)

功能：
1) 基于 allele_table.with_flags.tsv：
   - 计算 fullset 频率：AC_full = AC_anc + AC_cult + AC_wild；AN_full 同理；fa_full = AC_full/AN_full
   - 分箱：先抽出 singleton(AC_full==1)，其余按 fa_full：<0.5%、0.5–1%、1–5%、>5%
   - 目标集合（targets）：
       ancients64, anc_nat26(若有), anc_cult38(若有),
       anc_admix4(若有), anc_Min43(若有), anc_Zhu17(若有),
       wild33
   - 覆盖方式（covers）：
       covered_by_cultivated (= in_cult)
       covered_by_wild       (= in_wild)
       covered_by_cultivated+wild (= in_cult OR in_wild)
   - 输出所有你要求的 target×cover 组合（若某target缺少对应列则跳过）

2) 单棵古树 IRR_allele（默认以“栽培覆盖”为口径；可选“栽培或野生覆盖”）：
   - 对每个古树 i，取其“稀有古树等位”（anc_count ≤ max_occ），权重 w=-log10(max(fa_full, eps))
   - IRR(i)   = Σ w * (1 - cover_flag) ；cover_flag取 in_cult（默认）或 (in_cult|in_wild)
   - IRR_norm = IRR / Σ w（0–1）

输入（最关键的列由 02_allele_count.py 产生）：
  --allele-table allele_table.with_flags.tsv
  --gt-tsv       all147.gt.tsv
  --samples-order samples.order.txt
  --ancients ancients64.list
  --anc-nat anc_nat26.list        (可选)
  --anc-cult anc_cult38.list      (可选)
  --anc-admix anc_admix4.list     (可选)
  --anc-min   anc_Min43.list      (可选)
  --anc-zhu   anc_Zhu17.list      (可选)
  --cultivated cult50.list
  --wild wild33.list

可选：
  --epsilon 1e-3
  --max-occ 2
  --irr-coverage cult|cultwild   (默认cult)
  --chunksize 200000
  --out-prefix gpen_acwac_full
"""

import argparse
import numpy as np
import pandas as pd
import sys
import os

# ----------------- args -----------------
def parse_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--allele-table', required=True, help='allele_table.with_flags.tsv')
    p.add_argument('--gt-tsv',       required=True, help='all147.gt.tsv (%CHROM %POS %REF %ALT [GT×N])')
    p.add_argument('--samples-order', required=True, help='samples.order.txt')

    # groups / targets
    p.add_argument('--ancients', required=True, help='ancients64.list')
    p.add_argument('--anc-nat',  default=None, help='anc_nat26.list (optional)')
    p.add_argument('--anc-cult', default=None, help='anc_cult38.list (optional)')
    p.add_argument('--anc-admix', default=None, help='anc_admix4.list (optional)')
    p.add_argument('--anc-min',   default=None, help='anc_Min43.list (optional)')
    p.add_argument('--anc-zhu',   default=None, help='anc_Zhu17.list (optional)')
    p.add_argument('--cultivated', required=True, help='cult50.list')
    p.add_argument('--wild',       required=True, help='wild33.list')

    p.add_argument('--epsilon', type=float, default=1e-3, help='epsilon for -log10(freq) weights')
    p.add_argument('--max-occ', type=int, default=2, help='anc_count ≤ max_occ defines rare-in-ancients for IRR')
    p.add_argument('--irr-coverage', choices=['cult','cultwild'], default='cult',
                   help='IRR coverage flag: "cult" uses in_cult; "cultwild" uses (in_cult OR in_wild)')
    p.add_argument('--chunksize', type=int, default=200000, help='rows per chunk for reading GT TSV')
    p.add_argument('--out-prefix', default='gpen_acwac_full', help='output prefix')
    return p.parse_args()

def read_list(path):
    if path is None: return []
    return [x.strip() for x in open(path) if x.strip()]

# ----------------- helpers -----------------
def bin_label_full(ac_full, fa_full):
    if pd.notna(ac_full) and ac_full == 1:
        return 'singleton'
    if pd.isna(fa_full):
        return 'NA'
    if fa_full < 0.005:
        return '<0.5%'
    elif fa_full < 0.01:
        return '0.5–1%'
    elif fa_full < 0.05:
        return '1–5%'
    else:
        return '>5%'

def alt_count(gt: str):
    """ GT -> ALT=1 (0/1/2); any missing allele => NaN. Compatible with 0|1, 0/1, 1/1, 0/0, ./., 0/., ./0, single 0/1 """
    if gt is None:
        return np.nan
    gt = str(gt).strip()
    if gt == '' or gt == '.' or gt == './.' or gt == '.|.':
        return np.nan
    gt = gt.replace('|','/')
    parts = gt.split('/')
    if len(parts)==1:
        a = parts[0]
        if a == '.': return np.nan
        return 1.0 if a=='1' else 0.0
    a,b = parts[0], parts[1]
    if a=='.' or b=='.': return np.nan
    c = 0
    if a=='1': c+=1
    if b=='1': c+=1
    return float(c)

# ----------------- main -----------------
def main():
    args = parse_args()
    eps = args.epsilon

    # 读 allele_table
    df = pd.read_csv(args.allele_table, sep='\t',
                     dtype={'CHR':str,'POS':str,'REF':str,'ALT':str})

    # 必要列检查（来自 02_allele_count.py）
    base_needed = ['AC_anc','AN_anc','AC_cult','AN_cult','AC_wild','AN_wild',
                   'in_anc','in_cult','in_wild','anc_count']
    for c in base_needed:
        if c not in df.columns:
            raise SystemExit(f"[ERROR] missing column in allele_table: {c}")

    # fullset 频率与权重
    df['AC_full'] = df[['AC_anc','AC_cult','AC_wild']].sum(axis=1)
    df['AN_full'] = df[['AN_anc','AN_cult','AN_wild']].sum(axis=1)
    with np.errstate(divide='ignore', invalid='ignore'):
        df['fa_full'] = np.where(df['AN_full']>0, df['AC_full']/df['AN_full'], np.nan)
    df['bin'] = [bin_label_full(ac, f) for ac, f in zip(df['AC_full'], df['fa_full'])]
    df['w'] = -np.log10(np.clip(df['fa_full'].fillna(0.0).values, eps, None))

    # 组合出的覆盖标志
    df['in_cultwild'] = ((df['in_cult']==1) | (df['in_wild']==1)).astype(int)

    # 目标集合 S（按是否“在该集合出现过”）
    # 基本：
    S_defs = {'ancients64': (df['in_anc'] == 1)}
    # 可选子群（这些列在 02_allele_count.py 新增）
    if 'in_anc_nat' in df.columns:
        S_defs['anc_nat26']  = (df['in_anc_nat'] == 1)
    if 'in_anc_cult' in df.columns:
        S_defs['anc_cult38'] = (df['in_anc_cult'] == 1)
    if 'in_anc_admix' in df.columns:
        S_defs['anc_admix4'] = (df['in_anc_admix'] == 1)
    if 'in_anc_min' in df.columns:
        S_defs['anc_Min43']  = (df['in_anc_min'] == 1)
    if 'in_anc_zhu' in df.columns:
        S_defs['anc_Zhu17']  = (df['in_anc_zhu'] == 1)
    # wild33 作为 target：用 in_wild
    S_defs['wild33'] = (df['in_wild'] == 1)

    # 覆盖方定义（三种）
    cover_defs = {
        'covered_by_cultivated':   (df['in_cult'] == 1),
        'covered_by_wild':         (df['in_wild'] == 1),
        'covered_by_cultivated+wild': (df['in_cultwild'] == 1),
    }

    # ---------- AC / wAC：输出所有要求的组合 ----------
    out_rows = []
    for target, mask_S in S_defs.items():
        S = df[mask_S].copy()
        n_AS = len(S)
        if n_AS == 0:
            sys.stderr.write(f"[WARN] Target {target} has 0 alleles; skip\n")
            continue

        for cover_name, cover_mask in cover_defs.items():
            # 用户清单中只要求 wild33 与 covered_by_cultivated，其它 cover 可按需过滤
            if target == 'wild33' and cover_name != 'covered_by_cultivated':
                continue

            covered = cover_mask.loc[S.index].astype(bool)

            # overall
            AC_overall  = covered.mean()
            wAC_overall = float((S['w'] * covered).sum() / S['w'].sum()) if S['w'].sum()>0 else np.nan

            # by-bin（避免 groupby.apply 警告）
            tmp = S.assign(covered=covered.astype(float))
            tmp['w_cov'] = tmp['w'] * tmp['covered']
            bybin = tmp.groupby('bin', dropna=False).agg(
                n_alleles = ('covered','size'),
                AC        = ('covered','mean'),
                w_sum     = ('w','sum'),
                w_cov_sum = ('w_cov','sum')
            ).reset_index()
            bybin['wAC'] = bybin['w_cov_sum'] / bybin['w_sum']
            bybin = bybin[['bin','n_alleles','AC','wAC']]

            # 写每个组合的 by-bin 明细
            bybin_path = f"{args.out_prefix}.ac_wac_bybin.{target}.{cover_name}.csv"
            bybin.to_csv(bybin_path, index=False)

            # 汇总行（overall）
            out_rows.append({
                'target': target,
                'cover': cover_name,
                'n_AS': n_AS,
                'AC_overall': AC_overall,
                'wAC_overall': wAC_overall
            })
            # 同时把分箱作为“扩展行”附加到同一个 summary CSV（便于一次读取）
            for _, r in bybin.iterrows():
                out_rows.append({
                    'target': target,
                    'cover': cover_name,
                    'bin': r['bin'],
                    'n_in_bin': int(r['n_alleles']),
                    'AC_bin': float(r['AC']),
                    'wAC_bin': float(r['wAC']),
                })

    summary_df = pd.DataFrame(out_rows)
    summary_df.to_csv(f"{args.out_prefix}.ac_wac_summary.csv", index=False)
    sys.stderr.write(f"[OK] AC/wAC done -> {args.out_prefix}.ac_wac_summary.csv + per-bin CSVs for all targets/covers\n")

    # ---------- IRR_allele（逐古树；coverage 口径可切换） ----------
    # IRR 只针对古树（ancients64.list）逐个体计算；wild33 不参与 IRR。
    irr_cov_flag = 'in_cult' if args.irr_coverage == 'cult' else 'in_cultwild'

    # 准备 meta：key -> anc_count, fa_full, cover_flag, w
    meta_df = df[['CHR','POS','REF','ALT','anc_count',irr_cov_flag,'w']].copy()
    for c in ['CHR','POS','REF','ALT']:
        meta_df[c] = meta_df[c].astype(str)
    meta_df['key'] = meta_df['CHR']+'|'+meta_df['POS']+'|'+meta_df['REF']+'|'+meta_df['ALT']
    meta_df = meta_df.set_index('key')

    ancient_ids = read_list(args.ancients)
    anc_nat_ids  = set(read_list(args.anc_nat))  if args.anc_nat  else set()
    anc_cult_ids = set(read_list(args.anc_cult)) if args.anc_cult else set()

    # 只解析古树列
    samples_order = [s.strip() for s in open(args.samples_order) if s.strip()]
    sample_to_col = {s:i+4 for i,s in enumerate(samples_order)}
    anc_col_idx = [sample_to_col[s] for s in ancient_ids if s in sample_to_col]
    if not anc_col_idx:
        raise SystemExit("[ERROR] no ancient IDs found in samples.order")

    accum = {tid: {'num':0.0, 'den':0.0} for tid in ancient_ids}

    usecols = list(range(0, 4)) + list(range(4, 4+len(samples_order)))
    reader = pd.read_csv(args.gt_tsv, sep='\t', header=None, usecols=usecols,
                         dtype=str, chunksize=args.chunksize, low_memory=True)

    chunk_k = 0
    for chunk in reader:
        chunk_k += 1
        sys.stderr.write(f"[INFO] IRR pass, chunk {chunk_k}, rows={len(chunk)}\n")

        # 古树 ALT 计数
        A = chunk.iloc[:, anc_col_idx].map(alt_count).to_numpy(dtype=float)  # (nrow, n_anc)
        keys = (chunk.iloc[:,0].astype(str)+'|'+chunk.iloc[:,1].astype(str)+'|'
                +chunk.iloc[:,2].astype(str)+'|'+chunk.iloc[:,3].astype(str)).values

        meta = meta_df.reindex(keys)
        ok_mask = meta['w'].notna().values
        if ok_mask.sum() == 0:
            continue

        anc_count_v = meta['anc_count'].values
        cover_v     = meta[irr_cov_flag].values.astype(float)  # 0/1
        w_v         = meta['w'].values

        # anc_count 数值化
        anc_count_num = pd.to_numeric(pd.Series(anc_count_v), errors='coerce').fillna(0).to_numpy()
        rare_mask = (anc_count_num <= args.max_occ)

        for j, tid in enumerate(ancient_ids):
            carr = A[:, j] > 0
            mask = carr & ok_mask & rare_mask
            if not np.any(mask):
                continue
            w_sel = w_v[mask]
            c_sel = cover_v[mask]  # 0/1
            accum[tid]['num'] += float((w_sel * (1.0 - c_sel)).sum())
            accum[tid]['den'] += float(w_sel.sum())

    rows = []
    for tid in ancient_ids:
        num = accum[tid]['num']
        den = accum[tid]['den']
        irr = num
        irr_norm = (num/den) if den > 0 else np.nan
        # 标注子组（便于后续分组可视化）
        grp = 'anc_nat26' if tid in anc_nat_ids else ('anc_cult38' if tid in anc_cult_ids else 'anc_other')
        rows.append({'id': tid, 'group': grp, 'IRR': irr, 'IRR_norm01': irr_norm})

    irr_df = pd.DataFrame(rows).sort_values(['group','IRR_norm01'], ascending=[True, False])
    irr_df.to_csv(f"{args.out_prefix}.irr_per_tree.csv", index=False)
    sys.stderr.write(f"[OK] IRR done -> {args.out_prefix}.irr_per_tree.csv (coverage base: {args.irr_coverage})\n")

if __name__ == '__main__':
    main()

