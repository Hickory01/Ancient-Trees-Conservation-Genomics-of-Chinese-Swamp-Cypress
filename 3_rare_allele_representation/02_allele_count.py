#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
02_allele_count.py  (no-unrelated + ancient sublineages)

从 all147.gt.tsv（由 bcftools query 导出：'%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n'）分块读取，
对给定的样本分组（必需：ancients / cultivated / wild；可选：anc_nat26 / anc_cult38 /
anc_admix4 / anc_Min43 / anc_Zhu17）计算每个 ALT 的：
  - in_*（该群体是否至少一位携带 ALT）
  - *_count（携带 ALT 的人数：ALT计数>0 的样本数）
  - AC_* / AN_*（ALT 等位总数 / 分母 = 2 * 非缺失基因型数）
  - fa_* = AC_* / AN_*（等位基因频率）
并输出“等位元素表”（每行一个 ALT）。

用法示例：
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
"""

import argparse
import sys
import os
import pandas as pd
import numpy as np

def parse_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--gt-tsv', required=True, help="由 bcftools query 导出的全体GT表（%CHROM %POS %REF %ALT [GT×N]）")
    p.add_argument('--samples-order', required=True, help="VCF中的样本顺序（bcftools query -l 导出）")

    # 必需主群
    p.add_argument('--ancients',   required=True, help="ancients64.list")
    p.add_argument('--cultivated', required=True, help="cult50.list")
    p.add_argument('--wild',       required=True, help="wild33.list")

    # 可选子群：历史与自然
    p.add_argument('--anc-nat',  default=None, help="可选：自然孑遗26名单（anc_nat26.list）")
    p.add_argument('--anc-cult', default=None, help="可选：历史人为栽培38名单（anc_cult38.list）")

    # 可选子群：遗传种群（新增）
    p.add_argument('--anc-admix', default=None, help="可选：古树-混合谱系（anc_admix4.list）")
    p.add_argument('--anc-min',   default=None, help="可选：古树-Minjiang（anc_Min43.list）")
    p.add_argument('--anc-zhu',   default=None, help="可选：古树-Zhujiang（anc_Zhu17.list）")

    p.add_argument('--out', default='allele_table.with_flags.tsv', help="输出文件（TSV）")
    p.add_argument('--chunksize', type=int, default=200000, help="分块大小（行）")
    p.add_argument('--sep', default='\t', help="输入文件分隔符")
    return p.parse_args()

def read_list(path):
    s = set()
    if path is None:
        return s
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                s.add(line)
    return s

def alt_count(gt: str):
    """
    将GT转换为ALT=1的等位计数：
      0/0, 0|0 -> 0
      0/1, 1/0, 0|1, 1|0 -> 1
      1/1, 1|1 -> 2
    含缺失等位（任一为 '.'）或完全缺失（'.', './.', '.|.'） -> np.nan（整个位点当缺失，不计入AN）
    兼容单等位字符串：'0' -> 0, '1' -> 1
    """
    if gt is None:
        return np.nan
    gt = str(gt).strip()
    if gt == '' or gt == '.' or gt == './.' or gt == '.|.':
        return np.nan

    gt = gt.replace('|', '/')
    parts = gt.split('/')

    # 单等位（少见；当作haploid或异常编码）
    if len(parts) == 1:
        a = parts[0]
        if a == '.':
            return np.nan
        return 1.0 if a == '1' else 0.0

    a, b = parts[0], parts[1]
    # 任一等位缺失 -> 整个基因型视作缺失
    if a == '.' or b == '.':
        return np.nan

    c = 0
    if a == '1':
        c += 1
    if b == '1':
        c += 1
    return float(c)

def build_col_index(samples_order, group_set):
    """
    给定 VCF 样本顺序（list[str]）和一个组（set[str]），
    返回该组样本在 all147.gt.tsv 中的列索引列表（注意：前4列是 CHR POS REF ALT，从第5列开始是GT）。
    """
    sample_to_col = {s: i+4 for i, s in enumerate(samples_order)}
    cols = []
    missing = []
    for s in group_set:
        if s in sample_to_col:
            cols.append(sample_to_col[s])
        else:
            missing.append(s)
    if missing:
        sys.stderr.write(f"[WARN] {len(missing)} IDs not found in samples.order: {missing[:5]}{' ...' if len(missing)>5 else ''}\n")
    cols.sort()
    return cols

def group_stats(gtm_float: np.ndarray, col_idx: list[int]):
    """
    对一个群体（列索引集合）在当前chunk计算：
      AC = ALT等位总数
      AN = 2 * 非缺失基因型数
      carriers = 携带ALT的“人数”（ALT计数>0）
    其中 gtm_float 是当前chunk的 GT 矩阵（numpy浮点数组；NaN 表示缺失），行是位点，列是样本。
    col_idx 使用 all147.gt.tsv 的列号，从 4 起（第5列是样本1）。
    """
    if not col_idx:
        # 该群体未提供名单
        nrow = gtm_float.shape[0]
        return (np.zeros(nrow), np.zeros(nrow), np.zeros(nrow))

    sub = gtm_float[:, [i-4 for i in col_idx]]  # shift，因为 gtm 从第5列开始
    ok = ~np.isnan(sub)
    AC = np.nansum(sub, axis=1)
    AN = 2.0 * ok.sum(axis=1)
    carriers = np.nansum((sub > 0).astype(float), axis=1)
    return AC, AN, carriers

def main():
    args = parse_args()

    # 读取样本顺序
    samples_order = [s.strip() for s in open(args.samples_order) if s.strip()]
    n_samples = len(samples_order)
    sys.stderr.write(f"[INFO] Loaded {n_samples} samples from {args.samples_order}\n")

    # 读取各群体名单
    anc_all = read_list(args.ancients)
    cult    = read_list(args.cultivated)
    wild    = read_list(args.wild)

    # 可选子群
    anc_nat   = read_list(args.anc_nat)  if args.anc_nat  else set()
    anc_cul   = read_list(args.anc_cult) if args.anc_cult else set()
    anc_admix = read_list(args.anc_admix) if args.anc_admix else set()
    anc_min   = read_list(args.anc_min)   if args.anc_min   else set()
    anc_zhu   = read_list(args.anc_zhu)   if args.anc_zhu   else set()

    # 构建列索引
    anc_cols   = build_col_index(samples_order, anc_all)
    cult_cols  = build_col_index(samples_order, cult)
    wild_cols  = build_col_index(samples_order, wild)

    nat_cols   = build_col_index(samples_order, anc_nat)   if anc_nat   else []
    acul_cols  = build_col_index(samples_order, anc_cul)   if anc_cul   else []
    admix_cols = build_col_index(samples_order, anc_admix) if anc_admix else []
    min_cols   = build_col_index(samples_order, anc_min)   if anc_min   else []
    zhu_cols   = build_col_index(samples_order, anc_zhu)   if anc_zhu   else []

    # 输入列：前4列 + 全部样本GT列
    usecols = list(range(0, 4)) + list(range(4, 4 + n_samples))

    # 分块读取
    reader = pd.read_csv(
        args.gt_tsv,
        sep=args.sep,
        header=None,
        usecols=usecols,
        chunksize=args.chunksize,
        dtype=str,
        low_memory=True
    )

    out_path = args.out
    # 若已存在旧文件，先删
    if os.path.exists(out_path):
        os.remove(out_path)

    chunk_idx = 0
    for chunk in reader:
        chunk_idx += 1
        sys.stderr.write(f"[INFO] Processing chunk #{chunk_idx}, rows={len(chunk)}\n")

        # 将 GT 列映射为 ALT=1 的等位计数（float；NaN表示缺失）
        gtm = chunk.iloc[:, 4:].map(alt_count).to_numpy(dtype=float)

        # 各群体统计
        AC_anc,  AN_anc,  N_anc  = group_stats(gtm, anc_cols)
        AC_cul,  AN_cul,  N_cul  = group_stats(gtm, cult_cols)
        AC_wld,  AN_wld,  N_wld  = group_stats(gtm, wild_cols)

        # 可选子群体
        if nat_cols:
            AC_nat, AN_nat, N_nat = group_stats(gtm, nat_cols)
        if acul_cols:
            AC_acul, AN_acul, N_acul = group_stats(gtm, acul_cols)
        if admix_cols:
            AC_admix, AN_admix, N_admix = group_stats(gtm, admix_cols)
        if min_cols:
            AC_min, AN_min, N_min = group_stats(gtm, min_cols)
        if zhu_cols:
            AC_zhu, AN_zhu, N_zhu = group_stats(gtm, zhu_cols)

        # 组装输出 DataFrame
        df = pd.DataFrame({
            'CHR': chunk.iloc[:, 0].values,
            'POS': chunk.iloc[:, 1].values,
            'REF': chunk.iloc[:, 2].values,
            'ALT': chunk.iloc[:, 3].values,

            'anc_count':  N_anc,
            'cult_count': N_cul,
            'wild_count': N_wld,
            'in_anc':  (N_anc  > 0).astype(int),
            'in_cult': (N_cul  > 0).astype(int),
            'in_wild': (N_wld  > 0).astype(int),

            'AC_anc':  AC_anc,  'AN_anc':  AN_anc,
            'AC_cult': AC_cul,  'AN_cult': AN_cul,
            'AC_wild': AC_wld,  'AN_wild': AN_wld,
        })

        # 可选子群体输出（自然孑遗/历史栽培）
        if nat_cols:
            df['anc_nat_count'] = N_nat
            df['in_anc_nat']    = (N_nat > 0).astype(int)
            df['AC_anc_nat']    = AC_nat
            df['AN_anc_nat']    = AN_nat
        if acul_cols:
            df['anc_cult_count'] = N_acul
            df['in_anc_cult']    = (N_acul > 0).astype(int)
            df['AC_anc_cult']    = AC_acul
            df['AN_anc_cult']    = AN_acul

        # 可选子群体输出（谱系：admix / Minjiang / Zhujiang）
        if admix_cols:
            df['anc_admix_count'] = N_admix
            df['in_anc_admix']    = (N_admix > 0).astype(int)
            df['AC_anc_admix']    = AC_admix
            df['AN_anc_admix']    = AN_admix
        if min_cols:
            df['anc_min_count'] = N_min
            df['in_anc_min']    = (N_min > 0).astype(int)
            df['AC_anc_min']    = AC_min
            df['AN_anc_min']    = AN_min
        if zhu_cols:
            df['anc_zhu_count'] = N_zhu
            df['in_anc_zhu']    = (N_zhu > 0).astype(int)
            df['AC_anc_zhu']    = AC_zhu
            df['AN_anc_zhu']    = AN_zhu

        # 频率列（避免除0：AN=0 -> fa=NaN）
        base_groups = ['anc', 'cult', 'wild']
        for g in base_groups:
            num = df[f'AC_{g}'].astype(float)
            den = df[f'AN_{g}'].astype(float)
            with np.errstate(divide='ignore', invalid='ignore'):
                df[f'fa_{g}'] = np.where(den > 0, num/den, np.nan)

        # 子群体频率（若存在）
        if nat_cols:
            with np.errstate(divide='ignore', invalid='ignore'):
                df['fa_anc_nat'] = np.where(df['AN_anc_nat'] > 0, df['AC_anc_nat'] / df['AN_anc_nat'], np.nan)
        if acul_cols:
            with np.errstate(divide='ignore', invalid='ignore'):
                df['fa_anc_cult'] = np.where(df['AN_anc_cult'] > 0, df['AC_anc_cult'] / df['AN_anc_cult'], np.nan)
        if admix_cols:
            with np.errstate(divide='ignore', invalid='ignore'):
                df['fa_anc_admix'] = np.where(df['AN_anc_admix'] > 0, df['AC_anc_admix'] / df['AN_anc_admix'], np.nan)
        if min_cols:
            with np.errstate(divide='ignore', invalid='ignore'):
                df['fa_anc_min'] = np.where(df['AN_anc_min'] > 0, df['AC_anc_min'] / df['AN_anc_min'], np.nan)
        if zhu_cols:
            with np.errstate(divide='ignore', invalid='ignore'):
                df['fa_anc_zhu'] = np.where(df['AN_anc_zhu'] > 0, df['AC_anc_zhu'] / df['AN_anc_zhu'], np.nan)

        # 写出（首块写表头，后续追加）
        mode = 'w' if chunk_idx == 1 else 'a'
        header = True if chunk_idx == 1 else False
        df.to_csv(out_path, sep='\t', index=False, mode=mode, header=header)

    sys.stderr.write(f"[DONE] Wrote output to: {out_path}\n")

if __name__ == '__main__':
    main()

