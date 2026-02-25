#!/usr/bin/env bash
set -euo pipefail

# ============================
# Usage
# ============================
if [ $# -lt 3 ]; then
    echo "Usage: $0 <id_list> <chr_list> <ref.fa> [threads_for_realSFS]"
    echo "Example: $0 id.list chr11.list ref.fa 4"
    exit 1
fi

ID_LIST=$1            # 每行一个样本ID（与 .bam 的前缀一致，如 AF1A）
CHR_LIST=$2           # 每行一个 scaffold/chr 名称（与 .fai 一致）
REF=$3                # 引用，实质上 realSFS 不用，但留作一致性校验
P=${4:-4}             # realSFS -P 线程数（默认4）

OUTDIR="angsd_out"    # 与上一步 ANGSD 输出目录一致
[[ -d "$OUTDIR" ]] || { echo "❌ ERROR: $OUTDIR not found. Run ANGSD step first."; exit 1; }

SFS_SUMMARY="${OUTDIR}/heterozygosity_per_scaffold.tsv"
ID_SUMMARY="${OUTDIR}/heterozygosity_per_individual.tsv"

# 头文件
echo -e "id\tscaf\tsfs0\tsfs1\tsfs2\ttotal\thet" > "$SFS_SUMMARY"

# ============ 第1阶段：逐 (id,scaf) 跑 realSFS 产生 .sfs ============
echo ">>> Phase 1: running realSFS for each (id, scaffold) ..."
while read -r id; do
    id=$(echo "$id" | xargs); [[ -z "$id" ]] && continue

    while read -r scaf; do
        scaf=$(echo "$scaf" | xargs); [[ -z "$scaf" ]] && continue

        idx="${OUTDIR}/${id}_${scaf}.saf.idx"
        sfs="${OUTDIR}/${id}_${scaf}.sfs"

        if [[ ! -f "$idx" ]]; then
            echo "⚠️  WARNING: missing $idx ; skip $id $scaf" >&2
            continue
        fi

        # realSFS 输出一行 3 个数（单个体 1D-SFS：0/1/2）
        # 不限制并发：按需开启很多 realSFS；每个 realSFS 内部用 -P $P
        realSFS "$idx" -P "$P" > "$sfs" &

    done < "$CHR_LIST"
done < "$ID_LIST"

# 等待所有 realSFS 完成
wait
echo ">>> Phase 1 done."

# ============ 第2阶段：计算每对 (id,scaf) 的 H 并写入明细表 ============
echo ">>> Phase 2: computing heterozygosity per (id, scaffold) ..."
while read -r id; do
    id=$(echo "$id" | xargs); [[ -z "$id" ]] && continue

    while read -r scaf; do
        scaf=$(echo "$scaf" | xargs); [[ -z "$scaf" ]] && continue

        sfs="${OUTDIR}/${id}_${scaf}.sfs"
        [[ -f "$sfs" ]] || continue

        # 读取 S0 S1 S2；无论是概率型还是计数型，S1/sum 都成立
        read -r S0 S1 S2 < "$sfs" || continue

        # bc/awk 兼容科学计数法
        total=$(awk -v a="$S0" -v b="$S1" -v c="$S2" 'BEGIN{print a+b+c}')
        het=$(awk -v b="$S1" -v t="$total" 'BEGIN{ if(t>0){printf "%.10f", b/t}else{print "NA"} }')

        echo -e "${id}\t${scaf}\t${S0}\t${S1}\t${S2}\t${total}\t${het}" >> "$SFS_SUMMARY"

    done < "$CHR_LIST"
done < "$ID_LIST"

echo ">>> Per-scaffold heterozygosity written to: $SFS_SUMMARY"

# ============ 第3阶段：按个体汇总（把各 scaffold 的 SFS 相加后再算 H） ============
echo ">>> Phase 3: aggregating per individual ..."
echo -e "id\tsum_sfs0\tsum_sfs1\tsum_sfs2\ttotal\thet" > "$ID_SUMMARY"

while read -r id; do
    id=$(echo "$id" | xargs); [[ -z "$id" ]] && continue

    # 累加该 id 下所有 scaffold 的 SFS
    sums=$(awk -v ID="$id" 'BEGIN{s0=0;s1=0;s2=0}
        NR>1 && $1==ID { s0+= $3; s1+= $4; s2+= $5 }
        END{printf("%.10f %.10f %.10f\n", s0, s1, s2)}' "$SFS_SUMMARY")

    S0sum=$(echo "$sums" | awk '{print $1}')
    S1sum=$(echo "$sums" | awk '{print $2}')
    S2sum=$(echo "$sums" | awk '{print $3}')

    total=$(awk -v a="$S0sum" -v b="$S1sum" -v c="$S2sum" 'BEGIN{print a+b+c}')
    het=$(awk -v b="$S1sum" -v t="$total" 'BEGIN{ if(t>0){printf "%.10f", b/t}else{print "NA"} }')

    echo -e "${id}\t${S0sum}\t${S1sum}\t${S2sum}\t${total}\t${het}" >> "$ID_SUMMARY"
done < "$ID_LIST"

echo ">>> Per-individual heterozygosity written to: $ID_SUMMARY"
echo "✅ All done."

