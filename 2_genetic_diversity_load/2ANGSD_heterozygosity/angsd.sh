#!/usr/bin/env bash
set -euo pipefail

# ============================
# Usage
# ============================
if [ $# -lt 3 ]; then
    echo "Usage: $0 <id_list> <chr_list> <ref.fa>"
    echo "Example: $0 id.list chr11.list ref.fa"
    exit 1
fi

ID_LIST=$1         # 个体ID列表文件（每行一个）
CHR_LIST=$2        # scaffold/chr 列表文件（每行一个）
REF=$3             # 参考基因组

OUTDIR="angsd_out"
mkdir -p "$OUTDIR"

# ============================
# 主循环
# ============================

# 检查 ID 列表
if [ ! -s "$ID_LIST" ]; then
    echo "❌ ERROR: ID list ($ID_LIST) is empty or missing."
    exit 1
fi

# 检查 Scaffold 列表
if [ ! -s "$CHR_LIST" ]; then
    echo "❌ ERROR: Chr list ($CHR_LIST) is empty or missing."
    exit 1
fi

while read -r id; do
    id=$(echo "$id" | xargs)
    [[ -z "$id" ]] && continue

    bam="${id}.bam"
    if [ ! -f "$bam" ]; then
        echo "⚠️ WARNING: BAM file not found for ID '$id' ($bam)" >&2
        continue
    fi

    while read -r scaf; do
        scaf=$(echo "$scaf" | xargs)
        [[ -z "$scaf" ]] && continue

        out_prefix="${OUTDIR}/${id}_${scaf}"
        echo ">>> Running ANGSD for ${id} on ${scaf} ..."

        angsd -i "$bam" \
              -ref "$REF" \
              -anc "$REF" \
              -r "$scaf" \
              -GL 1 \
              -doSaf 1 \
              -baq 1 \
              -C 50 \
              -minQ 20 \
              -minMapQ 30 \
              -remove_bads 1 \
              -only_proper_pairs 1 \
              -out "$out_prefix" &

    done < "$CHR_LIST"

done < "$ID_LIST"

# 等待所有任务完成
wait
echo "✅ 所有 scaffold × 个体 的 ANGSD 任务完成"

