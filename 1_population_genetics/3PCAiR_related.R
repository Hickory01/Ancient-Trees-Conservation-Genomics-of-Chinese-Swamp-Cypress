#!/usr/bin/env Rscript
### =========================================================
### PC-AiR：KING-robust → PC-AiR；输出 2D + 3D（Overall）
### - 2D：Overall / Zhujiang / Minjiang（PC1 vs PC2，轴标题含%）
### - 3D：Overall（PC1–PC3，轴标题含%）
### - 图片：PNG，600 dpi，6×5 cm
### =========================================================

suppressPackageStartupMessages({
  library(SNPRelate)
  library(GENESIS)
  library(GWASTools)
  library(parallel)
  library(ggplot2)
  library(scatterplot3d)   # 3D 散点
})

## ========= 1) 输入文件 =========
bed.fn <- "Gpen147_mis20.bed"
bim.fn <- "Gpen147_mis20.bim"
fam.fn <- "Gpen147_mis20.fam"
gdsfile <- "Gpen147_mis20.gds"

## 群组 ID 列表（逐行一个样本 ID）
read_idlist <- function(f) unique(scan(f, what = character(), quiet = TRUE))
admixed_ids   <- read_idlist("admixed7.list")
zhujiang_ids  <- read_idlist("Zhujiang48.list")
minjiang_ids  <- read_idlist("Minjiang92.list")
ancient_ids   <- read_idlist("ancient64.list")
wild_ids      <- read_idlist("wild33.list")
plant_ids     <- read_idlist("plant50.list")   # Cultivated

## ========= 2) PLINK → GDS =========
snpgdsBED2GDS(bed.fn = bed.fn, bim.fn = bim.fn, fam.fn = fam.fn, out.gdsfn = gdsfile)
snpgdsSummary(gdsfile)

## ========= 3) KING-robust（不做 LD，不限常染）=========
gds <- snpgdsOpen(gdsfile)
all_snps <- read.gdsn(index.gdsn(gds, "snp.id"))
ibd <- snpgdsIBDKING(
  gdsobj         = gds,
  sample.id      = NULL,
  snp.id         = all_snps,
  autosome.only  = FALSE,
  remove.monosnp = TRUE,
  type           = "KING-robust",
  verbose        = TRUE
)
KINGmat <- ibd$kinship
rownames(KINGmat) <- colnames(KINGmat) <- ibd$sample.id
snpgdsClose(gds)

## ========= 4) PC-AiR =========
geno_reader <- GdsGenotypeReader(filename = gdsfile)
genoData <- GenotypeData(geno_reader)
all_snps <- GWASTools::getSnpID(genoData)

mypcair <- pcair(
  gdsobj        = genoData,
  kinobj        = KINGmat,
  divobj        = KINGmat,
  snp.include   = all_snps,
  autosome.only = FALSE,
  num.cores     = max(1, parallel::detectCores(logical = FALSE)),
  verbose       = TRUE
)

## 提取 PC 得分与方差解释比例
pcs <- as.data.frame(mypcair$vectors)
colnames(pcs) <- paste0("PC", seq_len(ncol(pcs)))
pcs$ID <- rownames(mypcair$vectors)
pcs <- pcs[, c("ID", paste0("PC", seq_len(ncol(mypcair$vectors))))]

eigvals <- mypcair$values
propVar <- eigvals / sum(eigvals)
xlab_pc <- sprintf("PC1 (%.2f%%)", 100 * propVar[1])
ylab_pc <- sprintf("PC2 (%.2f%%)", 100 * propVar[2])
zlab_pc <- sprintf("PC3 (%.2f%%)", 100 * propVar[3])

## 保存 PC 表
write.table(pcs, "PC-AiR_scores.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## ========= 5) 标注 Type / 拆分子集 =========
pcs$Type <- ifelse(pcs$ID %in% ancient_ids, "Ancient",
                   ifelse(pcs$ID %in% wild_ids, "Wild",
                          ifelse(pcs$ID %in% plant_ids, "Cultivated", NA)))
pcs$Type <- factor(pcs$Type, levels = c("Ancient","Wild","Cultivated"))

overall_df  <- pcs
zhujiang_df <- subset(pcs, ID %in% zhujiang_ids)
minjiang_df <- subset(pcs, ID %in% minjiang_ids)

## ========= 6) 统一美学 =========
cols_type <- c(Ancient = "#ff7f00", Wild = "#6299ff", Cultivated = "#2e963a")
border_lwd <- 0.25

theme_uni <- theme_bw(base_size = 4) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = border_lwd),
    axis.line        = element_blank(),
    axis.title       = element_text(size = 4, color = "black"),
    axis.text        = element_text(size = 4, color = "black"),
    legend.position  = "none",
    plot.margin      = margin(1, 1, 1, 1, "mm")
  )

## ========= 7) 2D 图（PC1 vs PC2，含百分比轴标题）=========
p_overall_2d <- ggplot(overall_df, aes(PC1, PC2, color = Type)) +
  geom_point(size = 0.9, alpha = 0.85, stroke = 0) +
  scale_color_manual(values = cols_type, drop = FALSE) +
  labs(x = xlab_pc, y = ylab_pc) +
  theme_uni

p_zhujiang_2d <- ggplot(zhujiang_df, aes(PC1, PC2, color = Type)) +
  geom_point(size = 0.9, alpha = 0.85, stroke = 0) +
  scale_color_manual(values = cols_type, drop = FALSE) +
  labs(x = xlab_pc, y = ylab_pc) +
  theme_uni

p_minjiang_2d <- ggplot(minjiang_df, aes(PC1, PC2, color = Type)) +
  geom_point(size = 0.9, alpha = 0.85, stroke = 0) +
  scale_color_manual(values = cols_type, drop = FALSE) +
  labs(x = xlab_pc, y = ylab_pc) +
  theme_uni

## 保存 2D PNG（600 dpi，6×5 cm）
ggsave("PC-AiR_overall_2D.png",   p_overall_2d,  width = 6, height = 5, units = "cm", dpi = 600)
ggsave("PC-AiR_zhujiang_2D.png",  p_zhujiang_2d, width = 6, height = 5, units = "cm", dpi = 600)
ggsave("PC-AiR_minjiang_2D.png",  p_minjiang_2d, width = 6, height = 5, units = "cm", dpi = 600)

## ========= 8) Overall 的 3D 图（plot3D + jitter + 灰轴 + 黑框 + 细线）=========
suppressPackageStartupMessages({
  library(plot3D)
})

# 提取绘图数据
overall_3d <- overall_df[complete.cases(overall_df[, c("PC1","PC2","PC3","Type")]), ]
type_levels  <- levels(pcs$Type)
type_palette <- cols_type[type_levels]
myColors     <- cols_type[as.character(overall_3d$Type)]

## —— 添加轻微扰动（防止点重叠）——
set.seed(2025)
jitter_frac <- 0.025   # 扰动比例，可 0.005–0.015
dx <- diff(range(overall_3d$PC1)); if (!is.finite(dx) || dx == 0) dx <- 1
dy <- diff(range(overall_3d$PC2)); if (!is.finite(dy) || dy == 0) dy <- 1
dz <- diff(range(overall_3d$PC3)); if (!is.finite(dz) || dz == 0) dz <- 1

PC1_j <- overall_3d$PC1 + rnorm(nrow(overall_3d), 0, dx * jitter_frac)
PC2_j <- overall_3d$PC2 + rnorm(nrow(overall_3d), 0, dy * jitter_frac)
PC3_j <- overall_3d$PC3 + rnorm(nrow(overall_3d), 0, dz * jitter_frac)

## —— 导出 PNG：600 dpi，6×5 cm —— 
png("PC-AiR_overall_3D_plot3D_jitter_v4.png",
    width = 10, height = 9, units = "cm", res = 600, bg = "white")
par(mar = c(2, 2, 0.5, 0.5), family = "sans")  # 紧凑边距 + sans 字体

## —— 绘图 —— 
plot3D::scatter3D(
  x = PC1_j,
  y = PC2_j,
  z = PC3_j,
  pch = 21,              # 实心圆
  cex = 0.6,             # 点大小
  col = "black",         # 点边框为黑
  lwd = 0.1,            # 边框更细
  bg  = myColors,        # 填充色 = type 颜色
  xlab = xlab_pc, ylab = ylab_pc, zlab = zlab_pc,
  ticktype = "simple",   # 简洁刻度
  bty = "b2",            # 黑框立体盒
  box = TRUE,            # 显示坐标轴框
  theta = 45,            # 旋转角度
  phi = 35,              # 俯仰角
  d = 10,                # 透视深度
  colkey = FALSE,        # 不显示颜色条
  cex.axis = 0.5,        # 刻度 5pt
  cex.lab  = 0.5,        # 标签 5pt
  col.panel = "white",  # 面板底色更浅
  col.grid  = "grey99",  # 网格浅灰
  lwd.grid  = 0.005,       # 网格细线
  col.axis  = "grey40",  # 坐标轴（箭头）灰色
  lwd.axis  = 0.1,       # 坐标轴线更细
  col.lab   = "black",   # 轴标题文字黑色
  col.box   = "black",   # 立体盒子黑边框
  lwd.box   = 0.1       # 盒子边线细一点
)

## —— 底部水平图例 —— 
legend("bottom",
       legend = type_levels,
       pch = 21,
       pt.cex = 0.9,
       cex = 0.7,
       pt.bg = type_palette,
       col = "black",
       lwd = 0.25,
       bg = "white",
       bty = "n",
       horiz = TRUE)

dev.off()



