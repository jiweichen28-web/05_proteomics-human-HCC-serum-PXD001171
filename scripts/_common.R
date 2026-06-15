# _common.R —— 项目共享配置与函数 (R 版)
# 人类血清蛋白质组: 肝细胞癌(HCC) vs 肝硬化 (PRIDE PXD001171, MaxQuant LFQ)
# 严格遵循固定绘图规范: 8pt 字体, theme_pub_bw, 每图 PDF + 300ppi PNG, COL 配色
# 图内文字一律英文(文件名+标签), 代码注释中文; 路径全部相对项目根
# 本文件自动定位项目根并 setwd, 故从任意目录运行均可:
#   Rscript scripts/00_setup_and_load.R

suppressMessages({library(ggplot2); library(data.table); library(matrixStats)})

# ── 自动定位项目根: 取本脚本路径 -> 上溯一级(scripts -> 根), 再 setwd ──
.find_proj_root <- function() {
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f) == 0) {
    for (i in sys.nframe():1) { of <- sys.frame(i)$ofile; if (!is.null(of)) { f <- of; break } }
  }
  if (length(f) == 0 || is.na(f[1]) || f[1] == "") return(NULL)
  d <- dirname(normalizePath(f[1], winslash="/", mustWork=FALSE))
  normalizePath(file.path(d, ".."), winslash="/", mustWork=FALSE)
}
.root <- .find_proj_root()
if (!is.null(.root) && dir.exists(file.path(.root, "data"))) setwd(.root)

# ── 输出目录 ──
FIG_DIR <- "results/figures"; TAB_DIR <- "results/tables"
for (d in c(FIG_DIR, TAB_DIR)) if (!dir.exists(d)) dir.create(d, recursive=TRUE)

# ── 数据集常量 ──
ACCESSION <- "PXD001171"
SPECIES   <- "Homo sapiens"
# MaxQuant 结果: 两队列各一份 txt.zip 解压后置于 data/GU/ 和 data/TU/
# GU=Georgetown(美国, 57 HCC/59 cirrhosis), TU=Tanta(埃及, 40 HCC/49 cirrhosis)
PG_FILES  <- c(GU = "data/GU/proteinGroups.txt", TU = "data/TU/proteinGroups.txt")

# 实验分组: 疾病状态 HCC(肝癌) vs Cirrhosis(肝硬化, 对照)
# 分组来源: MaxQuant 样本名前缀 / experimentalDesign / 原文样本表(P2 校验后确定映射)
CONTROL <- "Cirrhosis"; CASE <- "HCC"
# 图内英文标签
GROUP_LABEL_EN <- c("Cirrhosis"="Cirrhosis", "HCC"="HCC")
# 队列(批次)标签
COHORT_LABEL <- c("GU"="Georgetown (US)", "TU"="Tanta (EG)")

# MaxQuant 过滤标志列 (值为 "+" 表示命中, 需剔除)
# 注: 本数据集用旧版 MaxQuant(1.4), 污染物列名是 "Contaminant" 而非 "Potential contaminant"
FILTER_FLAGS <- c("Reverse", "Contaminant", "Potential contaminant", "Only identified by site")

# 差异判定阈值
FC_THR <- log2(1.5); FDR_THR <- 0.05

# ── 绘图风格 (8pt, 黑色描边, 无网格; 图内文字一律英文) ──
FONT_SIZE <- 8
PT <- FONT_SIZE / 2.845          # geom 字体 pt -> mm
LW <- 0.5 / 2.1333               # 线宽 0.5pt -> linewidth
PLOT_FAMILY <- "sans"            # 图内英文用无衬线
# 2 组对比配色: 对照(肝硬化)=深蓝, 病例(HCC)=正红
COL_GRP <- c("Cirrhosis"="#00468A", "HCC"="#EC0000")
# 队列配色
COL_COHORT <- c("GU"="#42B540", "TU"="#925E9F")
# 固定配色板(多类别取色)
COL7 <- c("#00468A","#EC0000","#42B540","#0099B4","#925E9F","#FCAE91","#AC002A")

theme_pub_bw <- function() {
  theme_bw(base_size=FONT_SIZE) +
    theme(
      text         =element_text(size=FONT_SIZE, color="black", face="plain", family=PLOT_FAMILY),
      plot.title   =element_text(size=FONT_SIZE, color="black", face="plain", family=PLOT_FAMILY, hjust=0.5),
      axis.title   =element_text(size=FONT_SIZE, color="black", face="plain", family=PLOT_FAMILY),
      axis.text    =element_text(size=FONT_SIZE, color="black", face="plain", family=PLOT_FAMILY),
      legend.title =element_text(size=FONT_SIZE, color="black", face="plain", family=PLOT_FAMILY),
      legend.text  =element_text(size=FONT_SIZE, color="black", face="plain", family=PLOT_FAMILY),
      strip.text   =element_text(size=FONT_SIZE, color="black", face="plain", family=PLOT_FAMILY),
      panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
      panel.border =element_rect(color="black", fill=NA, linewidth=LW),
      axis.line    =element_blank(),
      axis.ticks   =element_line(color="black", linewidth=LW),
      panel.background =element_rect(fill="transparent", color=NA),
      plot.background  =element_rect(fill="transparent", color=NA),
      legend.background=element_rect(fill="transparent", color=NA),
      legend.key   =element_rect(fill="transparent", color=NA))
}

# 同时保存 PDF 矢量 + PNG 300ppi (图内英文, 普通设备即可)
savefig <- function(p, name, w=5, h=4) {
  ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), p, width=w, height=h)
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), p, width=w, height=h, dpi=300)
  cat(sprintf("   [fig] %s.pdf / .png\n", name))
}

save_table <- function(df, name, row.names=TRUE) {
  path <- file.path(TAB_DIR, paste0(name, ".csv"))
  write.csv(df, path, row.names=row.names, fileEncoding="UTF-8")
  cat(sprintf("   [tab] %s.csv  (%d x %d)\n", name, nrow(df), ncol(df)))
  invisible(path)
}

# ── 数据加载 ──

# 加载单个队列的 proteinGroups.txt -> 蛋白 x 样本 LFQ 矩阵 (已过滤污染/反库, 0->NA)
# 行名优先 Gene names(人类有真基因名), 缺失回退 Majority protein IDs
# 返回 list(X=矩阵, anno=注释)
load_one_cohort <- function(pg_file, quant_prefer = c("LFQ intensity ", "iBAQ ", "Intensity ")) {
  pg <- fread(pg_file, sep="\t", check.names=FALSE)
  n0 <- nrow(pg)
  for (col in FILTER_FLAGS) if (col %in% names(pg)) pg <- pg[get(col) != "+"]
  cat(sprintf("   [%s] 过滤污染/反库/仅位点: %d -> %d 蛋白组\n", basename(dirname(pg_file)), n0, nrow(pg)))
  prefix <- NULL; qcols <- NULL
  for (pref in quant_prefer) {
    cols <- grep(paste0("^", pref), names(pg), value=TRUE)
    cols <- setdiff(cols, trimws(pref))
    if (length(cols)) { prefix <- pref; qcols <- cols; break }
  }
  if (is.null(qcols)) stop("proteinGroups.txt 未找到定量列")
  samples <- sub(paste0("^", prefix), "", qcols)
  X <- as.matrix(sapply(pg[, qcols, with=FALSE], as.numeric))
  ids  <- sub(";.*", "", pg[["Majority protein IDs"]])
  # 基因名: 本旧版 MaxQuant 无 Gene names 列, 但 Fasta headers 是完整 UniProt 头,
  # 内含 "GN=XXX" -> 正则抽真基因名(人类相对 04 硅藻的关键优势: 蛋白可映射到基因)
  hdr <- if ("Fasta headers" %in% names(pg)) pg[["Fasta headers"]] else pg[["Majority protein IDs"]]
  gene <- sub(".*GN=(\\S+).*", "\\1", hdr)
  gene[!grepl("GN=", hdr)] <- ""            # 无 GN 标记的置空
  # UniProt accession(从 sp|ACC|NAME 抽), 供富集时映射 ENTREZ
  uniprot <- sub(".*?\\|([A-Z0-9]+)\\|.*", "\\1", ids)
  rn   <- ifelse(!is.na(gene) & gene != "", gene, uniprot)
  rownames(X) <- make.unique(rn); colnames(X) <- samples
  X[X == 0] <- NA
  anno <- data.frame(protein_group = pg[["Majority protein IDs"]],
                     gene = gene, uniprot = uniprot, n_peptides = pg[["Peptides"]],
                     row.names = rownames(X), stringsAsFactors=FALSE)
  list(X = X, anno = anno)
}

# Perseus 下移填补: 每样本缺失值从 N(mean-1.8sd, 0.3sd) 抽样 (模拟低丰度左删失)
perseus_impute <- function(X, width=0.3, downshift=1.8, seed=42) {
  set.seed(seed); Xi <- X
  for (s in seq_len(ncol(Xi))) {
    col <- Xi[, s]; obs <- col[!is.na(col)]
    na_idx <- is.na(col)
    if (any(na_idx) && length(obs) > 1) {
      mu <- mean(obs); sdv <- sd(obs)
      Xi[na_idx, s] <- rnorm(sum(na_idx), mu - downshift*sdv, width*sdv)
    }
  }
  Xi
}

# BH-FDR
bh_fdr <- function(p) p.adjust(p, method="BH")

# ── 手写 NIPALS PLS1 + VIP + LOO (无 ropls/mixOmics 依赖) ──
nipals_pls1 <- function(X, y, ncomp=2) {
  X <- as.matrix(X); y <- as.numeric(y); n <- nrow(X); p <- ncol(X)
  Xr <- X; yr <- y
  Tt <- matrix(0,n,ncomp); W <- matrix(0,p,ncomp); P <- matrix(0,p,ncomp); Q <- numeric(ncomp)
  for (a in 1:ncomp) {
    w <- t(Xr) %*% yr; w <- w/sqrt(sum(w^2))
    tt <- Xr %*% w; q <- sum(yr*tt)/sum(tt*tt); pp <- t(Xr) %*% tt/sum(tt*tt)
    Xr <- Xr - tt %*% t(pp); yr <- yr - q*tt
    Tt[,a] <- tt; W[,a] <- w; P[,a] <- pp; Q[a] <- q
  }
  list(scores=Tt, weights=W, loadings=P, q=Q, ncomp=ncomp)
}
pls_vip <- function(fit) {
  W <- fit$weights; Tt <- fit$scores; Q <- fit$q; p <- nrow(W); ncomp <- fit$ncomp
  ssy <- (Q^2)*colSums(Tt^2); total <- sum(ssy)
  sapply(1:p, function(i) {
    wgt <- sapply(1:ncomp, function(a) (W[i,a]/sqrt(sum(W[,a]^2)))^2)
    sqrt(p*sum(ssy*wgt)/total) })
}
pls_predict <- function(fit, Xnew) {
  B <- fit$weights %*% solve(t(fit$loadings) %*% fit$weights) %*% fit$q
  as.numeric(Xnew %*% B)
}
# 留一交叉验证准确率 (二分类: 是否 HCC)
plsda_loo_acc <- function(X, y, ncomp=2) {
  n <- nrow(X); correct <- 0
  for (i in 1:n) {
    mu <- colMeans(X[-i,]); sdv <- apply(X[-i,], 2, sd); sdv[sdv==0] <- 1
    Xtr <- scale(X[-i,], center=mu, scale=sdv)
    ym <- mean(y[-i]); fit <- nipals_pls1(Xtr, y[-i]-ym, ncomp)
    xte <- (X[i,]-mu)/sdv
    pred <- as.numeric(pls_predict(fit, matrix(xte,1)) + ym > 0.5)
    correct <- correct + (pred == y[i])
  }
  correct/n*100
}
