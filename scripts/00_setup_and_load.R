# 00_setup_and_load.R —— 读两队列 proteinGroups → 去污染/反库/仅位点 → 抽 LFQ →
#   0→NA → 合并 GU+TU 两队列(按基因名取交集) → 关联疾病分组与队列(批次) → 存原始矩阵
# 输入: data/GU/proteinGroups.txt, data/TU/proteinGroups.txt
# 输出: results/tables/00_raw_lfq_matrix.csv, 00_sample_metadata.csv, 00_protein_annotation.csv
# 注: 图内文字一律英文; 行名用人类真基因名(GN= 抽取, 对照 04 只能用 Phatr3 ID)。

.this <- (function(){a<-commandArgs(FALSE);f<-sub("^--file=","",a[grepl("^--file=",a)])
  if(length(f))return(dirname(normalizePath(f[1],winslash="/")))
  for(i in sys.nframe():1){o<-sys.frame(i)$ofile;if(!is.null(o))return(dirname(normalizePath(o,winslash="/")))};"."})()
suppressWarnings(suppressMessages(source(file.path(.this, "_common.R"))))

main <- function() {
  cat("=== 00 读取两队列 + 合并 ===\n")
  # 1) 各队列加载(已过滤污染/反库, 0->NA, 抽基因名)
  cohorts <- lapply(names(PG_FILES), function(coh) {
    r <- load_one_cohort(PG_FILES[[coh]])
    # 行名加队列前缀防同名基因冲突前先记录原基因
    list(coh=coh, X=r$X, anno=r$anno)
  })
  names(cohorts) <- names(PG_FILES)

  # 2) 以基因名为键合并; 仅保留两队列都检出的基因(交集), 行名=基因
  gene_of <- lapply(cohorts, function(c) c$anno$gene)
  shared <- Reduce(intersect, lapply(gene_of, function(g) unique(g[g != ""])))
  cat(sprintf("   共享基因(交集): %d\n", length(shared)))

  # 每队列: 同基因多行(蛋白组)按均值合并到基因层
  collapse_to_gene <- function(cobj, genes) {
    X <- cobj$X; g <- cobj$anno$gene
    mat <- sapply(genes, function(gn) {
      idx <- which(g == gn)
      if (length(idx) == 1) X[idx, ] else colMeans(X[idx, , drop=FALSE], na.rm=TRUE)
    })
    t(mat)  # 基因 x 样本
  }
  Xg <- lapply(cohorts, collapse_to_gene, genes=shared)
  # NaN(整组全NA取均值) -> NA
  Xg <- lapply(Xg, function(m){ m[is.nan(m)] <- NA; m })

  # 3) 横向拼接两队列样本
  X <- do.call(cbind, Xg)
  cat(sprintf("   合并矩阵: %d 基因 x %d 样本\n", nrow(X), ncol(X)))

  # 4) 样本元数据: 疾病分组(从样本名 HCC_/CIRR_) + 队列(批次)
  samp <- colnames(X)
  cohort <- rep(names(Xg), sapply(Xg, ncol))
  group  <- ifelse(grepl("HCC", samp), CASE, CONTROL)
  meta <- data.frame(sample=samp, group=group, cohort=cohort,
                     row.names=samp, stringsAsFactors=FALSE)
  cat("   分组 x 队列:\n"); print(table(meta$group, meta$cohort))

  # 5) 蛋白注释(以 GU 为准, 补 TU 独有)
  anno_all <- do.call(rbind, lapply(cohorts, function(c) c$anno))
  anno <- anno_all[!duplicated(anno_all$gene) & anno_all$gene %in% shared, ]
  rownames(anno) <- anno$gene; anno <- anno[shared, ]

  # 6) 落盘
  save_table(as.data.frame(X), "00_raw_lfq_matrix")
  save_table(meta, "00_sample_metadata", row.names=FALSE)
  save_table(anno, "00_protein_annotation")
  cat("=== 00 完成 ===\n")
}

if (sys.nframe() == 0) main()
