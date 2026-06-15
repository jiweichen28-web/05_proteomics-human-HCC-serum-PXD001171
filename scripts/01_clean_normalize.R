# 01_clean_normalize.R —— 有效值过滤(组内≥50%检出) → log2 → 批次校正(ComBat, 两队列) →
#   Perseus 下移填补 → QC 图(检出数/箱线/填补密度/队列分离)
# 输入: results/tables/00_raw_lfq_matrix.csv, 00_sample_metadata.csv
# 输出: results/tables/01_clean_log2_imputed.csv, 01_clean_log2_observed.csv
#       results/figures/01_QC_*.{pdf,png}
# 注: 人血清是异质队列 + 两中心 → 批次效应突出, 需 ComBat(对照 04 等基因型可忽略批次)。

.this <- (function(){a<-commandArgs(FALSE);f<-sub("^--file=","",a[grepl("^--file=",a)])
  if(length(f))return(dirname(normalizePath(f[1],winslash="/")))
  for(i in sys.nframe():1){o<-sys.frame(i)$ofile;if(!is.null(o))return(dirname(normalizePath(o,winslash="/")))};"."})()
suppressWarnings(suppressMessages({source(file.path(.this, "_common.R")); library(sva)}))

main <- function() {
  cat("=== 01 清洗 / 批次校正 / 填补 ===\n")
  X <- as.matrix(read.csv(file.path(TAB_DIR, "00_raw_lfq_matrix.csv"), row.names=1, check.names=FALSE))
  meta <- read.csv(file.path(TAB_DIR, "00_sample_metadata.csv"), row.names=1)
  meta <- meta[colnames(X), ]
  n0 <- nrow(X)

  # 1) 有效值过滤: 每基因在至少一组内检出 >=50%
  grps <- unique(meta$group)
  valid <- sapply(grps, function(g) {
    cols <- rownames(meta)[meta$group == g]
    rowSums(!is.na(X[, cols, drop=FALSE])) >= 0.5*length(cols)
  })
  keep <- rowSums(valid) >= 1
  X <- X[keep, , drop=FALSE]
  cat(sprintf("   有效值过滤(组内>=50%%检出): %d -> %d 基因\n", n0, nrow(X)))

  # 2) log2 (LFQ 已 MaxLFQ 归一, 不再二次归一)
  Xlog <- log2(X)
  save_table(as.data.frame(Xlog), "01_clean_log2_observed")

  # 3) Perseus 下移填补(先填补再 ComBat: ComBat 不接受 NA)
  Ximp <- perseus_impute(Xlog)

  # 4) 批次校正: ComBat 去两队列(GU/TU)系统差异, 保护疾病分组
  mod <- model.matrix(~ group, data=meta)
  Xcb <- ComBat(dat=Ximp, batch=meta$cohort, mod=mod, par.prior=TRUE, prior.plots=FALSE)
  cat("   ComBat 批次校正完成(batch=队列, 保护 group)\n")
  save_table(as.data.frame(Xcb), "01_clean_log2_imputed")

  # ── QC1: 各样本检出基因数(填补前), 按队列上色 ──
  ndet <- colSums(!is.na(Xlog))
  d1 <- data.frame(sample=names(ndet), n=ndet,
                   cohort=meta[names(ndet),"cohort"], group=meta[names(ndet),"group"])
  p1 <- ggplot(d1, aes(reorder(sample,n), n, fill=cohort)) +
    geom_col(width=0.8) +
    scale_fill_manual(values=COL_COHORT, name="Cohort") +
    labs(x="Samples (sorted)", y="Proteins quantified",
         title="Quantified proteins per sample") +
    theme_pub_bw() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
  savefig(p1, "01_QC_valid_values", w=5.5, h=3.4)

  # ── QC2: log2 强度箱线(填补前), 按队列上色 ──
  long <- data.frame(value=as.vector(Xlog),
                     sample=rep(colnames(Xlog), each=nrow(Xlog)),
                     cohort=rep(meta$cohort, each=nrow(Xlog)))
  p2 <- ggplot(long, aes(sample, value, fill=cohort)) +
    geom_boxplot(outlier.size=0.2, linewidth=LW) +
    scale_fill_manual(values=COL_COHORT, name="Cohort") +
    labs(x="Samples", y="log2 LFQ intensity",
         title="Per-sample intensity (observed)") +
    theme_pub_bw() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
  savefig(p2, "01_QC_boxplot", w=5.5, h=3.4)

  # ── QC3: 填补值左移密度(观测 vs 填补) ──
  obs_v <- as.vector(Xlog[!is.na(Xlog)])
  imp_v <- as.vector(Ximp[is.na(Xlog)])
  d3 <- rbind(data.frame(value=obs_v, type="Observed"),
              data.frame(value=imp_v, type="Imputed"))
  p3 <- ggplot(d3, aes(value, fill=type, color=type)) +
    geom_density(alpha=0.4, linewidth=LW) +
    scale_fill_manual(values=c("Observed"="#00468A","Imputed"="#EC0000"), name=NULL) +
    scale_color_manual(values=c("Observed"="#00468A","Imputed"="#EC0000"), name=NULL) +
    labs(x="log2 intensity", y="Density",
         title="Perseus downshift imputation") + theme_pub_bw()
  savefig(p3, "01_QC_density_impute", w=4.5, h=3.4)

  # ── QC4: ComBat 前后 PCA, 看队列(批次)是否被消除 ──
  pca_df <- function(mat, tag) {
    pc <- prcomp(t(mat), center=TRUE, scale.=TRUE)
    ve <- round(100*pc$sdev^2/sum(pc$sdev^2), 1)
    data.frame(PC1=pc$x[,1], PC2=pc$x[,2], cohort=meta$cohort, group=meta$group,
               stage=tag, ve1=ve[1], ve2=ve[2])
  }
  d4 <- rbind(pca_df(Ximp, "Before ComBat"), pca_df(Xcb, "After ComBat"))
  p4 <- ggplot(d4, aes(PC1, PC2, color=cohort, shape=group)) +
    geom_point(size=1.5, alpha=0.8) +
    scale_color_manual(values=COL_COHORT, name="Cohort") +
    scale_shape_manual(values=c("Cirrhosis"=1,"HCC"=19), name="Group") +
    facet_wrap(~stage, scales="free") +
    labs(title="Batch effect: cohort before/after ComBat") + theme_pub_bw()
  savefig(p4, "01_QC_batch_pca", w=6.5, h=3.2)

  cat("=== 01 完成 ===\n")
}

if (sys.nframe() == 0) main()
