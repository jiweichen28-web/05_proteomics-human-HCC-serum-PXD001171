# 02_pca_plsda.R —— PCA(按疾病/队列分别上色, 排查批次是否主导) + PLS-DA(手写 NIPALS,
#   LOO 交叉验证 + VIP) → 得分图/Scree/PLS-DA 图
# 输入: results/tables/01_clean_log2_imputed.csv, 00_sample_metadata.csv
# 输出: results/tables/02_PCA_scores.csv, 02_PLSDA_VIP.csv
#       results/figures/02_PCA_scores.{pdf,png}(疾病/队列双版), 02_PCA_scree.*, 02_PLSDA_scores.*
# 注: 人群 PC1 常是性别/队列而非疾病, 必须按协变量上色看(对照 04 PC1 干净就是处理)。

.this <- (function(){a<-commandArgs(FALSE);f<-sub("^--file=","",a[grepl("^--file=",a)])
  if(length(f))return(dirname(normalizePath(f[1],winslash="/")))
  for(i in sys.nframe():1){o<-sys.frame(i)$ofile;if(!is.null(o))return(dirname(normalizePath(o,winslash="/")))};"."})()
suppressWarnings(suppressMessages(source(file.path(.this, "_common.R"))))

main <- function() {
  cat("=== 02 PCA + PLS-DA ===\n")
  X <- as.matrix(read.csv(file.path(TAB_DIR, "01_clean_log2_imputed.csv"), row.names=1, check.names=FALSE))
  meta <- read.csv(file.path(TAB_DIR, "00_sample_metadata.csv"), row.names=1)
  meta <- meta[colnames(X), ]

  # ── PCA (样本 x 基因, 标准化) ──
  pca <- prcomp(t(X), center=TRUE, scale.=TRUE)
  ve <- round(100*pca$sdev^2/sum(pca$sdev^2), 1)
  sc <- data.frame(sample=colnames(X), PC1=pca$x[,1], PC2=pca$x[,2], PC3=pca$x[,3],
                   group=meta$group, cohort=meta$cohort)
  save_table(sc, "02_PCA_scores", row.names=FALSE)
  cat(sprintf("   PCA 方差: PC1 %.1f%% / PC2 %.1f%% / PC3 %.1f%%\n", ve[1], ve[2], ve[3]))

  # 得分图 A: 按疾病上色
  pA <- ggplot(sc, aes(PC1, PC2, color=group)) +
    geom_point(size=1.8, alpha=0.8) +
    scale_color_manual(values=COL_GRP, labels=GROUP_LABEL_EN, name="Group") +
    geom_vline(xintercept=0, linetype="dashed", linewidth=LW, color="grey60") +
    geom_hline(yintercept=0, linetype="dashed", linewidth=LW, color="grey60") +
    labs(x=sprintf("PC1 (%.1f%%)", ve[1]), y=sprintf("PC2 (%.1f%%)", ve[2]),
         title="PCA scores - colour = disease group") + theme_pub_bw()
  savefig(pA, "02_PCA_scores", w=4.6, h=3.6)

  # 得分图 B: 按队列上色(排查批次是否残留)
  pB <- ggplot(sc, aes(PC1, PC2, color=cohort)) +
    geom_point(size=1.8, alpha=0.8) +
    scale_color_manual(values=COL_COHORT, labels=COHORT_LABEL, name="Cohort") +
    geom_vline(xintercept=0, linetype="dashed", linewidth=LW, color="grey60") +
    geom_hline(yintercept=0, linetype="dashed", linewidth=LW, color="grey60") +
    labs(x=sprintf("PC1 (%.1f%%)", ve[1]), y=sprintf("PC2 (%.1f%%)", ve[2]),
         title="PCA scores - colour = cohort") + theme_pub_bw()
  savefig(pB, "02_PCA_scores_cohort", w=4.6, h=3.6)

  # Scree
  scree <- data.frame(PC=factor(paste0("PC",1:10), levels=paste0("PC",1:10)), ve=ve[1:10])
  pS <- ggplot(scree, aes(PC, ve)) +
    geom_col(width=0.7, fill=COL_GRP[["HCC"]], color="black", linewidth=LW) +
    labs(x=NULL, y="Variance explained (%)", title="Scree plot") + theme_pub_bw()
  savefig(pS, "02_PCA_scree", w=4.2, h=3.0)

  # ── PLS-DA (响应=是否 HCC), 标准化后手写 NIPALS ──
  y <- as.numeric(meta$group == CASE)
  Xs <- scale(t(X))
  fit <- nipals_pls1(Xs, y - mean(y), ncomp=2)
  acc <- plsda_loo_acc(t(X), y, ncomp=2)
  cat(sprintf("   PLS-DA LOO 准确率: %.1f%%\n", acc))
  vip <- pls_vip(fit); names(vip) <- rownames(X)
  vdf <- data.frame(gene=rownames(X), VIP=vip)
  vdf <- vdf[order(-vdf$VIP), ]
  save_table(vdf, "02_PLSDA_VIP", row.names=FALSE)
  cat(sprintf("   VIP>1 基因数: %d\n", sum(vip > 1)))

  sc2 <- data.frame(LV1=fit$scores[,1], LV2=fit$scores[,2],
                    group=meta$group, cohort=meta$cohort)
  pP <- ggplot(sc2, aes(LV1, LV2, color=group)) +
    geom_point(size=1.8, alpha=0.8) +
    scale_color_manual(values=COL_GRP, labels=GROUP_LABEL_EN, name="Group") +
    geom_vline(xintercept=0, linetype="dashed", linewidth=LW, color="grey60") +
    geom_hline(yintercept=0, linetype="dashed", linewidth=LW, color="grey60") +
    labs(x="PLS component 1", y="PLS component 2",
         title=sprintf("PLS-DA scores (LOO accuracy %.1f%%)", acc)) + theme_pub_bw()
  savefig(pP, "02_PLSDA_scores", w=4.6, h=3.6)

  cat("=== 02 完成 ===\n")
}

if (sys.nframe() == 0) main()
