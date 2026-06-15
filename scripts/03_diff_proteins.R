# 03_diff_proteins.R —— 差异蛋白: HCC vs 肝硬化。limma 线性模型(~ group + cohort 校正批次) +
#   log2FC + BH-FDR, 合并 VIP → 火山图/显著蛋白热图/top 蛋白箱线
# 输入: 01_clean_log2_imputed.csv, 00_sample_metadata.csv, 02_PLSDA_VIP.csv
# 输出: results/tables/03_diff_main.csv
#       results/figures/03_volcano.{pdf,png}, 03_heatmap.*, 03_top_boxplot.*
# 注: 人异质队列用 limma + 协变量校正(对照 04 两组干净直接 Welch t)。

.this <- (function(){a<-commandArgs(FALSE);f<-sub("^--file=","",a[grepl("^--file=",a)])
  if(length(f))return(dirname(normalizePath(f[1],winslash="/")))
  for(i in sys.nframe():1){o<-sys.frame(i)$ofile;if(!is.null(o))return(dirname(normalizePath(o,winslash="/")))};"."})()
suppressWarnings(suppressMessages({source(file.path(.this, "_common.R"))
  library(limma); library(ggrepel); library(pheatmap)}))

main <- function() {
  cat("=== 03 差异蛋白 (HCC vs 肝硬化) ===\n")
  X <- as.matrix(read.csv(file.path(TAB_DIR, "01_clean_log2_imputed.csv"), row.names=1, check.names=FALSE))
  meta <- read.csv(file.path(TAB_DIR, "00_sample_metadata.csv"), row.names=1)
  meta <- meta[colnames(X), ]
  vip <- read.csv(file.path(TAB_DIR, "02_PLSDA_VIP.csv"))

  # ── limma: ~ group + cohort, 队列作协变量校正批次 ──
  grp <- factor(meta$group, levels=c(CONTROL, CASE))   # Cirrhosis 为参照
  coh <- factor(meta$cohort)
  design <- model.matrix(~ grp + coh)
  fit <- eBayes(lmFit(X, design))
  res <- topTable(fit, coef="grpHCC", number=Inf, sort.by="none")
  # log2FC = HCC - Cirrhosis (limma logFC 即此)
  out <- data.frame(gene=rownames(res), log2FC=res$logFC, P=res$P.Value,
                    FDR=res$adj.P.Val, row.names=rownames(res))
  out$VIP <- vip$VIP[match(out$gene, vip$gene)]
  out$sig <- out$FDR < FDR_THR & abs(out$log2FC) >= FC_THR
  out$direction <- ifelse(out$log2FC > 0, "up", "down")
  out <- out[order(out$FDR), ]
  save_table(out, "03_diff_main", row.names=FALSE)
  nsig <- sum(out$sig); nup <- sum(out$sig & out$direction=="up")
  cat(sprintf("   显著差异蛋白(FDR<%.2f & |log2FC|>=log2(1.5)): %d (上调%d/下调%d)\n",
      FDR_THR, nsig, nup, nsig-nup))

  # ── 火山图 ──
  out$lab <- ifelse(out$sig, out$gene, "")
  topn <- head(out$gene[out$sig][order(out$FDR[out$sig])], 12)
  out$lab[!(out$gene %in% topn)] <- ""
  pv <- ggplot(out, aes(log2FC, -log10(FDR), color=interaction(sig, direction))) +
    geom_point(size=1, alpha=0.7) +
    scale_color_manual(values=c("TRUE.up"="#EC0000","TRUE.down"="#00468A",
                                "FALSE.up"="grey75","FALSE.down"="grey75"), guide="none") +
    geom_vline(xintercept=c(-FC_THR, FC_THR), linetype="dashed", linewidth=LW, color="grey50") +
    geom_hline(yintercept=-log10(FDR_THR), linetype="dashed", linewidth=LW, color="grey50") +
    geom_text_repel(aes(label=lab), size=PT, color="black", family=PLOT_FAMILY,
                    max.overlaps=20, min.segment.length=0, segment.size=LW) +
    labs(x="log2 fold change (HCC / Cirrhosis)", y="-log10 FDR",
         title="Differential proteins: HCC vs Cirrhosis") + theme_pub_bw()
  savefig(pv, "03_volcano", w=5.2, h=4.2)

  # ── 显著蛋白热图(z-score, 样本按疾病注释) ──
  sigg <- out$gene[out$sig]
  if (length(sigg) >= 2) {
    z <- t(scale(t(X[sigg, , drop=FALSE])))
    ord <- order(meta$group, meta$cohort)
    cann <- data.frame(Group=meta$group, Cohort=meta$cohort, row.names=rownames(meta))
    ann_col <- list(Group=COL_GRP, Cohort=COL_COHORT)
    ph <- pheatmap(z[, ord], cluster_cols=FALSE, show_colnames=FALSE,
                   annotation_col=cann, annotation_colors=ann_col,
                   color=colorRampPalette(c("#00468A","white","#EC0000"))(50),
                   breaks=seq(-2,2,length.out=51), fontsize=FONT_SIZE,
                   main=sprintf("Significant proteins (%d, z-score)", length(sigg)), silent=TRUE)
    pdf(file.path(FIG_DIR,"03_heatmap.pdf"), width=6, height=5); grid::grid.newpage(); grid::grid.draw(ph$gtable); dev.off()
    png(file.path(FIG_DIR,"03_heatmap.png"), width=6, height=5, units="in", res=300); grid::grid.newpage(); grid::grid.draw(ph$gtable); dev.off()
    cat("   [fig] 03_heatmap.pdf / .png\n")
  }

  # ── top 蛋白箱线(按疾病) ──
  topg <- head(out$gene[out$sig][order(out$FDR[out$sig])], 6)
  if (length(topg) >= 1) {
    long <- do.call(rbind, lapply(topg, function(g)
      data.frame(gene=g, value=X[g,], group=meta$group)))
    long$gene <- factor(long$gene, levels=topg)
    pb <- ggplot(long, aes(group, value, fill=group)) +
      geom_boxplot(outlier.size=0.3, linewidth=LW) +
      facet_wrap(~gene, scales="free_y", nrow=2) +
      scale_fill_manual(values=COL_GRP, labels=GROUP_LABEL_EN, name="Group") +
      labs(x=NULL, y="log2 LFQ intensity", title="Top differential proteins") +
      theme_pub_bw() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
    savefig(pb, "03_top_boxplot", w=5.5, h=4)
  }

  cat("=== 03 完成 ===\n")
}

if (sys.nframe() == 0) main()
