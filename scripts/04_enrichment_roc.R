# 04_enrichment_roc.R —— 功能富集 + 诊断性能 (人类是模式物种, 这是相对 04 最大的不同)
#   (A) GO/KEGG 富集: 显著蛋白 → clusterProfiler(org.Hs.eg.db) → 气泡图/条形图
#   (B) ROC/AUC: top 差异蛋白对 HCC vs 肝硬化的诊断性能 → pROC 曲线
# 输入: results/tables/03_diff_main.csv, 01_clean_log2_imputed.csv, 00_sample_metadata.csv
# 输出: results/tables/04_GO_enrichment.csv, 04_KEGG_enrichment.csv, 04_roc_auc.csv
#       results/figures/04_GO_dotplot.{pdf,png}, 04_KEGG_barplot.*, 04_roc_curve.*
# 注: 人蛋白映射到基因 → 真 GO/KEGG(对照 04 非模式物种被迫砍掉假富集改 on/off 专题)。

.this <- (function(){a<-commandArgs(FALSE);f<-sub("^--file=","",a[grepl("^--file=",a)])
  if(length(f))return(dirname(normalizePath(f[1],winslash="/")))
  for(i in sys.nframe():1){o<-sys.frame(i)$ofile;if(!is.null(o))return(dirname(normalizePath(o,winslash="/")))};"."})()
suppressWarnings(suppressMessages({source(file.path(.this, "_common.R"))
  library(clusterProfiler); library(org.Hs.eg.db); library(pROC)}))

main <- function() {
  cat("=== 04 富集 + ROC ===\n")
  out  <- read.csv(file.path(TAB_DIR, "03_diff_main.csv"))
  X    <- as.matrix(read.csv(file.path(TAB_DIR, "01_clean_log2_imputed.csv"), row.names=1, check.names=FALSE))
  meta <- read.csv(file.path(TAB_DIR, "00_sample_metadata.csv"), row.names=1)
  meta <- meta[colnames(X), ]

  # 背景=本数据所有定量基因; 前景=显著差异基因
  bg   <- out$gene
  sigg <- out$gene[out$sig == TRUE | out$sig == "TRUE"]
  cat(sprintf("   背景基因 %d, 显著基因 %d\n", length(bg), length(sigg)))
  bg_eg  <- bitr(bg,  "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID
  sig_eg <- bitr(sigg,"SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID

  # ── (A) GO 富集 (BP) ──
  # 注: 血清深度耗竭后蛋白少(161), 若再限制 universe 到这 161 个则统计功效过低;
  # 改用全基因组背景(标准做法), 让真实肝病/补体/急性期通路浮现(对照 04 非模式物种根本做不了)
  ego <- enrichGO(sig_eg, OrgDb=org.Hs.eg.db, keyType="ENTREZID",
                  ont="BP", pvalueCutoff=0.05, qvalueCutoff=0.1, readable=TRUE)
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    god <- as.data.frame(ego); save_table(god, "04_GO_enrichment", row.names=FALSE)
    topgo <- head(god[order(god$p.adjust), ], 12)
    # 长 GO 名: 先截断超长(>55字符)再换行, 避免某条目占太多行与邻行重叠
    desc <- ifelse(nchar(topgo$Description) > 55,
                   paste0(substr(topgo$Description, 1, 52), "..."), topgo$Description)
    topgo$lab <- factor(stringr::str_wrap(desc, 30),
                        levels=rev(stringr::str_wrap(desc, 30)))
    pG <- ggplot(topgo, aes(Count, lab)) +
      geom_point(aes(color=p.adjust, size=Count)) +
      scale_color_gradient(low="#EC0000", high="#00468A", name="p.adjust") +
      scale_size_continuous(range=c(2,5), name="Count") +
      scale_x_continuous(breaks=scales::pretty_breaks()) +
      labs(x="Gene count", y=NULL, title="GO BP enrichment (significant proteins)") +
      theme_pub_bw()
    savefig(pG, "04_GO_dotplot", w=7.5, h=5)
    cat(sprintf("   GO BP 富集: %d 条目\n", nrow(god)))
  } else cat("   GO 富集无显著条目(血清蛋白背景小, 正常)\n")

  # ── KEGG 富集 (全基因组背景) ──
  ek <- tryCatch(enrichKEGG(sig_eg, organism="hsa", pvalueCutoff=0.05),
                 error=function(e) NULL)
  if (!is.null(ek) && nrow(as.data.frame(ek)) > 0) {
    kd <- as.data.frame(ek); save_table(kd, "04_KEGG_enrichment", row.names=FALSE)
    topk <- head(kd[order(kd$p.adjust), ], 12)
    topk$Description <- factor(topk$Description, levels=rev(topk$Description))
    pK <- ggplot(topk, aes(Count, Description, fill=p.adjust)) +
      geom_col(width=0.7, color="black", linewidth=LW) +
      scale_fill_gradient(low="#EC0000", high="#00468A", name="p.adjust") +
      labs(x="Gene count", y=NULL, title="KEGG enrichment (significant proteins)") +
      theme_pub_bw()
    savefig(pK, "04_KEGG_barplot", w=6.5, h=4)
    cat(sprintf("   KEGG 富集: %d 条目\n", nrow(kd)))
  } else cat("   KEGG 富集无显著条目\n")

  # ── (B) ROC/AUC: top 差异蛋白单标志物诊断 HCC ──
  y <- as.numeric(meta$group == CASE)
  topg <- head(out$gene[out$sig == TRUE | out$sig == "TRUE"], 4)
  rocs <- lapply(topg, function(g) roc(y, X[g,], quiet=TRUE, direction="auto"))
  names(rocs) <- topg
  aucs <- sapply(rocs, function(r) as.numeric(r$auc))
  rdf <- data.frame(gene=topg, AUC=round(aucs,3))
  save_table(rdf, "04_roc_auc", row.names=FALSE)
  cat("   top 标志物 AUC:\n"); print(rdf)

  # ROC 曲线叠加
  roc_long <- do.call(rbind, lapply(topg, function(g) {
    r <- rocs[[g]]
    data.frame(gene=sprintf("%s (AUC %.2f)", g, as.numeric(r$auc)),
               spec=rev(r$specificities), sens=rev(r$sensitivities))
  }))
  pR <- ggplot(roc_long, aes(1-spec, sens, color=gene)) +
    geom_line(linewidth=LW*2) +
    geom_abline(slope=1, intercept=0, linetype="dashed", linewidth=LW, color="grey60") +
    scale_color_manual(values=COL7[1:length(topg)], name=NULL) +
    coord_equal() +
    labs(x="1 - Specificity", y="Sensitivity",
         title="Single-protein ROC (HCC vs Cirrhosis)") +
    theme_pub_bw() + theme(legend.position=c(0.65,0.25))
  savefig(pR, "04_roc_curve", w=4.4, h=4.2)

  cat("=== 04 完成 ===\n")
}

if (sys.nframe() == 0) main()
