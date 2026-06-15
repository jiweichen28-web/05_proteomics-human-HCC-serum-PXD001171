# 06_gsea_hallmark.R —— 通路富集进阶: MSigDB Hallmark ORA + GSEA preranked
#   补 04 的 GO/KEGG ORA 之不足: (1) ORA 阈值依赖强, GSEA 用全部基因排序信息更稳;
#   (2) Hallmark 对人类疾病解释更友好(凝血/补体/EMT/炎症等)。
#   注: 原计划也做 Reactome(ReactomePA), 但 reactome.db(454MB) 本网络下载失败,
#       脚本保留 Reactome 分支(若日后装上自动生成), 当前以 Hallmark 为主。
#
# 输入: results/tables/03_diff_main.csv
# 输出: results/tables/06_hallmark_ora.csv, 06_gsea_hallmark.csv (+若装了ReactomePA: 06_reactome_*.csv)
#       results/figures/06_hallmark_bar.{pdf,png}, 06_gsea_hallmark_dot.*
# 注: 血清深度耗竭后仅 161 基因, GSEA 统计偏薄, 结果作补充并诚实标注; 图内英文。

.this <- (function(){a<-commandArgs(FALSE);f<-sub("^--file=","",a[grepl("^--file=",a)])
  if(length(f))return(dirname(normalizePath(f[1],winslash="/")))
  for(i in sys.nframe():1){o<-sys.frame(i)$ofile;if(!is.null(o))return(dirname(normalizePath(o,winslash="/")))};"."})()
suppressWarnings(suppressMessages({source(file.path(.this, "_common.R"))
  library(clusterProfiler); library(org.Hs.eg.db)}))

# 通用条形图: 通路名(截断换行) x -log10(p.adjust)
bar_enrich <- function(df, title, name, n=12) {
  df <- head(df[order(df$p.adjust), ], n)
  desc <- ifelse(nchar(df$Description) > 50, paste0(substr(df$Description,1,47),"..."), df$Description)
  df$lab <- factor(stringr::str_wrap(desc, 34), levels=rev(stringr::str_wrap(desc, 34)))
  p <- ggplot(df, aes(-log10(p.adjust), lab, fill=Count)) +
    geom_col(width=0.7, color="black", linewidth=LW) +
    scale_fill_gradient(low="#FCAE91", high="#AC002A", name="Count") +
    labs(x="-log10 adj.P", y=NULL, title=title) + theme_pub_bw()
  savefig(p, name, w=7.2, h=4.2)
}

main <- function() {
  cat("=== 06 GSEA + Reactome/Hallmark ===\n")
  out <- read.csv(file.path(TAB_DIR, "03_diff_main.csv"))
  # SYMBOL -> ENTREZ
  map <- bitr(out$gene, "SYMBOL", "ENTREZID", org.Hs.eg.db)
  out$entrez <- map$ENTREZID[match(out$gene, map$SYMBOL)]
  sig_eg <- out$entrez[(out$sig==TRUE | out$sig=="TRUE") & !is.na(out$entrez)]
  cat(sprintf("   显著基因映射到 ENTREZ: %d\n", length(sig_eg)))

  # preranked 统计量: signed -log10(P) = sign(log2FC) * -log10(P)
  rk <- out[!is.na(out$entrez), ]
  stat <- sign(rk$log2FC) * -log10(rk$P)
  names(stat) <- rk$entrez
  stat <- sort(stat[!duplicated(names(stat))], decreasing=TRUE)
  cat(sprintf("   GSEA 排序基因: %d (signed -log10P)\n", length(stat)))

  # ── (A) Hallmark ORA (msigdbr) ──
  if (requireNamespace("msigdbr", quietly=TRUE)) {
    library(msigdbr)
    hm <- msigdbr(species="Homo sapiens", category="H")
    t2g <- hm[, c("gs_name","entrez_gene")]; t2g$gs_name <- sub("^HALLMARK_","",t2g$gs_name)
    eh <- enricher(sig_eg, TERM2GENE=as.data.frame(t2g), pvalueCutoff=0.1, qvalueCutoff=0.25)
    if (!is.null(eh) && nrow(as.data.frame(eh))>0) {
      hd <- as.data.frame(eh); save_table(hd, "06_hallmark_ora", row.names=FALSE)
      bar_enrich(hd, "MSigDB Hallmark ORA (significant proteins)", "06_hallmark_bar")
      cat(sprintf("   Hallmark ORA: %d 条\n", nrow(hd)))
    } else cat("   Hallmark ORA 无显著条目\n")

    # GSEA on Hallmark
    gh <- tryCatch(GSEA(stat, TERM2GENE=as.data.frame(t2g), pvalueCutoff=0.25, eps=0),
                   error=function(e) NULL)
    if (!is.null(gh) && nrow(as.data.frame(gh))>0) {
      gd <- as.data.frame(gh); save_table(gd[,c("ID","setSize","NES","pvalue","p.adjust")], "06_gsea_hallmark", row.names=FALSE)
      cat(sprintf("   Hallmark GSEA: %d 条 (top NES: %s)\n", nrow(gd),
          paste(head(gd$ID[order(-abs(gd$NES))],3), collapse=", ")))
      # Hallmark GSEA 点图: NES 正=HCC上调通路, 负=下调
      gd$lab <- factor(gd$ID, levels=gd$ID[order(gd$NES)])
      pgh <- ggplot(gd, aes(NES, lab, color=p.adjust, size=setSize)) + geom_point() +
        geom_vline(xintercept=0, linewidth=LW, color="grey60") +
        scale_color_gradient(low="#EC0000", high="#00468A", name="p.adjust") +
        scale_size_continuous(range=c(3,6), name="Set size") +
        labs(x="NES (HCC vs Cirrhosis)", y=NULL, title="MSigDB Hallmark GSEA") + theme_pub_bw()
      savefig(pgh, "06_gsea_hallmark_dot", w=6.2, h=3.2)
    } else cat("   Hallmark GSEA 无显著条目(161基因偏薄, 正常)\n")
  } else cat("   msigdbr 未安装, 跳过 Hallmark\n")

  # ── (B) Reactome ORA + GSEA (ReactomePA) ──
  if (requireNamespace("ReactomePA", quietly=TRUE)) {
    library(ReactomePA)
    er <- tryCatch(enrichPathway(sig_eg, organism="human", pvalueCutoff=0.1, readable=TRUE),
                   error=function(e) NULL)
    if (!is.null(er) && nrow(as.data.frame(er))>0) {
      rd <- as.data.frame(er); save_table(rd, "06_reactome_ora", row.names=FALSE)
      bar_enrich(rd, "Reactome pathway ORA (significant proteins)", "06_reactome_bar")
      cat(sprintf("   Reactome ORA: %d 条\n", nrow(rd)))
    } else cat("   Reactome ORA 无显著条目\n")

    gr <- tryCatch(gsePathway(stat, organism="human", pvalueCutoff=0.25, eps=0),
                   error=function(e) NULL)
    if (!is.null(gr) && nrow(as.data.frame(gr))>0) {
      grd <- as.data.frame(gr)
      save_table(grd[,c("ID","Description","setSize","NES","pvalue","p.adjust")], "06_gsea_reactome", row.names=FALSE)
      # GSEA 点图
      gt <- head(grd[order(grd$p.adjust),], 12)
      desc <- ifelse(nchar(gt$Description)>50, paste0(substr(gt$Description,1,47),"..."), gt$Description)
      gt$lab <- factor(stringr::str_wrap(desc,34), levels=rev(stringr::str_wrap(desc,34)))
      pg <- ggplot(gt, aes(NES, lab, color=p.adjust, size=setSize)) + geom_point() +
        geom_vline(xintercept=0, linewidth=LW, color="grey60") +
        scale_color_gradient(low="#EC0000", high="#00468A", name="p.adjust") +
        scale_size_continuous(range=c(2,5), name="Set size") +
        labs(x="NES (HCC vs Cirrhosis)", y=NULL, title="Reactome GSEA") + theme_pub_bw()
      savefig(pg, "06_gsea_dot", w=7.2, h=4.2)
      cat(sprintf("   Reactome GSEA: %d 条\n", nrow(grd)))
    } else cat("   Reactome GSEA 无显著条目\n")
  } else cat("   ReactomePA 未安装, 跳过 Reactome\n")

  cat("=== 06 完成 ===\n")
}

if (sys.nframe() == 0) main()
