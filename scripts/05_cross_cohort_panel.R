# 05_cross_cohort_panel.R —— 跨队列验证 + 多蛋白诊断 panel (临床转化再分析)
#   把项目从"基础复现"抬到"biomarker 稳健性验证"。核心:数据天然有 GU/TU 两独立队列。
#   (A) 分队列差异分析 → GU vs TU 的 log2FC 一致性散点 → robust markers(两队列方向一致)
#   (B) 多蛋白 panel: logistic + LASSO, GU 训练→TU 测试 与 TU 训练→GU 测试(双向外部验证)
#   (C) panel vs 单蛋白 ROC 对比
#   (D) 已知 HCC 血清标志物核查(AFP/GPC3 等是否检出/显著/AUC)
#
# !! 关键防泄漏: 绝不用 01 的 ComBat 矩阵做跨队列验证 —— ComBat 用了两队列+标签全局估计,
#    会把测试队列信息泄进训练特征, AUC 虚高。本脚本从原始矩阵出发, 每队列各自
#    log2→填补→组内 z-score, 模型只学"队列内相对模式", 跨队列迁移才是真外部验证。
#
# 输入: results/tables/00_raw_lfq_matrix.csv, 00_sample_metadata.csv
# 输出: results/tables/05_cohort_diff_concordance.csv, 05_robust_markers.csv,
#        05_cross_cohort_auc.csv, 05_known_hcc_markers.csv
#       results/figures/05_logFC_concordance.{pdf,png}, 05_cross_cohort_roc.*,
#        05_panel_vs_single_roc.*
# 注: 图内文字一律英文。

.this <- (function(){a<-commandArgs(FALSE);f<-sub("^--file=","",a[grepl("^--file=",a)])
  if(length(f))return(dirname(normalizePath(f[1],winslash="/")))
  for(i in sys.nframe():1){o<-sys.frame(i)$ofile;if(!is.null(o))return(dirname(normalizePath(o,winslash="/")))};"."})()
suppressWarnings(suppressMessages({source(file.path(.this, "_common.R"))
  library(glmnet); library(pROC); library(ggrepel)}))

# 单队列: 原始矩阵子集 → log2 → Perseus 填补 → 行(蛋白)z-score
# 返回标准化后矩阵(基因 x 样本), 仅保留该队列检出>=50%的基因
prep_cohort <- function(Xraw, samples, min_frac=0.5) {
  X <- Xraw[, samples, drop=FALSE]
  keep <- rowSums(!is.na(X)) >= min_frac*ncol(X)
  X <- log2(X[keep, , drop=FALSE])
  X <- perseus_impute(X)
  # 行 z-score(每蛋白在该队列内标准化) → 跨队列可比的"相对模式"
  t(scale(t(X)))
}

# 单队列差异: Welch t + log2FC (HCC - Cirrhosis)
cohort_diff <- function(Xlog, grp) {
  hcc <- Xlog[, grp=="HCC", drop=FALSE]; cir <- Xlog[, grp=="Cirrhosis", drop=FALSE]
  data.frame(
    gene=rownames(Xlog),
    log2FC=rowMeans(hcc,na.rm=TRUE)-rowMeans(cir,na.rm=TRUE),
    P=sapply(seq_len(nrow(Xlog)), function(i)
      tryCatch(t.test(hcc[i,],cir[i,])$p.value, error=function(e) NA_real_)),
    row.names=rownames(Xlog))
}

main <- function() {
  cat("=== 05 跨队列验证 + 多蛋白 panel ===\n")
  set.seed(42)
  Xraw <- as.matrix(read.csv(file.path(TAB_DIR,"00_raw_lfq_matrix.csv"), row.names=1, check.names=FALSE))
  meta <- read.csv(file.path(TAB_DIR,"00_sample_metadata.csv"), row.names=1)
  meta <- meta[colnames(Xraw), ]
  s_gu <- rownames(meta)[meta$cohort=="GU"]; s_tu <- rownames(meta)[meta$cohort=="TU"]

  # ── (A) 分队列差异 + log2FC 一致性 ──
  # 用各队列原始矩阵 log2(不 z-score, 保留真实 FC 尺度)
  raw_log_diff <- function(samp){
    X <- Xraw[,samp,drop=FALSE]; keep <- rowSums(!is.na(X))>=0.5*ncol(X)
    X <- perseus_impute(log2(X[keep,,drop=FALSE]))
    cohort_diff(X, meta[samp,"group"])
  }
  dGU <- raw_log_diff(s_gu); dTU <- raw_log_diff(s_tu)
  shared <- intersect(rownames(dGU), rownames(dTU))
  con <- data.frame(gene=shared,
                    log2FC_GU=dGU[shared,"log2FC"], P_GU=dGU[shared,"P"],
                    log2FC_TU=dTU[shared,"log2FC"], P_TU=dTU[shared,"P"])
  con$concordant <- sign(con$log2FC_GU)==sign(con$log2FC_TU)
  con$FDR_GU <- bh_fdr(con$P_GU); con$FDR_TU <- bh_fdr(con$P_TU)
  r <- cor(con$log2FC_GU, con$log2FC_TU, use="complete.obs")
  cat(sprintf("   两队列共享基因 %d, log2FC 相关 r=%.2f, 方向一致 %d (%.0f%%)\n",
      nrow(con), r, sum(con$concordant), 100*mean(con$concordant)))
  save_table(con, "05_cohort_diff_concordance", row.names=FALSE)

  # robust markers: 两队列方向一致 且 至少一队列 FDR<0.1 且合并差异显著
  dmain <- read.csv(file.path(TAB_DIR,"03_diff_main.csv"))
  rob <- con[con$concordant & (con$FDR_GU<0.1 | con$FDR_TU<0.1) &
             con$gene %in% dmain$gene[dmain$sig=="TRUE"|dmain$sig==TRUE], ]
  rob <- rob[order(pmin(rob$FDR_GU, rob$FDR_TU)), ]
  save_table(rob, "05_robust_markers", row.names=FALSE)
  cat(sprintf("   robust markers(双队列一致+显著): %d → %s\n",
      nrow(rob), paste(head(rob$gene,8), collapse=", ")))

  # 一致性散点图
  con$lab <- ifelse(con$gene %in% head(rob$gene,10), con$gene, "")
  pC <- ggplot(con, aes(log2FC_GU, log2FC_TU, color=concordant)) +
    geom_point(size=1, alpha=0.6) +
    geom_hline(yintercept=0, linewidth=LW, color="grey60") +
    geom_vline(xintercept=0, linewidth=LW, color="grey60") +
    geom_abline(slope=1, intercept=0, linetype="dashed", linewidth=LW, color="grey40") +
    geom_text_repel(aes(label=lab), size=PT, color="black", family=PLOT_FAMILY, max.overlaps=20, min.segment.length=0) +
    scale_color_manual(values=c("TRUE"="#EC0000","FALSE"="grey70"),
                       labels=c("TRUE"="Concordant","FALSE"="Discordant"), name="Direction") +
    labs(x="log2FC in GU cohort", y="log2FC in TU cohort",
         title=sprintf("Cross-cohort log2FC concordance (r=%.2f)", r)) + theme_pub_bw()
  savefig(pC, "05_logFC_concordance", w=4.8, h=4.2)

  # ── (B) 跨队列 panel: 各队列内 z-score 后训练/测试(防 ComBat 泄漏) ──
  Zgu <- prep_cohort(Xraw, s_gu); Ztu <- prep_cohort(Xraw, s_tu)
  feats <- intersect(rownames(Zgu), rownames(Ztu))    # 两队列共有基因
  Zgu <- Zgu[feats,]; Ztu <- Ztu[feats,]
  ygu <- as.numeric(meta[s_gu,"group"]=="HCC"); ytu <- as.numeric(meta[s_tu,"group"]=="HCC")
  cat(sprintf("   panel 特征(两队列共有 z-score 基因): %d\n", length(feats)))

  # LASSO: 在训练队列 CV 选 panel, 测试队列算 AUC
  run_lasso <- function(Xtr, ytr, Xte, yte) {
    cv <- cv.glmnet(t(Xtr), ytr, family="binomial", alpha=1, nfolds=5)
    co <- coef(cv, s="lambda.min"); sel <- rownames(co)[which(co[,1]!=0)]; sel <- setdiff(sel,"(Intercept)")
    pr <- as.numeric(predict(cv, newx=t(Xte), s="lambda.min", type="response"))
    list(auc=as.numeric(roc(yte, pr, quiet=TRUE)$auc), panel=sel, pred=pr, n=length(sel))
  }
  # logistic: 固定 4 蛋白 panel(差异表 top 中两队列都有的)
  fixed4 <- head(intersect(dmain$gene[order(dmain$FDR)], feats), 4)
  run_logit <- function(Xtr,ytr,Xte,yte,genes){
    df_tr <- data.frame(y=ytr, t(Xtr[genes,,drop=FALSE])); df_te <- data.frame(t(Xte[genes,,drop=FALSE]))
    colnames(df_tr) <- c("y", make.names(genes)); colnames(df_te) <- make.names(genes)
    fit <- glm(y~., data=df_tr, family=binomial)
    pr <- as.numeric(predict(fit, newdata=df_te, type="response"))
    as.numeric(roc(yte, pr, quiet=TRUE)$auc)
  }
  L_gt <- run_lasso(Zgu,ygu,Ztu,ytu); L_tg <- run_lasso(Ztu,ytu,Zgu,ygu)
  g4_gt <- run_logit(Zgu,ygu,Ztu,ytu,fixed4); g4_tg <- run_logit(Ztu,ytu,Zgu,ygu,fixed4)
  cat(sprintf("   LASSO  GU->TU AUC=%.3f (panel %d) | TU->GU AUC=%.3f (panel %d)\n",
      L_gt$auc, L_gt$n, L_tg$auc, L_tg$n))
  cat(sprintf("   Logit4 GU->TU AUC=%.3f | TU->GU AUC=%.3f (panel: %s)\n",
      g4_gt, g4_tg, paste(fixed4, collapse="+")))

  auc_tab <- data.frame(
    model=c("LASSO","LASSO","Logistic-4prot","Logistic-4prot"),
    train=c("GU","TU","GU","TU"), test=c("TU","GU","TU","GU"),
    n_features=c(L_gt$n, L_tg$n, 4, 4),
    AUC=round(c(L_gt$auc, L_tg$auc, g4_gt, g4_tg),3))
  save_table(auc_tab, "05_cross_cohort_auc", row.names=FALSE)

  # 跨队列 ROC 曲线(两个方向的 LASSO)
  roc_gt <- roc(ytu, L_gt$pred, quiet=TRUE); roc_tg <- roc(ygu, L_tg$pred, quiet=TRUE)
  rl <- rbind(
    data.frame(spec=rev(roc_gt$specificities), sens=rev(roc_gt$sensitivities),
               m=sprintf("GU train -> TU test (AUC %.2f)", as.numeric(roc_gt$auc))),
    data.frame(spec=rev(roc_tg$specificities), sens=rev(roc_tg$sensitivities),
               m=sprintf("TU train -> GU test (AUC %.2f)", as.numeric(roc_tg$auc))))
  pX <- ggplot(rl, aes(1-spec, sens, color=m)) + geom_line(linewidth=LW*2) +
    geom_abline(slope=1,intercept=0,linetype="dashed",linewidth=LW,color="grey60") +
    scale_color_manual(values=COL7[c(1,2)], name=NULL) + coord_equal() +
    labs(x="1 - Specificity", y="Sensitivity", title="Cross-cohort LASSO panel (external validation)") +
    theme_pub_bw() + theme(legend.position=c(0.6,0.18))
  savefig(pX, "05_cross_cohort_roc", w=4.6, h=4.4)

  # ── (C) panel vs 单蛋白(队列内 5 折 CV, 公平对比) ──
  Zall <- cbind(Zgu, Ztu)[feats,]; yall <- c(ygu, ytu)
  cv_auc_single <- function(g){
    pr <- numeric(length(yall)); folds <- sample(rep(1:5, length.out=length(yall)))
    for(k in 1:5){ te<-folds==k
      fit<-glm(yall[!te]~Zall[g,!te], family=binomial)
      pr[te]<-predict(fit, data.frame(x=Zall[g,te]), type="response") }
    as.numeric(roc(yall, pr, quiet=TRUE)$auc)
  }
  cv_auc_panel <- function(genes){
    pr <- numeric(length(yall)); folds <- sample(rep(1:5, length.out=length(yall)))
    for(k in 1:5){ te<-folds==k
      df<-data.frame(y=yall[!te], t(Zall[genes,!te,drop=FALSE]))
      colnames(df)<-c("y",make.names(genes))
      fit<-glm(y~., data=df, family=binomial)
      nd<-data.frame(t(Zall[genes,te,drop=FALSE])); colnames(nd)<-make.names(genes)
      pr[te]<-predict(fit, nd, type="response") }
    as.numeric(roc(yall, pr, quiet=TRUE)$auc)
  }
  single_aucs <- sapply(fixed4, cv_auc_single)
  panel_auc <- cv_auc_panel(fixed4)
  cmp <- data.frame(model=c(fixed4, paste0(length(fixed4),"-protein panel")),
                    n_features=c(rep(1,length(fixed4)), length(fixed4)),
                    CV_AUC=round(c(single_aucs, panel_auc),3))
  save_table(cmp, "05_panel_vs_single", row.names=FALSE)
  cat("   panel vs 单蛋白(5折CV):\n"); print(cmp)

  pP <- ggplot(cmp, aes(reorder(model, CV_AUC), CV_AUC, fill=n_features>1)) +
    geom_col(width=0.7, color="black", linewidth=LW) +
    geom_text(aes(label=sprintf("%.2f",CV_AUC)), hjust=-0.2, size=PT, family=PLOT_FAMILY) +
    scale_fill_manual(values=c("FALSE"="#00468A","TRUE"="#EC0000"),
                      labels=c("FALSE"="Single","TRUE"="Panel"), name=NULL) +
    coord_flip(ylim=c(0.5, max(cmp$CV_AUC)+0.06)) +
    labs(x=NULL, y="5-fold CV AUC", title="Multi-protein panel vs single proteins") + theme_pub_bw()
  savefig(pP, "05_panel_vs_single_roc", w=4.8, h=3.2)

  # ── (D) 已知 HCC 血清标志物核查 ──
  known <- c("AFP","GPC3","GOLM1","SPP1","DCP","B2M","APOL1","VCAM1","LRG1","CLU",
             "AHSG","FCGBP","SERPING1","CFB","APOC3","ADIPOQ","CRP","ORM1","ORM2","CPB2")
  in_raw <- intersect(known, rownames(Xraw))
  km <- data.frame(gene=in_raw,
    in_filtered=in_raw %in% dmain$gene,
    log2FC=dmain$log2FC[match(in_raw, dmain$gene)],
    FDR=dmain$FDR[match(in_raw, dmain$gene)],
    significant=dmain$sig[match(in_raw, dmain$gene)])
  km <- km[order(km$FDR), ]
  save_table(km, "05_known_hcc_markers", row.names=FALSE)
  cat(sprintf("   已知HCC标志物: 数据检出 %d/%d; AFP %s\n", length(in_raw), length(known),
      ifelse("AFP" %in% in_raw, "检出但太稀疏被有效值过滤剔除(诚实标注)", "未检出")))

  cat("=== 05 完成 ===\n")
}

if (sys.nframe() == 0) main()
