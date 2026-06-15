# 00_install_packages.R —— 依赖安装 (仅首次)
# 运行: Rscript scripts/00_install_packages.R

cran <- c("data.table", "matrixStats", "ggplot2", "ggrepel", "pheatmap", "RColorBrewer",
          "pROC", "glmnet", "msigdbr", "stringr", "scales")
# 人类是模式物种: 有真 GO/KEGG 注释库, 可做真富集 (对照 04 非模式物种做不了)
bioc <- c("limma", "clusterProfiler", "org.Hs.eg.db", "sva", "fgsea")  # sva: ComBat 批次校正; fgsea: GSEA
# 可选(深入分析 06): Reactome 通路富集; reactome.db 体积大(454MB), 网络差时可跳过
bioc_optional <- c("ReactomePA")

for (p in cran) if (!requireNamespace(p, quietly=TRUE)) install.packages(p, repos="https://cloud.r-project.org")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")
for (p in bioc) if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, update=FALSE, ask=FALSE)
# 可选包: 失败不阻断(reactome.db 大, 网络差可能下不动; 06 脚本会自动跳过 Reactome)
for (p in bioc_optional) if (!requireNamespace(p, quietly=TRUE))
  tryCatch(BiocManager::install(p, update=FALSE, ask=FALSE),
           error=function(e) cat(sprintf("  [可选] %s 安装失败, 06 脚本将跳过 Reactome 分支\n", p)))

cat("依赖检查完成。\n")
