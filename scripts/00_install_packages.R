# 00_install_packages.R —— 依赖安装 (仅首次)
# 运行: Rscript scripts/00_install_packages.R

cran <- c("data.table", "matrixStats", "ggplot2", "ggrepel", "pheatmap", "RColorBrewer", "pROC")
# 人类是模式物种: 有真 GO/KEGG 注释库, 可做真富集 (对照 04 非模式物种做不了)
bioc <- c("limma", "clusterProfiler", "org.Hs.eg.db", "sva")  # sva: 两队列批次校正(ComBat)

for (p in cran) if (!requireNamespace(p, quietly=TRUE)) install.packages(p, repos="https://cloud.r-project.org")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")
for (p in bioc) if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, update=FALSE, ask=FALSE)

cat("依赖检查完成。\n")
