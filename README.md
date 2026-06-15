# 人类血清蛋白质组：肝细胞癌 vs 肝硬化（HCC vs Cirrhosis）

LC-MS/MS 鸟枪蛋白质组，比较肝细胞癌（HCC）与肝硬化患者血清的蛋白质组差异，寻找诊断标志物。数据来自 ProteomeXchange PXD001171（MaxQuant 标记自由定量），含两个独立队列（美国 Georgetown + 埃及 Tanta）。本项目为其再分析。

标签：`proteomics` `LC-MS/MS` `MaxQuant` `LFQ` `human` `HCC` `biomarker` `ROC` `GO-enrichment` `batch-correction` `R`

## 数据来源

| 字段 | 内容 |
|------|------|
| 数据库 | ProteomeXchange / PRIDE [PXD001171](https://www.ebi.ac.uk/pride/archive/projects/PXD001171) |
| 物种 | *Homo sapiens*（NCBI TaxID 9606）|
| 平台 | Thermo LTQ Orbitrap Velos；MaxQuant 检索，LFQ 定量 |
| 分组 | HCC（肝细胞癌）vs Cirrhosis（肝硬化，对照）|
| 队列 | GU = Georgetown（美国，57 HCC / 59 cirrhosis）；TU = Tanta（埃及，40 HCC / 49 cirrhosis）；共 205 |
| 样本 | 血清，Agilent Plasma 7 去高丰度，胰酶消化 |
| 原文 | Tsai TH 等, *Proteomics* 2015, 15(13):2369-81. DOI [10.1002/pmic.201400364](https://doi.org/10.1002/pmic.201400364) · PMID [25778709](https://pubmed.ncbi.nlm.nih.gov/25778709) |

> 数据由用户自 PRIDE 下载 `GU_MaxQuant_txt.zip` + `TU_MaxQuant_txt.zip`，解压后将各自的 `proteinGroups.txt` 置于 `data/GU/` 与 `data/TU/`。原始谱图 `.mzXML`（98 个）不需要、不入库。

## 目录结构

```
05_proteomics-human-HCC-serum-PXD001171/
├── data/                       # MaxQuant 文本结果（大文件 .mzXML/.zip 不入库）
│   ├── GU/proteinGroups.txt
│   └── TU/proteinGroups.txt
├── scripts/                    # _common.R + 00_install + 00–04
├── results/
│   ├── figures/                # 每图 .pdf + .png(300ppi)，图内文字一律英文
│   └── tables/                 # 原始/清洗矩阵、差异表、VIP、富集、ROC
├── README.md  LICENSE  .gitignore  PROJECT_STATUS.md  sessionInfo.txt
```

## 分析脚本

```
00_install_packages.R  依赖安装（仅首次）
_common.R              路径自定位 + 绘图风格 + 两队列 proteinGroups 加载 + Perseus/统计工具
00_setup_and_load.R    读两队列 proteinGroups → 去污染/反库/仅位点 → 抽 LFQ → 0→NA → 合并 → 关联分组/队列
01_clean_normalize.R   有效值过滤 → log2 → ComBat 批次校正 → Perseus 填补 → QC 图
02_pca_plsda.R         PCA（疾病/队列双上色）+ PLS-DA（LOO + VIP）
03_diff_proteins.R     HCC vs 肝硬化差异（limma + 队列校正）→ 火山图/热图/箱线
04_enrichment_roc.R    GO/KEGG 富集（clusterProfiler）+ top 蛋白 ROC/AUC 诊断性能
```

运行（从项目任意目录，脚本自动定位项目根，按编号逐个跑）：

```bash
Rscript scripts/00_install_packages.R   # 仅首次，装依赖
Rscript scripts/00_setup_and_load.R     # 读两队列 → 合并 → LFQ 矩阵
Rscript scripts/01_clean_normalize.R    # 过滤 → log2 → ComBat → 填补 → QC 图
Rscript scripts/02_pca_plsda.R          # PCA + PLS-DA（LOO + VIP）
Rscript scripts/03_diff_proteins.R      # 差异蛋白 → 火山/热图/箱线
Rscript scripts/04_enrichment_roc.R     # GO/KEGG 富集 + ROC
```

`_common.R` 自动定位项目根，从任意目录运行均可。PLS-DA 为手写 NIPALS 实现。

## 主要结果

**鉴定与合并（两队列）**

- 两队列各自 MaxQuant `proteinGroups.txt`：GU 274 → 244 蛋白组、TU 290 → 261（去 Reverse/Contaminant/仅位点）。血清经高丰度蛋白去除后鉴定数本就不多。
- 按基因名取两队列**交集 212 个**，再按「至少一组内 ≥50% 检出」过滤至 **161 个基因 × 205 样本**（HCC 97 / Cirrhosis 108；GU 116 / TU 89）。合并矩阵缺失约 34%，Perseus 下移填补（mean−1.8·sd, 0.3·sd，`seed=42`）。
- 旧版 MaxQuant 无 `Gene names` 列，基因名从 `Fasta headers` 的 `GN=` 正则抽取（GU 230/244 命中）——人是模式物种，蛋白可映射到基因，这正是与项目 04（非模式硅藻只能用 `Phatr3_J` 编号）的关键差异。

**批次效应（两中心 → ComBat）**

- 两队列是不同国家/仪器批次，PCA 上**校正前 PC1 完全按队列分开**（技术批次主导变异）；ComBat 按队列校正、保护疾病分组后两队列混匀（见 `01_QC_batch_pca`）。这是人多中心数据特有的一步，项目 04 等基因型培养无此问题。

**总体结构（PCA / PLS-DA）**

- PCA：PC1 仅 **22.3%** / PC2 7.7% / PC3 5.0%。对比项目 04 硅藻 PC1=59.4% 完全分开两组——**人异质性远高，PC1 不再是干净的疾病轴**，疾病信号被个体差异稀释。
- PLS-DA（手写 NIPALS）LOO 准确率 **66.3%**，VIP>1 共 **72 个**。远低于 04 的 100%——真实人群血清诊断本就难，这个数字是诚实的。

**差异蛋白（HCC vs 肝硬化，limma `~ group + cohort` 队列校正）**

- 阈值 FDR<0.05 且 |log2FC|≥log2(1.5)：共 **19 个显著（上调 4 / 下调 15）**。
- 最显著下调：`B2M`（β2-微球蛋白，FDR 6.1e-5）、`VCAM1`、`EFEMP1`、`LRG1`、`FCGBP`；上调：`APOL1`、`CLU`、`CFB`、`APOC3`。均为可解释的肝病/炎症/急性期血清蛋白。

**功能富集（人是模式物种 → 真 GO/KEGG，04 做不到）**

- GO BP 富集 **111 条**（以免疫球蛋白介导免疫、补体/体液免疫、B 细胞免疫为主）；KEGG **5 条**（含胆固醇代谢 hsa04979、PPAR 信号 hsa03320 等）。用全基因组背景。
- 对照项目 04：硅藻无基因注释，被迫放弃富集改做 on/off 开关蛋白专题。**这是"人 vs 非模式物种"最大的下游差异**。

**诊断性能（ROC/AUC）**

- top 单蛋白对 HCC vs 肝硬化的判别：`B2M` AUC **0.691**、`APOL1` 0.674、`EFEMP1` 0.667、`VCAM1` 0.665。中等——单蛋白血清诊断本就有限，临床通常需多标志物 panel 组合。

> **诚实说明**：本再分析仅用公开的 `proteinGroups.txt`，未做原文的 MRM 靶向验证；血清深度耗竭后蛋白数少（161），PLS-DA/AUC 偏中等是数据本身决定的，不是流程问题。结果方向（HCC 血清免疫球蛋白/急性期蛋白重塑）与原文一致。

## 图件输出

`results/figures/`，每图含 `.pdf` + `.png`(300ppi)，图内文字一律英文：

| 编号 | 内容 |
|------|------|
| 01 | 各样本检出蛋白数、log2 强度箱线、填补密度图、批次校正前后队列分离 |
| 02 | PCA 得分（疾病上色 / 队列上色两版）、Scree、PLS-DA 得分 |
| 03 | 火山图、显著蛋白热图、top 蛋白分组箱线 |
| 04 | GO 富集气泡图、KEGG 富集条形、top 标志物 ROC 曲线 |

## 环境

```
R 4.x：data.table, matrixStats, ggplot2, ggrepel, pheatmap, RColorBrewer, pROC
Bioconductor：limma（差异）、clusterProfiler + org.Hs.eg.db（GO/KEGG）、sva（ComBat 批次校正）
```

图内文字一律英文，`theme_pub_bw()` 用 sans 字体，无需中文字体。`sessionInfo.txt` 记录完整版本。

## 引用

> Tsai TH, Song E, Zhu R, Di Poto C, Wang M, Luo Y, Varghese RS, Tadesse MG, Ziada DH, Desai CS, Shetty K, Mechref Y, Ressom HW. LC-MS/MS-based serum proteomics for identification of candidate biomarkers for hepatocellular carcinoma. *Proteomics* 2015;15(13):2369-2381. https://doi.org/10.1002/pmic.201400364
