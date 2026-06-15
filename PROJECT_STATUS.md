# PROJECT_STATUS — 05 人类血清蛋白质组 (HCC vs 肝硬化, PXD001171)

> 本文件追踪 harness 8 阶段进度，供跨会话交接。图内文字一律英文，代码注释中文，README/笔记中文。

## 一句话定位
人类血清 LC-MS/MS 蛋白组（PRIDE PXD001171，MaxQuant LFQ），HCC vs 肝硬化两队列再分析。
**对照项目 04（三角褐指藻氮饥饿）**，演示"人 vs 非模式物种"的下游差异：真基因名、批次校正、真 GO/KEGG、ROC 诊断。

## 流程位置：harness 8 阶段
- P0 INTAKE ✅　P1 SCAFFOLD ✅　P2 DATA ✅　P3 IMPLEMENT ✅　P4 RESULTS ✅　P5 SHIP ✅　P6 LEARN ✅　P7 CLOSE ✅
- **🎉 05 项目全部完成（P0→P7）+ 深入分析扩展（05/06 脚本）**
- GitHub: https://github.com/jiweichen28-web/05_proteomics-human-HCC-serum-PXD001171 (public, commit 3ad46ba)
- 教训笔记: notes-lessons/No.5_lessons_human_proteomics.tex(+pdf, 12页, 不入库)
- 顶层 README 索引已加 05

## 深入分析扩展（基础复现之上，临床转化再分析）
- **05_cross_cohort_panel.R**：跨队列验证 + 多蛋白 panel
  - GU vs TU log2FC 一致性 r=0.73、77% 同向；16 robust markers(APOL1/B2M/FCGBP/VCAM1…)
  - ⚠️ 防 ComBat 标签泄漏：跨队列验证回原始矩阵、每队列组内 z-score(不用 01 的 ComBat 矩阵)
  - 4 蛋白逻辑回归 panel 外部验证 GU→TU AUC 0.73 / TU→GU 0.77；LASSO 0.70/0.68
  - panel vs 单蛋白(5折CV)：单蛋白 ~0.5 vs panel 0.73 → 组合 >> 单一
  - 已知 HCC 标志物核查 16/20 检出；AFP 检出但太稀疏被过滤(诚实标注)
- **06_gsea_hallmark.R**：Hallmark ORA(COAGULATION) + GSEA preranked(COAGULATION+1.53/COMPLEMENT+1.36/EMT−1.51)
  - Reactome 分支保留但本网络 reactome.db(454MB) 下不动 → 当前以 Hallmark 为主(诚实跳过)
- 产出累计：scripts 9 个、19 图(PDF+PNG)、18 表；sessionInfo 已含 glmnet/msigdbr/fgsea

## P3 实测关键结果（填 README 用，已逐脚本验证）
- 合并：212 共享基因 → 有效值过滤(组内≥50%) → **161 基因 × 205 样本**(HCC 97 / Cirrhosis 108; GU 116 / TU 89)
- 缺失：合并矩阵约 34%，Perseus 下移填补；ComBat 按队列校正(vision 核验：校正前两队列 PC1 分开、校正后混匀)
- PCA：PC1 仅 **22.3%** / PC2 7.7% / PC3 5.0%（对照 04 硅藻 PC1=59.4% → 人异质性远高，PC1 不再干净是疾病）
- PLS-DA LOO 准确率 **66.3%**，VIP>1 共 **72**（对照 04 = 100%，人难分得多）
- 差异蛋白(limma ~group+cohort, FDR<0.05 & |log2FC|≥log2(1.5))：**19 个(上调4/下调15)**
  - top: B2M↓(FDR 6e-5), APOL1↑, VCAM1↓, EFEMP1↓, CLU, LRG1, CRP, ADIPOQ —— 全是可解释的肝病/炎症血清蛋白
- **GO BP 富集 111 条**(免疫球蛋白/补体/B细胞免疫为主) + **KEGG 5 条**(胆固醇代谢 hsa04979、PPAR 等)
  → 这是相对 04 最大的不同：人能做真富集；04 非模式物种被迫砍掉假富集
- ROC 单标志物 AUC：B2M 0.691 / APOL1 0.674 / EFEMP1 0.667 / VCAM1 0.665（中等，符合单蛋白血清诊断预期）
- 14 图(PDF+PNG) + 11 表；sessionInfo.txt 已冻结(R 4.4.2)
- 已 vision 核验：批次校正PCA、火山图、GO气泡、ROC —— 全英文无乱码

## P3 踩坑(已解决，记入 No.5 教训)
1. 05 _common.R 初版漏了 NIPALS/VIP/LOO 函数 → 02 报"找不到 nipals_pls1" → 补回
2. 旧版 MaxQuant 无 Gene names 列 → 改从 Fasta headers 正则抽 "GN=" → 真基因名(GU 230/244)
3. 污染物列名 "Contaminant"(非 "Potential contaminant") → FILTER_FLAGS 增补
4. GO 富集限制 universe 到 161 基因 → 统计功效太低无结果 → 改全基因组背景(标准做法)
5. GO 气泡图首版数据点不渲染/超长 GO 名重叠 → 截断>55字符+换行+加宽画布修复(vision 反复核验)

## P2 校验结果（已通过）
- 两队列 proteinGroups.txt 到位（GU 8.25MB / TU 5.54MB），原始大文件+zip 全在但被 .gitignore 白名单挡掉
- 分组映射明确：experimentalDesignTemplate.txt 的 Experiment 列 = HCC_n / CIRR_n
- GU: 274→244 蛋白(过滤后), 116 样本(57 HCC/59 CIRR), 缺失 34.1%
- TU: 290→261 蛋白(过滤后), 89 样本(40 HCC/49 CIRR), 缺失 35.2%
- 两处实测修正(已写入 _common.R)：
  ① 旧版 MaxQuant 污染物列名是 "Contaminant"(非 "Potential contaminant") → FILTER_FLAGS 增补
  ② 无独立 Gene names 列, 但 Fasta headers 含 "GN=" → 正则抽真基因名(GU 230/244, TU 255/261)
     + 抽 UniProt accession 供富集映射。这正是"人 vs 硅藻"教学点：人能拿到真基因名
- 两队列共享基因 212 个 → 合并按基因名取交集可行

## 与 04 的关键差异（本项目要演示的点）
1. **行名真基因名**（Gene names 非空）→ 图可读、可做富集
2. **两队列 = 批次效应** → ComBat 校正；PCA 要按疾病/队列双上色排查
3. **真 GO/KEGG 富集**（org.Hs.eg.db + clusterProfiler）→ 对照 04 被迫砍掉假富集
4. **诊断 biomarker**：top 蛋白 ROC/AUC → 04 没有的临床下游
5. 差异用 **limma + 队列协变量校正** → 04 是干净两组直接 Welch t

## 数据计划（P2）
- 下载：`GU_MaxQuant_txt.zip` + `TU_MaxQuant_txt.zip`（PRIDE FTP，已确认含 proteinGroups.txt）
- 保存：解压后 `proteinGroups.txt` 分别置于 `data/GU/` 与 `data/TU/`
- 不下：98 个 `.mzXML` 原始谱图（不需要、不入库）
- 校验：两队列各自蛋白行数、样本列数、样本→分组(HCC/cirrhosis)映射

## 已建（P1 SCAFFOLD）
- scripts/：_common.R + 00_install_packages.R + 00–04 五个分析脚本（均为 stub，含一行任务说明）
- README.md 骨架（主要结果占位待 P4）
- .gitignore（白名单：仅放行两个 proteinGroups.txt，挡 .mzXML/.zip/大文件）
- LICENSE（MIT，仅代码）
- 本 PROJECT_STATUS.md

## 下一步（P2 DATA-HANDOFF）
给用户直链 + 精确保存路径，用户下载，回来校验行列数与分组映射，再进 P3。

## 数据映射待确认（P2 时定）
- 样本名 → HCC/Cirrhosis 的映射规则（来自 MaxQuant 样本名前缀 / experimentalDesign / 原文样本表）
- 两队列基因名取交集还是并集（建议交集，保证可比）
