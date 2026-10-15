# 遗留问题清单（known issues）

来源：2026-10 全面代码审查（算法正确性 / 统计方法 / 数据对齐 / 优化机会）。
本文件记录**尚未修复**的全部问题；已修复项见 `CHANGELOG.md` 的 `[Unreleased] / Fixed`
（含 P10 置换检验的两个有效性缺陷、模块 F/PD/P10 的数据错位、P9 秩归一化、PQ p=0、
P4 空模型等 —— 如曾引用旧版 P10 的 `p_perm`/`q_perm`/`burden_flag`/`lambda_perm`，
请用修复后的版本重跑）。

状态图例：`[ ]` 未开始。每项标注：位置（文件/函数）、问题、影响、建议修法。

优先级排序原则：**会被引用的数字的正确性风险 > 结果解释误导 > 静默污染/崩溃 >
防回归基建 > 性能整洁**。建议执行顺序：B1 → B2 → C1+C2 → B3+C3+C4 → D 档随缘穿插
→ E1 单独立项；若近期出图投稿，把 C2-1 提前。

---

## B 档 — 数字正确性（直接影响“哪些结果进论文”）

### [ ] B1. 显著性口径统一
同一份报告里“显著”的定义目前不一致，直接决定引用哪些峰/带。

1. **FDR 校正范围不一致**：`pairwise_model.run_module_p3`（P3_peak_models 全表一次
   BH）、`pairwise_amplitude._add_fdr`（P7 全表一次 BH）、`gs_score._test_groups`
   （GS 全表一次 BH，且 `gs_score.run_module_gs` 的 `GS_summary.tsv` 把混合 FDR 的
   计数写在只标第一个 comparison 的行里）vs `pairwise_regional._differential`、
   `pairwise_arm.run_module_p6`、`pairwise_permutation.run_module_p10`（按
   comparison 各自 BH）。多比较运行时两种口径并存。
   → 建议：统一为按 comparison 分别 BH（与 P4/P6/P10 对齐），GS summary 改 per
   comparison 行。
2. **`pairwise_regional.analyze_windows` 的 `sig_rows`（P4_significant_bands）**：
   `sig_amp`/`sig_del` 用原始 p<0.05，`sig_mw_fdr` 用 FDR<0.05，且 0.05 硬编码、
   不随 `rcfg["fdr"]`。→ 统一用 FDR + 可配置阈值。
3. **P10 双标尺**：`pairwise_permutation.run_module_p10` 的 `burden_flag` 用原始
   `p_perm` 对比 `pcfg["fdr"]`，而 `n_perm_sig` 用 `q_perm` —— 同一张表两种标尺。
   → 建议 flag 改用 `q_perm`（或两个字段并列明示）。
4. **硬编码 0.05**：`plotting/pairwise_p10.py`（灰/彩分色与标题，不随
   `pcfg["fdr"]`）、`pairwise_model.run_module_p3` 的 `n_crude_fdr05`/`n_adj_fdr05`、
   `confounding._sensitivity_rows` 的 `sig` 判定。→ 全部走配置阈值。

### [ ] B2. 协变量编码统一
1. **连续变量被 one-hot**：`pairwise_burden._rank_and_test` 与
   `pairwise_permutation._covariate_design` 用
   `pd.get_dummies(m[cov_names].astype("string"))` 把**所有**协变量（含年龄）按
   取值 one-hot；`pairwise_model._design_matrix`（P3）则把数值且 nunique>2 的列当
   连续变量。小样本下 P9/P10 参数爆炸且与 P3 不可比。
   → 建议：抽取统一的 `build_covariate_matrix()`（数值→连续、低基数→哑变量），
   三处共用。
2. **静默 NaN**：`pairwise_burden._rank_and_test` 的 `d = model.dropna(subset=[metric])`
   不同步删 `cov` 行，`_ols_cohort` 形状不匹配被 `except` 吞掉返回 NaN —— 协变量
   存在且指标含 NaN 时该 metric 的结果整列静默 NaN。→ dropna 对 y/cov 同步，或
   except 里至少打印一次警告。

### [ ] B3. P9/P6 事件语义与文档不符（影响具体数值）
1. **焦点事件排除口径**：`pairwise_burden._drop_broad_covered` 接收的 broad_spans
   在 `_cohort_scores` 里由 `arm_out`/`chrom_ev` 的**整臂/整染色体区间**构建；
   文档（docstring + docs/stats.md）承诺的是"合格片段（qualifying spans）"——
   `_arm_events` 算出的 `spans` 列从未被使用。实际效果：宽型臂上的焦点事件一律
   被排除。→ 用 spans 列或改文档；两者的 S_focal 数值不同。
2. **P6 归一化口径**：`pairwise_arm._broad_from_seg` 的 frac = 覆盖长度占比
   （`n_amp / n_cov`），文档说"臂长占比"；且 `>= min_frac` 而 P9 `_arm_events` 是
   `frac > min_frac`。→ 统一口径与边界。
3. **0.5/0.5 边缘**：`_broad_from_seg` 中 frac_amp = frac_del = 0.5 时同一臂同时
   记 Amp 和 Del（`broad_frac` 双计）。→ 排他判定。
4. **区间约定**：半开区间下 `start==end` 的单标记片段长度为 0，在
   `pairwise_regional._window_values`、`seg_qc.compute_fga`、
   `score_track.resample_track` 中被**静默丢弃**；`_window_values` 中一个 NaN
   `seg_mean` 片段会污染整个 样本×窗口 值（该样本被静默移出该窗口的检验）。
   → 约定写入文档 + NaN 片段按非 NaN 段加权。

### [ ] B4. 模块 G（paper_reproduce）语义对齐
1. **genome% 定义不一致**：实现是显著峰 wide-peak 并集 bp（amp/del 不互斥），
   `presets/reference_panels/*.yaml` 定义是"exclusively amplified/deleted
   segments 占比"；`paper_reproduce` 的 ±3pp 匹配门限比的是两种量。→ 对齐实现或
   对齐 preset 措辞。
2. top-region 的 Amp/Del 归属按列表成员判定（同一 band+min_freq 出现在两个列表时
   误判为 Amp）；`presets.load` 总是构造 `label == band` 的三元组，region label
   特性是死的。→ 遍历列表时记录 side。
3. `paper_reproduce`：无 experiment 队列时静默回退为 control 与自己比较。→ 报错。

---

## C 档 — 静默污染源与误导性输出

### [ ] C1. IO 与矩阵读取鲁棒性（"不报错但数据脏"一类）
1. **seg sample ID 不 strip**：`io_gistic.load_seg` 不 strip `Sample` 值，而
   `load_sample_seg_counts`（:448）和 all_lesions（:161）都 strip。ID 带空白/CR 时
   `seg_qc` 的 merge 静默产生 NaN、FGA 相关性悄悄少样本。→ 解析时统一 strip。
2. **all_lesions 元数据列精确匹配**：`io_gistic.load_all_lesions`（meta_cols 注释
   称"容忍变体"，实际精确匹配）——表头一旦改名，多出的列会作为"样本"混进
   `- CN values` 表，污染 P7 的 CN 层。→ 按前 9 列位置判定元数据列。
3. **marker_matrix 脆弱点**：`load_marker_matrix` 硬编码 `#refs#/p`（ir/jc）与
   `#refs#/o`（segval），应沿 `D/dat` 的 MATLAB 对象引用解析（换 GISTIC 版本/写入
   顺序可能静默读错对象）；`marker_window_matrix` 的 `np.empty(M)` 在首段不从
   marker 0 开始时注入未初始化内存。已用真实 `D.cap1.5.mat` 验证 ir/jc 为 0 基 ✓，
   但 tests/synth 同样写 0 基（循环验证）→ 建议保留一份真实文件布局断言或 fixture。
4. **加权口径**：`marker_matrix.marker_window_matrix` 的窗口均值是**标记数加权**，
   P4 `_window_values` 是**长度加权**，而 P4 marker 层（`pairwise_regional`
   模块 docstring，原 `pairwise_marker`）称两者 "directly comparable"。
   → 统一或改措辞。
5. **`pairwise_regional._load_refgene` 不做 chr 前缀归一**（`_load_cytobands` 做）：
   GISTIC `.mat` 里 `rg.chr` 为裸数字时 refgene 全部丢失、静默回退到 GISTIC 基因
   表。→ 与 `_load_cytobands` 同样归一。
6. `io_gistic._arm_key` 静默丢弃 Xp/Xq（与 `pairwise_arm._arm_label` 是两处独立
   实现，行为需一致并在文档写明"仅常染色体臂"）。

### [ ] C2. 绘图层修复批（按误导性排序）
1. **paper-loci 面板只画第一个 comparison**：`plotting/pairwise_p4.plot_p4_paper_loci`
   与 `plot_p4_paper_list_heatmap_*` 取 `comparison` 列的 `iloc[0]`，≥2 个非配对
   比较时其余数据被静默丢弃且文件名无比较标识 —— 论文图直接误导。→ 按比较分面
   或文件名带 comparison。
2. **MW 着色跨方向**：`plotting/pairwise_p4._plot_p4_genome`（criterion="mw"）
   Amp/Del 两面板共用非方向特异的 `p_mw_fdr` 着色，注解却声称"某队列更多"的方向性
   结论。→ 分方向表述或改用方向性检验着色。
3. **空面板崩溃**：`plotting/pairwise_p7.plot_p7_volcano`、`plotting/gs.plot_gs_volcano`
   两面板全空时 `info` 未定义 → UnboundLocalError；环脚注画在最后一个循环轴上。
4. **分面共享刻度标签**：`plotting/pairwise_p1.plot_p1_jaccard`、
   `plotting/benchmark_a.plot_jaccard` 的 tick 标签取 `iloc[0]`（moderate 层 Amp
   行），Del 面板被标注 Amp 的峰数。→ 按 (comparison, ptype) 生成标签。
5. **NaN 边缘崩溃**：`plotting/pairwise_p1.plot_p1_delta_track` 的
   `np.nanmax(...) or 1.0`（NaN 为真值，`or` 救不了）→ ylim ValueError；
   `plot_p1_set_summary` 的无-moderate 回退分支在重复索引上 `float()` → TypeError。
6. `plotting/pairwise_p1.plot_p1_delta_freq` 标注原始 p（表里有 BH `fdr` 列，与
   "BH-FDR annotations" 的图形规范不符）。
7. `plotting/benchmark_a`、`plotting/benchmark_c`：`np.polyfit` 的 NaN 守卫只查 x
   不查 y → ValueError。
8. `plotting/benchmark_g.plot_paper_genome_pct`：control-only 运行（过滤后 0 行）
   → IndexError。
9. `plotting/oncoprint`：cytoband 缺失时 `str(np.nan)` 为真值字符串 "nan"，行标签
   渲染成 `nan`（`_core._clean_label` 的回退永不触发）。
10. `plotting/pairwise_p9.plot_p9_dists` 的无-comparison 回退分支 KeyError /
    comp 名与 `A_vs_B` 不匹配导致 p 注解丢失。
11. `plotting/confounding` 的 tornado 图 `max(allx)` 空序列 ValueError。
12. `pairwise_p6`/`pairwise_p9`/`pairwise_pd` 的箱线图 tick 直接用 cohort id，
    其余图用 `theme.cohort_label` —— 统一。

### [ ] C3. P3 union-key 链式合并的重复计数
`pairwise_model.build_union_keys` 的传递链（A~B、B~C ⇒ A~B~C 一个 key）在 ≥3 队列
时可把**同一队列**的两个峰并进一个 key；`run_module_p3` 的 alt 表（key × sample）
不按 (key, sample) 聚合（对比 P10 的 `groupby("sample").mean()`），同一样本在同一
key 下贡献多行，`n_samples = len(y)` 把行数当样本数，独立性假设被破坏。
→ 按 (peak_key, sample) 聚合 + n_samples 用唯一样本数。

### [ ] C4. 配置契约统一
1. **两套契约**：benchmark 模块裸 `ana[...]` 缺键即 KeyError 并记 failed
   （`seg_qc`、`score_track`、`sample_calls`、`novelty`、`peak_overlap`、
   `gene_sets`），pairwise 模块全部有默认值。→ `load_config` 注入默认值或文档
   标明必填键。
2. `config._setup_output`：配置里显式设了 `output.tables/figures` 时 `--outdir`
   被静默忽略（与 CLI 帮助的"输出覆盖"语义相反）。→ 明确优先级并写入帮助文本。
3. `config.build_comparisons` 的 `all_pairs` 只检查 a→b 方向的 `paired_with`，
   `expB: {paired_with: expA}` 静默产生 unpaired。→ 双向检查。
4. ~~`cli` 帮助文本仍写 "pairwise: P1-P5"~~ **已修复**（见 `CHANGELOG.md`
   [Unreleased]）：现为 `benchmark: A-G vs control; pairwise: P1-P4, P6-P11`。
5. `run_manifest.git_commit` 在进程 CWD 而非配置/项目目录执行 git —— 可能记录
   无关仓库的 commit。

### [ ] C5. 其他小项
- `score_track`（模块 B）：`np.nan_to_num(na=0.0)` 把无覆盖 bin 当 0 混入
  top-diff 列表；输出文件名 `B_score_grid_1mb.tsv` 硬编码（`bin_size` 可配置）。
- `gene_sets`：`gene_long`/`set_long` 的基因行重复写出；与 `paper_reproduce` 的
  基因名清洗规则不一致（同一 driver 可被 G 找回而 D 找不回）。
- `scripts/build_refs.py`：`chr` 字段为数值（如 1.0）时生成 `chr1.0`。
- `confounding._as_fraction`：`v > 1` 视为计数 —— 比例恰为 1.0 时有歧义。
- `pairwise_burden.run_sensitivity` 的 verdict：全部 beta=0 时
  `sign_consistent=True`（语义边缘）。
- `io_gistic.load_genes`/`load_arm_qv` 对异常布局的容错路径未测试。

---

## D 档 — 方法学抛光与性能

### [ ] D1. 方法学
1. **P8 CV 泄漏**：`pairwise_multivariate._clean_features` 在全数据上做中位数填补
   + z 标准化后才进 CV（无监督泄漏，高维时轻度抬高 AUC）→ 清理步骤移入 sklearn
   Pipeline。另：docs/stats.md 称置换零模型"uses the same CV folds"，非配对时
   StratifiedKFold 按置换后标签重新分折（结论不变，改措辞）；配对 PERMANOVA 的
   统计量未按 pair 分块（置换是限制性的，p 有效但功效略损）。
2. **贪婪峰匹配**：`intervals.match_peaks` 按 overlap 贪心，计数可低于最大二分图
   匹配，Jaccard 略低估（两队列同向处理，可接受）→ 可作敏感性选项（Hungarian）。
3. P4 空模型计满 grid 窗口、主分析跳过 na/nb=0 窗口 —— 口径有细微差异。
4. P10 的 `cn_thr` 等复现参数无敏感性网格（P9 有 81 组合网格）。

### [ ] D2. 性能（按收益排序）
1. **参考表缓存**：`pairwise_regional._load_band_genes`/`_genes_by_band`（~4 万
   基因的 Python 循环）被 8 个模块各调一次，无缓存 → `lru_cache`（按 path/mtime）。
2. **marker matrix 运行内缓存**：1.6M 标记的 HDF5 在 P4 的 marker 层（P5_*）、
   P7/P8/P9/GS 中被反复解析
   （`_matrix_for_cohort` + `_sample_ids_for` 各一次）→ 按路径缓存，顺带消除双读
   的顺序契约隐患。
3. `novelty._aggregate_gscore`：O(peaks × scores) 全表扫描 → 按 (type, chrom)
   分组切片。
4. `plotting/_core._annotate_labeled_points`：每次调用 `fig.canvas.draw()`（600dpi
   全渲染）、循环内反复 `np.asarray`、O(k²) 碰撞扫描。
5. `plotting/benchmark_d` 排序键内 `mat.loc` 查找（O(G²)）。
6. `pairwise_regional._window_values` / `score_track.resample_track`：每段扫整条
   染色体的 bin 数组 → searchsorted。
7. `marker_matrix.marker_window_matrix`：每样本 Python 循环；`np.where(cnt>0,
   s/cnt)` 先除后选产生除零警告。
8. `pairwise_diag._lambda_rows` 每行 `_lambda` 算两次（PQ 已改一次）；
   `pairwise_diag` raw/centered 两遍 `_window_values` 的几何可复用。

### [ ] D3. 死代码清理
`genomes.get_genome`；`presets.loci_names`；`pairwise_model.Path_exists`、占位
`y`、`names_ext`、`mask`；`pairwise_p3` 的 `med` 与被忽略的 `cov_summary` 参数；
`theme.annotate_proportion`（无调用者）；`pairwise_p1` 的空守卫分支；`pairwise_arm`
的 `_ = (la, lb)`；`pairwise_burden._arm_events` 的 `spans` 列（若 B3-1 不采用则删）。

### [ ] D4. 测试余量
- 不等队列规模的端到端（当前 8×8 唯一布局）；空/退化 cohort 边界（0/1/2 样本、
  阈值下 0 峰）；`test_power_manifest.test_mde_matches_wald_power` 中与实现同源的
  代数往返部分。
- `pytest` 未设 `filterwarnings = error`（35 个运行时警告静默通过）——历史上的
  退化输入路径全靠它们暴露。
- P2 配对路径仍无执行覆盖（并入 E1）。

---

## E 档 — 功能排期（单独立项）

### [ ] E1. 真正的配对检验
现状：P4/P5/P6/P7/P9/P10/PD/GS 对 `kind=paired` 比较直接跳过（消息已明示"仅实现
非配对检验"）；P8 支持 pair 内置换 + GroupKFold（含 `_complete_pairs` 部分配对
处理）；P2 McNemar；P3 忽略样本内相关（已文档化）。

目标：P4/P9 增加 pair 内标签置换/符号秩版本；P10 配对版置换（注意
`_binary_matrix` 目前对共享样本 ID 去重保首列 —— 配对语义需重新设计为同一患者
两列或配对差值）；P6/P7/GS 明确配对或标注保守性（同对内正相关使非配对检验偏保守，
不会虚报显著但功效降低）；补齐配对路径测试。

### [x] E2. 绝对拷贝数扩展（已实现，Phase 1+2）
WGD/purity/ploidy——**ASCAT 3.x，R 后端**（`R/absolute_cn.R`；原定 FACETS
经源码核对不支持 SNP6 阵列输入，已切换，决策记录见
`docs/plan_absolute_cn.md` §3）。比较统计/table 模式 + 管线模式均已落地，
合成数据端到端通过；**待 Phase 3**：真实数据（CEL/Birdseed 下载、manifest
整理、`loh_frac_wgd` 敏感性定标）。
