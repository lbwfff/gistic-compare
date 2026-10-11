# GISTIC_compare（中文文档）

比较不同 GISTIC 2.0 结果之间的差异与相似性。**单一 CLI、两个 workflow**（配置字段切换），配置驱动、结果目录隔离、图表出版级。英文说明见 [README.md](README.md)。

| workflow | 模块 | 拓扑 | 用途 |
|---|---|---|---|
| **benchmark**（默认） | A–G | 星型：各实验 vs 对照 | 评估实验表现、相对对照的特有信号 |
| **pairwise** | P1–P10, PQ, PD, CF, GS, ABS | 两两实验比较（无对照） | P1–P3 检出/一致性（基于 GISTIC 峰/调用），P4/P5 区域级组间差异（P4 用 CBS 分段、P5 用 GISTIC marker 矩阵），P6 臂级 broad 景观，P7 altered 内部幅度，P8 多变量整体可分性，P9 三级负荷+wGII，P10 负担受控逐峰置换，PQ p 值质检，PD 批次/基线诊断，CF 混杂敏感性，GS 振幅差分，ABS 绝对拷贝数扩展（ASCAT，需原始数据，缺失自动跳过） |

## 安装

```bash
pip install gistic-compare       # PyPI（发布后）
# 或从源码：
pip install .
# 开发模式：
pip install -e ".[dev]"
```

Python ≥ 3.10；每个队列需要 hg19 的 GISTIC 2.0 输出。

## 快速上手

```bash
gistic-compare --list-modules
gistic-compare --config my_config.yaml --all
# 指定输出目录、只跑部分模块：
gistic-compare --config my_config.yaml --modules P1,P4 --outdir out/
```

- 结果写到 config 同级的 `results/<run_id>/{tables,figures}`（或 `--outdir`）；
- **每次运行自动写 `run_manifest.json`**（config 快照、包/依赖版本、git commit、各模块状态与耗时）；
- `--modules` 子集校验所属 workflow（benchmark 给 P1 会报错并提示）；
- 某模块失败不中断其余模块（除非 `--fail-fast`）；退出码 0 = 全部成功。

### 队列目录要求

每个队列目录下放 GISTIC 输出，命名为 `<prefix>.<后缀>`：

| 文件 | 用途 |
|---|---|
| `<prefix>.all_lesions.conf_90.txt` | A, C, P1–P3, P7/P10（peak 层 CN 值） |
| `<prefix>.amp_genes.conf_90.txt` / `.del_genes.conf_90.txt` | D, G, P3 |
| `<prefix>.scores.gistic` | B |
| `<prefix>.tumors.seg` | E, P4, P6, PD（臂级调用/基线）、P7–P9（回退） |
| `<prefix>.sample_seg_counts.txt` | E, P10（样本 QC） |
| `<prefix>.regions_track.conf_90.bed` | A（校验） |
| `<prefix>.D.cap1.5.mat` | P5, GS, P7–P9（可选） |
| `<prefix>.broad_values_by_arm.txt` | P6（可选；缺失时由 seg 推导） |
| `<prefix>.amp_qv.conf_90.txt` / `.del_qv.conf_90.txt` | P6（可选，仅展示） |

> **解析注意**：GISTIC `all_lesions.conf_90.txt` 中每个峰有两行（摘要行 + `- CN values` 行），解析器只取摘要行（0/1/2 调用）。

## 配置

### benchmark（`examples/config.benchmark.yaml`）

```yaml
workflow: benchmark
data_root: /path/to/data        # 队列目录的上级
run_id: exp_vs_luad
comparison_mode: control_vs_all # control_vs_all | all_pairs | explicit
control: luad
cohorts:
  cue: {label: CUE, path: CUE, prefix: cue}
  luad: {label: TCGA_LUAD, path: TCGA_LUAD, prefix: LUAD}
```

| mode | 生成的比较 | kind |
|---|---|---|
| **control_vs_all**（星型） | `{control}_vs_{exp}` | 全部 `unpaired` |
| `all_pairs` | 所有无序对 | `paired_with` 或样本交集 |
| `explicit` | 仅 `comparisons` 列表 | 按列表 |

### pairwise（`examples/config.pairwise.yaml`）

```yaml
workflow: pairwise
run_id: my_run
comparison_mode: experiment_pairs   # 两两自动展开；无 control
cohorts:
  expA: {label: ..., path: ..., prefix: ...}
  expB: {label: ..., path: ..., prefix: ...}
analysis:
  paired_frac: 0.8                  # 样本交集/min(n) ≥ 0.8 → paired
metadata:                           # 可选（P3 用）；缺省 → P3 跳过
  path: meta/sample_meta.tsv
  sample_id: sample
  covariates: [sex, age_group]
```

- **kind 自动判定**：加载后按双方 `calls` 样本 ID 交集，`交集 / min(n_a,n_b) ≥ paired_frac` → `paired`；
- pairwise **禁止**设置 `control`；benchmark **必须**设置 `control`。

## 论文/研究专属参数 → preset 文件

所有与具体研究绑定的数值都外置为 YAML preset（**不在代码里**）：

```yaml
analysis:
  regional:
    loci_set: presets/loci_sets/we_ea.yaml     # 命名位点 + band 标记
paper_benchmark:
  enabled: true
  panel: presets/reference_panels/we_lung.yaml  # 模块 G 参考值
```

复制 `presets/` 下的示例、改成自己研究的数值、在 config 里指向它即可。**不配置则完全通用**：无位点标注、无论文图、无硬编码假设。

## 参考表（cytoband / refgene）

hg19 的 cytoband 与 refgene 精简表随包分发（`src/gistic_compare/data/refs/`），config 未提供 `analysis.regional.cytoband_file` / `refgene_file` 时自动使用。从 GISTIC 的 `hg19.mat` 重新生成：

```bash
python scripts/build_refs.py --mat /path/to/hg19.mat
```

config 里也可指向 GISTIC `.mat` 或同列名的 TSV。P5/GS/P7–P9 用到各队列的 `<prefix>.D.cap1.5.mat`（GISTIC 自己产生）；缺失时 P5/GS 自动跳过并提示，P7–P9 回退到 seg 推导的窗口值。

## 模块

### benchmark（A–G）

| ID | 内容 |
|----|------|
| A | 峰区间 overlap、Jaccard、q/freq 相关 + bootstrap CI、置换 null、共有峰成员表 |
| B | scores.gistic → 1 Mb 网格重采样、轨迹相关热图、top 差异 bin |
| C | 配对样本 0/1/2 调用一致率 / Kappa / McNemar（仅 paired）；频率相关 |
| D | 基因 / cytoband 集合 Jaccard、driver presence（Amp≤30Mb / Del≤10Mb） |
| E | 每样本 FGA、segment QC、组间检验 |
| F | 实验特有信号占比（vs 对照）：count / freq / −log10(q) / G-score 四种权重 |
| G | 对照参考面板打分（**数值来自 preset**；未配置则跳过） |

### pairwise（P1–P10, PQ, PD, CF, GS, ABS）

> **定位区分**：P1–P3 回答**检出/一致性**——「两队列各自显著（GISTIC 峰）之后，峰清单/调用是否重合」，不是组间差异检验；低幅度/宽/低频的组间差异天生成不了峰。**P4–P7** 回答**组间差异**——「两队列在同一区域的 CN 分布是否不同」，不依赖 GISTIC 峰（P7 的 peak 层除外）。P8 是全局多变量摘要，**P9 是所有逐位点结果的前置背景层**，**P10 是负担受控的峰检验**，**PQ/PD 是质量闸门**（膨胀与批次/基线诊断）。

| ID | 定位 | 内容 |
|----|------|------|
| **P1** | 检出/一致性 | 匹配对 + Δq/Δfreq、单侧特有峰、三层 tier Jaccard、UpSet、逐样本 **OncoPrint**(union 峰 0/1/2 调用) |
| **P2** | 检出/一致性 | 样本级调用一致（仅 paired）：exact agreement、weighted κ、McNemar+FDR；unpaired 自动跳过 |
| **P3** | 建模 | 逐峰 `altered ~ cohort [+ 单协变量]` logistic、BH-FDR；无 meta 跳过 |
| **P4** | **组间差异（CBS 分段）** | `tumors.seg` 切固定窗（默认 1 Mb）→ 二值化 → Fisher + MW + FDR、区域合并、null 校准、**最小可检出差异（MDE）表** |
| **P5** | **组间差异（marker 矩阵）** | 与 P4 同构，输入换成 `D.cap1.5.mat` 重采样到同一网格 |
| **P6** | **臂级 broad 景观** | 每臂 Amp/Del 频率差 + Fisher/FDR、每样本非整倍体评分（broad 事件数/臂占比）+ MW；输入优先 `broad_values_by_arm.txt`，缺失时按 seg 推导（来源记录在表内） |
| **P7** | **altered 内部幅度** | 频率（Fisher）与**altered 样本内部的幅度**（MW + rank-biserial）双轴：高幅度比例（≥0.9）、freq×amp 四象限分类（freq_driven / amp_driven / both / neither）；窗口/band（D 矩阵）+ peak（all_lesions `- CN values` 行）三层 |
| **P8** | **多变量整体可分性** | PERMANOVA（伪 F + 置换 p，配对设计层内置换）、PCA、交叉验证 AUC（L2/弹性网 + 置换零分布）、bootstrap 稳定性选择；特征默认 1 Mb 窗（与峰检出无关），union peaks 模式有循环偏倚仅作描述 |
| **P9** | **三级负荷评分 + wGII** | 每样本 focal/arm/chromosome 三级事件计数（幅度加权 \|C\|，合并样本内秩归一化）+ 不依赖 GISTIC 的 wGII（倍性自适应基线）；回答「哪个队列整体更乱、差异来自 focal 还是 broad」——是解读所有逐位点结果的前置背景层 |
| **P10** | **负担受控逐峰置换** | union 峰 × 样本二值矩阵上**同时固定行和（每样本峰数）与列和（每峰频率）**的置换检验；对照 naive p 直接量化负担假阳性（burden_flag）；含着丝粒/端粒峰与高 segment 样本 QC |
| **PQ** | **p 值质检** | 自动读取本次运行的 P3/P4/P5/GS/P7/P10 全部 p 值集合：Q-Q 图 + 基因组膨胀因子 λ=median(z²)/0.4549（λ>1.1 = 系统性偏倚，先解决再引用结果） |
| **PD** | **λ 膨胀的批次/基线诊断** | 每样本基线（中性段 seg_mean 中位数）检查 + Δmedian 空间结构轨迹 + **样本内中心化敏感性**：P4 式检验在原始/中心化两套数据各跑一遍，λ 对比 + 每窗口可信度分层（robust / baseline_only / centered_only / ns） |
| **ABS** | **绝对拷贝数扩展（可选）** | ASCAT（SNP6 logR/BAF，R 包装 `R/absolute_cn.R`）：WGD 频率（Fisher+OR）、纯度/绝对倍性/LOH/CN-LOH（MW+rank-biserial）、CN 空间三级负荷（P9 机制移植）、纯度校正逐峰模型；无原始数据自动跳过 |
| **CF** | 混杂敏感性 | 只读 P4/P5 表 + 每队列边际分布：E-value、最小解释效应 δ*、偏差上界/去偏下界（**定界**而非校正） |
| **GS** | 振幅差分 | 逐窗/逐 band G-score 组间差 + 置换检验（normal 或 MC） |

**统计引用层级（P4/P5）**：逐窗 FDR 仅探索（1 Mb 窗高度相关）；**cytoband 级**为预设单元、较严谨；**`*_paper_loci`（preset 位点）**为先验假设检验、可确证性引用；`*_differential_regions` 有选择偏差，仅候选。

**CF 假设与简化**：只有每队列边际 → 只能**定界**；二值混杂偏差上界 `|Δp|·δ`、`δ* = D/Δp`，超过 `plausible_max_binary`（默认 0.20）判稳健；age 仅中位数近似；多混杂可加保守近似；EGFR/KRAS/TP53 可能是中介/碰撞变量，过校正风险。

## 绘图规范

由 `theme.py` 驱动（绘图实现按模块拆在 `src/gistic_compare/plotting/`）：

- 双栏宽 7.0 in 上限；Arial（回退 Helvetica/DejaVu）；`pdf.fonttype=42`（文字可编辑）；
- 同时输出 PDF + PNG（600 dpi）；Okabe-Ito 色盲友好；未列出队列按「与已用色距离最远」从 tab10 取色（可用 `output.figure.cohort_colors` 覆盖）；
- Δ 方向一律用**队列颜色**（不借用 Amp/Del 色）；相关矩阵检测到负相关自动切 `RdBu_r`；
- 统计标注统一 `r = 0.49 (95% CI 0.10–0.75, n = 13), p = 0.012`；多重比较用 `theme.bh_adjust`（BH）；比例数据用 Wilson 区间（`theme.bar_with_ci`）；小样本叠加原始散点。

## 关键参数

- `q_thresholds: [0.25, 0.05]` — 主 0.25，敏感性 0.05
- 峰匹配三层（同 chr + 同类型 + overlap≥50bp）：**moderate**（主，单侧 rec≥50%）/ **strict**（双侧 rec≥50%）/ **loose**（任意 overlap≥50bp）
- q 可比性：组内 rank-normalized −log10(q)（Blom z）
- 置换：保留长度多重集 + 每染色体峰计数；p=(k+1)/(n+1)
- `genes.max_wide_len_mb_amp: 30` / `max_wide_len_mb_del: 10`
- P4：`bin_size: 1000000`、`amp_thr: 0.1`、`del_thr: -0.1`、`fdr: 0.05`、`null_splits: 20`、`seed: 42`
- P6：`broad.amp_thr: 0.3`、`del_thr: -0.3`、`min_frac: 0.5`（seg 推导时的臂调用阈值）
- P7：`amplitude.amp_thr: 0.1`、`del_thr: -0.1`、`high_amp: 0.9`、`high_del: -0.9`（对齐 GISTIC 自己的振幅阈值 0/0.1/0.9）
- P8：`multivariate.features: grid|cytoband|peaks`、`n_perm: 999`、`n_auc_perm: 100`、`n_boot: 50`、`penalty: l2|elasticnet`
- P9：`burden.arm_thr: 0.25`、`arm_min_frac: 0.70`、`arm_match_diff: 0.10`、`focal_thr: 0.25`、`focal_broad_overlap_max: 0.5`、`sensitivity: true`（81 组合参数网格）
- P10：`permtest.cn_thr: 0.25`、`n_perm: 10000`、`stat: zprop|wald`、`max_segments: 2000`
- PQ：无需配置（自动发现本次运行的 p 值集合）
- PD：`diag.neutral_band: 0.25`、`smooth_windows: 20`
- ABS：`absolute_cn.tool: ascat`、`raw_dir: raw_snp6/`、`manifest: absolute_manifest.tsv`（列：sample, cohort, logr, baf [, germline_logr, germline_baf | het_calls]；路径相对 manifest；het_calls = 样本自己的 Birdseed 基因型 AA/AB/BB，无 matched normal 时的 tumor-only het 代理）、`cn_thr_wgd: 3.0`、`loh_frac_wgd: 0.30`（敏感性 0.2–0.5）、`min_purity: 0.30`、`ploidy_wgd: 3.7`、`penalty: 70`、`gamma: 0.55`；表缺失但 R+ASCAT 可用时自动走管线模式（Rscript R/absolute_cn.R），否则跳过

## 测试

```bash
pytest        # 97 个测试：解析/配置/区间/统计 + 两条 workflow 端到端（合成数据）+ 各模块单测；
              # ABS 的 R/ASCAT 管线测试在无 Rscript/ASCAT 时自动跳过
ruff check src tests
```

测试自带的合成 GISTIC 输出生成器在 `tests/synth.py`（含 MATLAB v7.3 `D.cap1.5.mat`），不依赖任何外部数据。

## 许可

MIT，见 [LICENSE](LICENSE)。学术使用请引用 [CITATION.cff](CITATION.cff)。
