# GISTIC_compare（中文文档）

比较不同 GISTIC 2.0 结果之间的差异与相似性。**单一 CLI、两个 workflow**（配置字段切换），配置驱动、结果目录隔离、图表出版级。英文说明见 [README.md](README.md)。

| workflow | 模块 | 拓扑 | 用途 |
|---|---|---|---|
| **benchmark**（默认） | A–G | 星型：各实验 vs 对照 | 评估实验表现、相对对照的特有信号 |
| **pairwise** | P1–P5, CF, GS | 两两实验比较（无对照） | P1–P3 检出/一致性（基于 GISTIC 峰/调用），P4/P5 区域级组间差异（P4 用 CBS 分段、P5 用 GISTIC marker 矩阵），CF 混杂敏感性，GS 振幅差分 |

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
- `--modules` 子集校验所属 workflow（benchmark 给 P1 会报错并提示）；
- 某模块失败不中断其余模块（除非 `--fail-fast`）；退出码 0 = 全部成功。

### 队列目录要求

每个队列目录下放 GISTIC 输出，命名为 `<prefix>.<后缀>`：

| 文件 | 用途 |
|---|---|
| `<prefix>.all_lesions.conf_90.txt` | A, C, P1–P3 |
| `<prefix>.amp_genes.conf_90.txt` / `.del_genes.conf_90.txt` | D, G, P3 |
| `<prefix>.scores.gistic` | B |
| `<prefix>.tumors.seg` | E, P4 |
| `<prefix>.sample_seg_counts.txt` | E |
| `<prefix>.regions_track.conf_90.bed` | A（校验） |
| `<prefix>.D.cap1.5.mat` | P5, GS（可选） |

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

config 里也可指向 GISTIC `.mat` 或同列名的 TSV。P5/GS 另需各队列的 `<prefix>.D.cap1.5.mat`（GISTIC 自己产生）；缺失时模块自动跳过并提示。

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

### pairwise（P1–P5, CF, GS）

> **定位区分**：P1–P3 回答**检出/一致性**——「两队列各自显著（GISTIC 峰）之后，峰清单/调用是否重合」，不是组间差异检验；低幅度/宽/低频的组间差异天生成不了峰。**P4/P5** 才回答**组间差异**——「两队列在同一区域的 CN 分布是否不同」，不依赖 GISTIC 峰。

| ID | 定位 | 内容 |
|----|------|------|
| **P1** | 检出/一致性 | 匹配对 + Δq/Δfreq、单侧特有峰、三层 tier Jaccard、UpSet、Δfreq Wilcoxon+FDR、1 Mb 轨迹相关 + Δ 轨迹 |
| **P2** | 检出/一致性 | 样本级调用一致（仅 paired）：exact agreement、weighted κ、McNemar+FDR；unpaired 自动跳过 |
| **P3** | 建模 | 逐峰 `altered ~ cohort [+ 单协变量]` logistic、BH-FDR；无 meta 跳过 |
| **P4** | **组间差异（CBS 分段）** | `tumors.seg` 切固定窗（默认 1 Mb）→ 二值化 → Fisher + MW + FDR、区域合并、null 校准 |
| **P5** | **组间差异（marker 矩阵）** | 与 P4 同构，输入换成 `D.cap1.5.mat` 重采样到同一网格 |
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

## 测试

```bash
pytest        # 27 个测试：解析/配置/区间/统计 + 两条 workflow 端到端（合成数据）
ruff check src tests
```

测试自带的合成 GISTIC 输出生成器在 `tests/synth.py`（含 MATLAB v7.3 `D.cap1.5.mat`），不依赖任何外部数据。

## 许可

MIT，见 [LICENSE](LICENSE)。学术使用请引用 [CITATION.cff](CITATION.cff)。
