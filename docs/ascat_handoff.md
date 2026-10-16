# ASCAT 交接说明（服务器执行 ↔ 本地下游分析）

> 状态：**ASCAT 执行在服务器，下游分析在本地/服务器均可**（2026-10 起）。
> 本文是两台机器之间的接口契约：服务器产出什么、放在哪里、列怎么定义，
> 以及本地包对表的校验方式。设计与实现背景见 `plan_absolute_cn.md`，
> 真实数据的排查记录见 `abs_real_data_status.md`。

## 0. 为什么要分两台机器

ASCAT 在 PC 上勉强：226 样本实测 **1h31m**（20 批 × 12 样本、4 并发、
`ncores=1` 的内存约束配置），粗分块会触发内存压缩风暴（见
`abs_real_data_status.md` §0）。因此 **ASCAT 移到服务器**，本包只做下游。

分工的原则：**跨机器边界只传两张表**。ASCAT 的算法、版本、阈值都留在
服务器侧；本包不关心它们怎么跑出来，只关心表是不是符合契约。

## 1. 分工一览

| 机器 | 做什么 | 代码 |
|---|---|---|
| **服务器** | CEL → Birdseed → 每样本 logR/BAF/het → ASCAT → 两张表 | `scripts/prepare_abs_inputs.py`、`R/absolute_cn.R`（可选用加速栈 `R/ascat_opt/`） |
| **本地** | 读两张表做队列间统计与画图（模块 ABS）+ HRD 分数 | `pairwise_absolute.py`（table 模式）、`plotting/pairwise_abs.py`、`scripts/run_hrd_scores.py` + `R/scarHRD/` |

本地默认**不需要 R + ASCAT**（`analysis.absolute_cn.auto_run: false`）：
表缺失时 ABS 写 `ABS_skipped.txt` 并给出指引，而不是尝试 shell 出去跑 R。
`scarHRD` 那一步仍需要本地有 R（只读表，不跑 ASCAT）。

## 2. 服务器运行步骤

1. **切分输入**（cohort 级宽表 → 每样本三件套 + manifest）：

   ```bash
   python scripts/prepare_abs_inputs.py --data-root GSEdata --out-dir GSEdata/raw_snp6 \
     --cohort cue:CUE/cue.lrr_baf.txt:CUE/birdseed.calls.txt \
     --cohort chb:CHB/CHB.lrr_baf.txt:CHB/birdseed.calls.txt
   ```

   产出 `<cohort>/<sample>__logr.tsv.gz`、`__baf.tsv.gz`、`__het.tsv.gz`
   与 `manifest.tsv`（列：`sample, cohort, logr, baf [, germline_logr,
   germline_baf | het_calls]`）。只保留常染色体 1–22。

2. **跑 ASCAT**（每队列一次；`--ascat-opt-dir` 默认 = 脚本旁的
   `R/ascat_opt/`，存在且 Rcpp 可用时自动加载 bit-identical 加速栈，
   否则回退 stock ASCAT）：

   ```bash
   Rscript R/absolute_cn.R --manifest GSEdata/raw_snp6/manifest.tsv \
     --cohort cue --out-dir GSEdata/CUE \
     --gc-file GSEdata/ref/GC_AffySNP6_102015.txt --ncores 8 \
     --cn-thr-wgd 3.0 --loh-frac-wgd 0.30 --min-purity 0.30 \
     --ploidy-wgd 3.7 --penalty 70 --gamma 0.55 --seed 42
   ```

   - **SNP6 必须给 `--gc-file`**：HAPMAP 归一化的 logR 不做 SNP6 GC 校正时
     网格搜索达不到 GoF 阈值（CUE/CHB pilot 实测）。
   - 逐样本失败不中断：失败样本写成 `fit_ok=0` + `note`。
   - 服务器上分批并行：按样本切 manifest 成多个块，多个进程各跑一块
     （Python 侧 `auto_run` 的 `jobs` 分块逻辑做的是同一件事；服务器上
     直接 shell 分批更直观）。

3. **回放表**：把 `absolute_cn_<cid>.tsv` + `absolute_seg_<cid>.tsv`
   放到 `<data_root>/<cohort>/`（即该队列 GISTIC 输出目录旁）；或集中放
   一处，用 `analysis.absolute_cn.tables_dir` 指过去。

4. （推荐）随表放一份 `absolute_meta.json`（见 §5）记录版本与阈值。

5. **本地跑下游**：

   ```bash
   gistic-compare --config config.local.yaml --modules ABS
   python scripts/run_hrd_scores.py --config config.local.yaml   # 可选 HRD
   ```

## 3. 表契约（冻结）

表头必须精确匹配（大小写、无多余空格），制表符分隔，`sample` 与
GISTIC 侧的样本名一致。

### 3.1 `absolute_cn_<cid>.tsv`（每样本一行）

**必需**（下游统计真正读取的列，缺一列即 skip 并指出缺哪列）：

| 列 | 含义 |
|---|---|
| `sample` | 样本 ID（与 GISTIC 侧一致） |
| `purity` | 肿瘤纯度（aberrant cell fraction），0–1 |
| `ploidy_abs` | 绝对倍性（平均总拷贝数） |
| `wgd_status` | 0/1，WGD 启发式 A 的判定 |
| `wgd_frac` | 常染色体"minor CN ≥ 1 且 total CN ≥ `cn_thr_wgd`"的长度占比 |
| `loh_frac` | 常染色体 LOH（`cn_minor == 0`）长度占比 |
| `cnloh_frac` | 拷贝中性 LOH（`cn_total == 2 且 cn_minor == 0`）长度占比 |
| `flag_low_purity` | purity < `min_purity`（默认 0.30） |

**可选**（有就记录进 `ABS_summary.tsv` / fit 统计，缺了只打印一行 note）：
`fit_ok`、`note`、`germline_mode`、`het_source`、`n_segments_allelic`、
`wgd_flag_b`、`wgd_disagreement`、`non_aberrant`。

### 3.2 `absolute_seg_<cid>.tsv`（等位基因分段）

必需：`sample, chrom, start, end, cn_total, cn_major, cn_minor`。

- `chrom`：`chr1`…`chr22`（**也接受 `1` / `1.0`**，读入时归一到 `chrN`）；
  必须落在配置的基因组内（hg19 `chr1`–`chr22`），否则 skip。
- 单位：`cn_*` 是**绝对拷贝数**（不是 log2）；`cn_minor ≤ cn_major`，
  `cn_total = cn_major + cn_minor`；`start/end` 是 hg19 bp。
- 只保留常染色体（`prepare_abs_inputs.py` 已经保证）。

### 3.3 本地校验

`pairwise_absolute._read_abs_table` 在读表时校验：必需列齐全、`chrom`
归一后合法、数值列可解析、表非空。**不通过 → `ABS_skipped.txt` 写明是哪个
表缺哪一列并指向本文档**，不会带病计算。校验常量在
`pairwise_absolute.py` 顶部（`_CN_REQUIRED` / `_SEG_REQUIRED`）。

## 4. 口径与解释（写图注/论文时照抄）

- `loh_frac` / `cnloh_frac` 是**长度占比**，不是样本比例。
- **tumour-only 时 `loh_frac`/`cnloh_frac` 是下界**：体细胞 LOH 位点在肿瘤
  自身基因型里表现为纯合，会掉出 het 集。`germline_mode` / `het_source`
  记录实际模式：`matched`（测得 germline BAF）/ `tumor_only` +
  `birdseed_calls`（肿瘤自身基因型代理）/ `tumor_only` + `inferred`。
- **WGD**：启发式 A（等位基因定义）为主，B（`ploidy_abs ≥ ploidy_wgd`）为
  参照，不一致记 `wgd_disagreement=1`。单快照无法可靠计数多轮 WGD，
  `wgd_status` 只作二元。
- `ploidy_abs` 同时是 CN 比较的混杂校正量。`purity ≥ 0.99` 的平坦基因组
  样本（`qc_tier=non-informative`）与 `fit_ok=0`/手动回退样本要在敏感性
  分析里单独看（CUE/CHB 实测 55 例 non-aberrant）。
- 低纯度样本（`flag_low_purity`）会被排除出 WGD 主检验的敏感性重跑，
  差异记在 `<comparison> [high-purity only]` 行。
- 队列间 CN/LOH 的比较检验清单见 `pairwise_absolute.py` 的模块 docstring
  与 `stats.md`。

## 5. 变更控制（服务器换 ASCAT/R 版本时）

- 表 schema 是**冻结**的。若必须改：先改 `pairwise_absolute.py` 顶部的
  必需列常量与本文档，再改服务器脚本，顺序不要反。
- 推荐随表放 `absolute_meta.json`（本地不解析，仅审计用，缺失不影响运行）：

  ```json
  {"ascat_version": "3.2.0", "r_version": "4.4.1", "genome": "hg19",
   "penalty": 70, "gamma": 0.55, "cn_thr_wgd": 3.0, "loh_frac_wgd": 0.30,
   "min_purity": 0.30, "ploidy_wgd": 3.7,
   "gc_file": "GSEdata/ref/GC_AffySNP6_102015.txt",
   "samples": 136, "fit_failed": 0, "generated": "2026-10-16"}
  ```

- 换版本后建议先跑一次本地 ABS，确认没有 skip（`ABS_skipped.txt` 不存在）
  再引用结果。

## 6. ASCAT 结果独有的分析（在包内做，不需要本地 ASCAT）

以下两条已在模块 ABS 内实现（表 + 图），不需要服务器上跑任何东西：

- **LOH × GISTIC-Del 交叉**（`ABS_loh_del_cross.tsv`，图 `ABS_loh_del_cross`）：
  每个 (样本 × Del 峰, q≤0.25) 分成 del+LOH / del-only / LOH-only /
  neither，且 LOH-only 再拆成 **拷贝中性**（`cn_total==2 & cn_minor==0`，
  无剂量变化 —— GISTIC 的 Del track 结构上看不见）与 CN-loss（有剂量丢失
  但低于 GISTIC 的幅度调用）两类。口径：`loh_cross_core`（默认 narrow
  Peak Limits）+ `loh_cross_cover`（默认 0.5 覆盖阈值）。
- **WGD 分层的 LOH 比较**（`ABS_purity_ploidy.tsv` 的 `<comp> [WGD+ only]` /
  `[WGD- only]` 行，图 `ABS_loh_wgd_strata`）：WGD 后 minor CN 普遍 ≥ 1，
  `loh_frac` 系统性下降，合并检验的方向由两队列的 WGD 构成决定（CUE/CHB
  实测反向）；引用 LOH 结论时读分层行，不读合并行。
- **HRD**：`scripts/run_hrd_scores.py`（vendored `R/scarHRD/`，MIT）读同一
  批表算 HRD-LOH / TAI / LST；同样受 purity/ploidy 影响，引用时给分层。

## 7. 相关文件

| 文件 | 角色 |
|---|---|
| `R/absolute_cn.R` | 服务器侧 ASCAT runner（写出两张表） |
| `R/ascat_opt/*` | runner 的 bit-identical 加速栈（可选，自动加载） |
| `scripts/prepare_abs_inputs.py` | 服务器侧输入切分 |
| `pairwise_absolute.py` | 下游统计（table 模式 + 可选 `auto_run` pipeline 模式） |
| `plotting/pairwise_abs.py` | ABS 图 |
| `scripts/run_hrd_scores.py` + `R/scarHRD/` | HRD 分数 |
| `tests/synth_absolute.py` | 合成表（不依赖 R，schema 回归用） |
| `plan_absolute_cn.md` / `abs_real_data_status.md` | 设计与真实数据排查记录 |
