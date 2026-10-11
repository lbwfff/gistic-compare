# 设计文档：绝对拷贝数扩展工作流（模块 ABS + 全套统计）

> 状态：**Phase 1 + 2 已实现（2026-10）**——`pairwise_absolute.py`
> （比较统计 + table mode）与 `R/absolute_cn.R`（ASCAT 管线模式，自动
> 降级）均已落地，合成数据端到端测试通过；Phase 3（真实数据）待数据
> 到位后启动（§8 检查清单、§9 开放问题）。
>
> 决策记录：工具选 **FACETS**（否决 HAPSEG+ABSOLUTE 遗留链与 ASCAT，理由 §3）；
> 范围选**全套统计**（用户确认，§5）；语言集成选 **subprocess + Rscript**
> （否决 rpy2，理由 §4.1）。

---

## 1. 目标

在现有"GISTIC 输出比较"工作流之外，增加一条可选的**绝对拷贝数扩展**：
原始 SNP6 数据（CEL → 基因型 + BAF）→ FACETS（纯度/倍性/等位基因特异性
CN/WGD）→ 队列间比较统计。回答的问题：

1. 两队列的 **WGD 频率**是否不同（手册 §4.4：有序结局 0/1/2+）？
2. **纯度、绝对倍性**是否不同（这同时是所有 CN 比较的混杂校正量）？
3. **LOH 占比、拷贝中性 LOH 负担**是否不同？
4. 绝对 CN 空间下的**三级负荷**（P9 的 CN 版）是否不同？
5. **纯度校正**后的逐峰模型（P3 扩展）——纯度衰减 CN 信号，只有绝对 CN
   阶段才能提供这个校正量。

与现有工作流的关系：GISTIC 输出比较模块（P1–P10/PQ/PD）完全不动；ABS 是
独立的可选扩展，明确标注需要原始数据。

## 2. 数据前置条件（实现前必须到位）

| # | 文件 | 说明 |
|---|---|---|
| 1 | 每队列的**原始 CEL** 或已处理的基因型调用 | Birdseed `.calls`/`.confidences`（含 BAF），或 apt-probeset-genotype 输出 |
| 2 | 每探针 **logR + BAF** 表（肿瘤） | CEL → apt/crlmm 生成；TCGA 式 "dat" 文件格式即可 |
| 3 | **匹配 normal**（CEL 或基因型） | FACETS 首选 tumor-normal；无匹配 normal 时 tumor-only 模式可靠性下降（开放问题 §9.1） |
| 4 | SNP6 注释 | SMAP/marker 位置表（项目已有 `genome.info.6.0_hg19` 同源表，可复用） |
| 5 | 样本清单 | sample ↔ CEL/BAF 文件映射 + 队列归属 + 组织来源（blood > solid normal > tumor，手册 §1.1 样本类型优先级） |

**当前缺口**：GSEdata 目录只有 GISTIC 输出（分段 logR），无 BAF/基因型。
用户确认原始 CEL 链路存在，待下载（§8 检查清单第 1 项）。

## 3. 工具决策记录

> **2026-10 修订**：原定 FACETS 为主动工。Phase 2 落地前核对源码发现
> **FACETS 只有测序 count 输入**（`readSnpMatrix` 读 snp-pileup 的
> tumor/normal 读数，0.5.14 与 0.6.2 均确认），并无下表所称的
> "原生支持 SNP6"——该前提有误，予以修正。CUE/CHB 是 SNP6 阵列 +
> Birdseed 数据，主工具切换为 **ASCAT**（阵列出身、直接消费
> logR/BAF/基因型调用）。FACETS 保留为将来拿到测序数据时的备选
> （`tool: facets` 走 count 路径），不作主链路。

| 候选 | 决定 | 理由 |
|---|---|---|
| HAPSEG + ABSOLUTE（原文） | ❌ | 2011–2012 遗留链：BEAGLE（Java）+ 按祖源选单倍型面板 + 老版 R 依赖，安装门槛过高；手册 §4.6 自己警告不同工具 WGD 判定可能不同 |
| **ASCAT**（Van Loo et al. 2010；VanLoo-lab/ascat 3.x） | ✅ **主** | SNP array 出身（Affy SNP6 是其原始平台之一）；直接输入每探针 tumor logR/BAF + germline 基因型（Birdseed .calls 直接给 het SNP）；无匹配 normal 时可走 pooled normal / tumor-only 变体；一次给出 purity（aberrant cell fraction）/ploidy/等位基因特异性分段 |
| FACETS（mskcc，0.6.2） | 备选 | 仅测序 count 输入（readSnpMatrix）；用于将来 WGS/WES/panel 数据的同启发式复算或工具间对照（手册 §4.6 #3 建议） |

## 4. 管线与模块设计

### 4.1 语言集成模式

Python 包装 + `subprocess` 调 `Rscript R/absolute_cn.R`。否决 rpy2（依赖
脆弱）。R 环境/包缺失时模块 **skip with reason**——与 P8 缺 sklearn、P5 缺
D 矩阵的处理一致。R 依赖：`ASCAT`（GitHub VanLoo-lab/ascat，≥3.x）。

### 4.2 CLI 与配置

```
gistic-compare --config cfg.yaml --modules ABS
```

```yaml
absolute_cn:
  tool: facets                  # 预留 ascat
  raw_dir: raw_snp6/            # CEL/基因型+BAF 文件目录（相对 data_root）
  manifest: absolute_manifest.tsv   # sample, cohort, cel/文件路径, 组织来源
  cn_thr_wgd: 3.0               # WGD 启发式 A 的总拷贝阈值
  loh_frac_wgd: 0.30            # WGD 启发式的基因组占比阈值（敏感性分析 0.2-0.5）
  min_purity: 0.30              # 手册 §4.5：纯度过低不可靠
  seed: 42
```

### 4.3 R 脚本契约（`R/absolute_cn.R`）

输入：manifest + 数据目录（CLI 参数）。输出：每样本一行 TSV
`absolute_cn_<cohort>.tsv`：

```
sample, purity, ploidy_abs, wgd_status(0/1/2+), wgd_frac,
loh_frac, cnloh_frac, n_segments_allelic, flag_low_purity
```

内部流程：读基因型+BAF/logR → FACETS `procSample` → purity/ploidy/
等位基因特异性分段 → WGD 启发式（§5）→ 写表。逐样本失败允许（写 flag
而非中断）；日志记录每样本的 FACETS `diplogR`/模型选择。

### 4.4 Python 比较统计（模块 ABS 的后半段）

读入两个队列的绝对 CN 表后：

| 统计 | 方法 | 输出表 |
|---|---|---|
| WGD 频率差 | 2×3 Fisher（0/1/2+）+ 有序 logistic（手册 §4.4） | `ABS_wgd.tsv` |
| purity / ploidy_abs | MW + rank-biserial + 分布图 | `ABS_purity_ploidy.tsv` |
| LOH / CN-LOH 负担 | MW + 分布图 | 同上 |
| 绝对 CN 三级负荷 | P9 机制在 CN 空间重算（事件阈值换算为拷贝数） | `ABS_burden.tsv` |
| 纯度校正逐峰 | P3 扩展：`altered ~ cohort + purity` | `ABS_purity_adjusted_peaks.tsv` |
| summary | 全部要点一行/比较 | `ABS_summary.tsv` |

### 4.5 图

`ABS_wgd_freq`（WGD 频率条形+星号）、`ABS_purity_ploidy`（分布对比）、
`ABS_burden_dists`（绝对 CN 负荷分布）。

## 5. WGD 判定启发式（实现时必须声明 + 敏感性分析）

FACETS 不直接输出 WGD 标签，按等位基因定义实现（手册 §4.6 建议交叉验证）：

- **启发式 A（主，等位基因定义）**：常染色体上"minor CN ≥ 1 且 total CN ≥
  3"（两条亲本拷贝都被复制）的长度占比 ≥ `loh_frac_wgd`（默认 0.30，
  敏感性 0.2–0.5）→ WGD+；占比、重复次数（0/1/2+）由绝对 CN 分段积分。
- **启发式 B（参照）**：ploidy_abs ≥ 3.7 的经典倍性规则，作敏感性对照。
- 两启发式不一致的样本在表中单列 `wgd_disagreement` 列。

## 6. 测试策略

- Python 侧比较统计：合成绝对 CN 表（两组已知 purity/ploidy/WGD 分布）
  纯 Python 测试，**CI 不需要 R**。
- R 脚本（ASCAT 管线）：`pytest.mark.skipif(Rscript 或 ASCAT 不可用)`；
  合成 marker 级 logR/BAF 输入（已知 WGD/纯度/倍性真相，logR 按
  ASCAT 的 gamma 约定压缩）验证 Python subprocess → R → harmonized 表
  → 比较统计的端到端链路。
- 降级路径：R 缺失 → `ABS_skipped.txt` + manifest skip 记录。

## 7. 风险与既有警告（全部来自 R05 §4.5/§4.6）

1. 低纯度（<0.3）样本的 WGD 判定不可靠 → `flag_low_purity` + 可选剔除。
2. 倍性/纯度可辨识性问题（高倍性+低纯度 vs 低倍性+高纯度）。
3. 不同工具 WGD 判定可能不同 → 启发式 A/B 对照 + （可选）ASCAT 交叉验证。
4. WGD 频率 sanity check：应落在 10%–60%（全 0 或全 100% = 算法问题）。
5. tumor-only 模式的可靠性折损（§9.1）。

## 8. 实施检查清单

- [ ] **数据**：下载/整理两队列的基因型+BAF（或 CEL），建 manifest
- [ ] **数据**：确认是否有匹配 normal（决定 tumor-normal vs pooled-normal / tumor-only）
- [x] **环境**：R ≥ 4.x + ASCAT 安装验证（VanLoo-lab/ascat 3.2.0；其
      未使用的 `copynumber` Import 已移除后安装——该包只出现在注释中）
- [x] 实现 `R/absolute_cn.R`（契约见 §4.3；ASCAT 3.2 API：
      `ascat.loadData` → `ascat.aspcf(penalty, seed)` → `ascat.runAscat(gamma)`）
- [x] 实现 `pairwise_absolute.py`（ABS 模块：R 包装 + 比较统计）
- [x] 绘图 `plotting/pairwise_abs.py`
- [x] CLI 注册 + `_skip_reasons` 降级 + config schema
- [x] 测试：合成比较统计（无 R）+ R 依赖端到端（skipif）
- [x] 文档：README 双语、CHANGELOG、stats.md 新增 ABS 章节（含 WGD 启发式声明）
- [ ] 敏感性：`loh_frac_wgd` 0.2/0.3/0.5；启发式 A vs B 对照表（机制已实现
      `wgd_flag_b`/`wgd_disagreement` 列，阈值敏感性待真实数据）

## 9. 开放问题

1. **匹配 normal 是否存在**？CUE/CHB 的 SNP6 运行是否包含 normal 样本？
   （决定 matched vs pooled-normal vs tumor-only）
   **2026-10 部分解决**：tumor-only 路径已实现且验证——het SNP 用样本
   自己的 Birdseed 基因型调用（`het_calls` manifest 列，经 ascat.gg 传
   ASCAT；SNP6 阵列测的就是种系基因型），或回退 ASCAT 内置推断；
   `Germline_LogR` 在 ASCAT 拟合阶段本来就不使用（源码确认），logR 参照
   由 HAPMAP 标准化承担（TCGA 标准做法）。代价：体细胞 LOH 位点在
   tumor 调用里呈纯合 → 从 het 集合消失 → `loh_frac`/`cnloh_frac` 为
   **下限估计**（`germline_mode`/`het_source` 列记录）；logR/WGD 影响
   小。若队列含 normal 样本，可抽做 pooled-normal 对照估偏差。
2. GEO 系列公开的是 Birdseed 调用（含 BAF）还是仅 CEL？决定第 2 步用
   现成调用还是自行 genotype。
3. 纯度校正逐峰（§4.4 最后一行）放 ABS 还是扩展 P3（两模块都要用纯度时
   建议读同一张 `absolute_cn_*.tsv`）。
4. WGD 启发式阈值默认值待 CUE/CHB 实跑后标定（先 0.30，敏感性定稿）。
