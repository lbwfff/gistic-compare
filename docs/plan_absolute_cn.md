# 设计文档：绝对拷贝数扩展工作流（模块 ABS + 全套统计）

> 状态：**Phase 1 已实现（table mode，2026-10）**——模块
> `pairwise_absolute.py` 读取两份 harmonized 表并完成全套比较统计；
> Phase 2（`R/absolute_cn.R` 的 FACETS pipeline 模式）与 Phase 3（真实
> 数据）待数据/环境到位后启动（§8 检查清单）。
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

| 候选 | 决定 | 理由 |
|---|---|---|
| HAPSEG + ABSOLUTE（原文） | ❌ | 2011–2012 遗留链：BEAGLE（Java）+ 按祖源选单倍型面板 + 老版 R 依赖，安装门槛过高；手册 §4.6 自己警告不同工具 WGD 判定可能不同 |
| **FACETS**（mskcc，Bioconductor 活跃维护） | ✅ | 原生支持 SNP6（为 TCGA 设计）；一次给出纯度/倍性/等位基因特异性分段/cellular prevalence；活跃维护 |
| ASCAT | 备选交叉验证 | SNP array 出身；如需双工具 WGD 交叉验证（手册 §4.6 #3 建议）可作第二意见 |

## 4. 管线与模块设计

### 4.1 语言集成模式

Python 包装 + `subprocess` 调 `Rscript R/absolute_cn.R`。否决 rpy2（依赖
脆弱）。R 环境/包缺失时模块 **skip with reason**——与 P8 缺 sklearn、P5 缺
D 矩阵的处理一致。R 依赖：`facets`（Bioconductor/GitHub mskcc）。

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
- R 脚本：`pytest.mark.skipif(Rscript 或 facets 不可用)`；合成一个极小的
  等位基因特异性输入（几十个 SNP，已知答案）验证端到端。
- 降级路径：R 缺失 → `ABS_skipped.txt` + manifest skip 记录。

## 7. 风险与既有警告（全部来自 R05 §4.5/§4.6）

1. 低纯度（<0.3）样本的 WGD 判定不可靠 → `flag_low_purity` + 可选剔除。
2. 倍性/纯度可辨识性问题（高倍性+低纯度 vs 低倍性+高纯度）。
3. 不同工具 WGD 判定可能不同 → 启发式 A/B 对照 + （可选）ASCAT 交叉验证。
4. WGD 频率 sanity check：应落在 10%–60%（全 0 或全 100% = 算法问题）。
5. tumor-only 模式的可靠性折损（§9.1）。

## 8. 实施检查清单

- [ ] **数据**：下载/整理两队列的基因型+BAF（或 CEL），建 manifest
- [ ] **数据**：确认是否有匹配 normal（决定 tumor-normal vs tumor-only）
- [ ] **环境**：R ≥ 4.x + `facets` 安装验证（`Rscript -e 'library(facets)'`）
- [ ] 实现 `R/absolute_cn.R`（契约见 §4.3）
- [ ] 实现 `pairwise_absolute.py`（ABS 模块：R 包装 + 比较统计）
- [ ] 绘图 `plotting/pairwise_abs.py`
- [ ] CLI 注册 + `_skip_reasons` 降级 + config schema
- [ ] 测试：合成比较统计（无 R）+ R 依赖端到端（skipif）
- [ ] 文档：README 双语、CHANGELOG、stats.md 新增 ABS 章节（含 WGD 启发式声明）
- [ ] 敏感性：`loh_frac_wgd` 0.2/0.3/0.5；启发式 A vs B 对照表

## 9. 开放问题

1. **匹配 normal 是否存在**？CUE/CHB 的 SNP6 运行是否包含 normal 样本？
   （决定 FACETS 跑法与可靠性声明）
2. GEO 系列公开的是 Birdseed 调用（含 BAF）还是仅 CEL？决定第 2 步用
   现成调用还是自行 genotype。
3. 纯度校正逐峰（§4.4 最后一行）放 ABS 还是扩展 P3（两模块都要用纯度时
   建议读同一张 `absolute_cn_*.tsv`）。
4. WGD 启发式阈值默认值待 CUE/CHB 实跑后标定（先 0.30，敏感性定稿）。
