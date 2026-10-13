# ABS 真实数据管线状态（Phase 3 留档）

> 状态：**管线机制已通；4 样本 pilot 中 2 个高畸变样本经手动回退拟合
> 成功（EA043/EA098: purity≈1.0、ploidy≈1.5、LOH 95-99%——形态待领域
> 审核）；2 个样本仍失败；全队列未运行**。本文档记录 2026-10 Phase 3
> 真实数据排查的全部问题、实测数据与下一步想法，供后续接手。机制/合成
> 数据部分见 `plan_absolute_cn.md` 与 CHANGELOG。

---

## 0. 实测耗时（修正此前所有估计；样本质量差异导致波动很大）

| 任务 | 实测 | 备注 |
|---|---|---|
| 数据切分（226 样本，7.6GB 流式） | ~45 分钟 | 一次成形 |
| **单样本完整管线**（loadData+GC+aspcf+runAscat） | **~3.5-7 分钟** | CN 清理后；污染时代曾 ~29 分钟 |
| 2 样本 chunk 并行（1 worker/样本） | ~7 分钟墙钟 | 含回退第二次 pass |
| 手动回退（每失败样本） | +1-2 分钟 | 网格扫描秒级 + 单样本手动拟合 |
| 每次 pilot 迭代（4 样本，2 队列） | 7-15 分钟墙钟 | +异常/休眠中断则更久 |
| 全队列预估（226 样本，jobs=4） | **~8-15 小时**（非此前估的 14-28h，也非 4h） | 约 1/4 样本会被预筛/早退 |

教训：此前"每轮 30-60 分钟"是污染时代（het 集 4 倍大）的数字；CN 清理后
单样本已降一个量级。但系统性**低估**仍反复出现——根因：样本质量差异
（噪声样本触发 segmentlength 3 次重试、回退第二次全量 pass）、诊断脚本
自身出错导致返工、以及一次机器休眠中断。

## 0.1 本次会话遇到的问题清单（按时间）

1. **CN 探针占位 BAF 污染 het 集**（birdseed specialSNPs 含 CN 探针，
   BAF=2.0 被当杂合）→ 81% 分段值变负 → GoF 永不达标。已修。
2. **缺 SNP6 GC 校正** → 已修（Zenodo 官方文件）。
3. **倍性网格 [1.5, 5.5] 拒收高 WGD 队列** → max_ploidy=8 已暴露。
4. **并行会话的 P11 半成品把 CLI 改坏**（cli.py 引用未完成的
   figures_p11）→ 绕过 CLI 直接调模块；该会话稍后自行提交完成
   （d292760）。教训：并行会话共享仓库时，每次运行前先 `git status`。
5. **机器休眠杀后台任务**（CHB 完成后 CUE 未跑，输出文件缺失）。
   教训：长任务前台跑，或 `caffeinate -dimsu` 包裹。
6. **回退链按名索引未命名 list**（ascat.aspcf_parallel 返回的
   Tumor_BAF_segmented 无名 → `[[s_name]]` 得 NULL → 手动拟合被静默
   跳过）→ 改位置索引 `[[match(s_name, samples)]]`。
7. **我自己的就地补丁写错列名**（新建全 NA 的 call 列）→ 一度误判
   数据损坏；实际数据未损坏。
8. **诊断脚本自身多次出错**（矩阵 names() vs rownames()、R 跨行字符串
   拼接缺逗号、awk 列统计把字段当样本）→ 多次产生误导信号。
9. **工具输出一度不可靠**（读文件返回陈旧/混乱内容）→ 反复用 git 与
   独立命令核实。

## 0.2 核心未解问题（L4 演进：角点收敛）

手动回退**机制上成功**（EA043/EA098 均产出拟合），但两个 CHB 高 FGA
样本的拟合**全部收敛到网格角点** (rho=1.0, psi≈1.45)：

| 样本 | purity | ploidy | LOH | wgd_frac |
|---|---|---|---|---|
| EA043 (FGA .52) | 1.00 | 1.487 | 99.4% | 0 |
| EA098 (FGA .51) | 1.000 | 1.521 | 95.0% | 0.0003 |

形态：purity 1.00 + 近单倍体 + 全基因组 LOH——即模型把样本解释为
"100% 纯度的近单倍体全 LOH 基因组"。对食管肿瘤不典型（文献近单倍体
克隆罕见），且 purity=1.00 在实体瘤中不合理。

**已排除**：数据问题（logR 信号与 GISTIC 一致：EA043 chr12 段实测
+0.666 vs GISTIC +0.69；BAF 直方图正常）、het 集对齐（het 调用位点
BAF 峰 0.507）、宽表内容（与手动成功运行 diff 一致）、max_ploidy 传参
（回显确认）、探针名 rownames。

**已验证的距离面事实**：EA043 清洗后（216K het 探针）40 段、99.6% 中性
（r≈0, b=0.5）；距离矩阵全局最小在 rho=1.0/psi≈1.5 角点——即距离面
**单调滑向角点**，不存在更优的内部极小。这不是搜索问题，是模型-数据
匹配问题：以 b=0.5 为主、logR 近 0 的分段结构在低 psi + 高 rho 下距离
趋零（模型可把中性段解释为任意纯度下的平衡拷贝）。

**候选解释与想法（未验证）**：
1. **中性段权重 0.05 的效应**：距离公式对 b=0.5 的段降权 20 倍—— Genome
   95% 中性时，全距离由少数畸变段主导；若畸变段在角点的距离恰不比内部
   高多少，角点即胜出。可检查：角点 vs 内部候选的逐段距离分解。
2. **logR 中心化偏移**：HAPMAP 标准化的 logR 中位 ≈ -0.05~-0.08（非 0），
   在 rho=1、psi=1.5 的模型下预期 logR = log2(psi/2) ≈ -0.58 ——
   观测 0 与两者都不完全一致；logR 的基线偏移可能系统性偏向低 psi 解。
   尝试：logR 减去全基因组中位（中心化）后再入 ASCAT。
3. **BAF 的镜像/相位问题**：fastAspcf 输出的分段 BAF 是 per-probe-side
   （0..1），make_segments 直接使用；若应在段级取 min(b, 1-b)（次等位
   基因侧）而未取，nA/nB 的指认会系统性错误。需对照 runASCAT 源码的
   下游用法逐段核实。
4. **EA003/EA043 的样本类型**：GEO 系列可能混有血样/正常——若 pilot 的
   高 FGA 样本实为正常对照，低纯度拟合失败是正确行为；但 EA043 的
   GISTIC 剖析（26 段 |mean|>0.1、含 +0.69 扩增）表明它确实是肿瘤。

**裁决手段（一次 ~10 分钟实验）**：对 EA043 出 sunrise 距离热图
（`ascat.plotSunrise(plot_d, psi, rho)`），直观看距离面在角点与内部
之间是否有可分辨的极小；若内部极小存在但被 percentzero 过滤，可考虑
放宽该单项（min_ploidy 语义）或直接采用手动模式 + GoF 阈值过滤。

## 1. 状态总结

| 环节 | 状态 |
|---|---|
| 数据切分（`scripts/prepare_abs_inputs.py`） | ✅ 完成：CUE 136 / CHB 90 肿瘤样本，与 GISTIC 队列规模精确一致；探针网格对齐；样本名与 GISTIC 峰表一致（纯度校正逐峰模型可直接按名 join） |
| 数据质量 | ✅ 正常：het 调用位点 BAF 峰 0.507、hom 双峰 0/1；LRR 极值形态为未压缩真 log2 |
| GC 校正（ASCAT 3 标准 SNP6 步骤） | ✅ 已接入（Zenodo 14008443 官方 SNP6 参考文件） |
| 自动拟合（runAscat 默认过滤器） | ❌ pilot 全部样本被拒（含 FGA 0.71 的高畸变样本） |
| 手动回退（网格最优 rho/psi 重拟合） | ✅ 机制已通：EA043/EA098 拟合成功（note=manual-fit）；**拟合形态（角点收敛）待领域审核** |
| 队列可拟合比例 | 待定：FGA(\|mean\|>0.15) 分布 p25=0.002 / med=0.077 / p75=0.257——约四分之一样本几乎零畸变（ASCAT 无法拟合是正确行为） |
| ASCAT 优化栈（R/ascat_opt/，13× 方案） | 已 vendor（ascat_fast.R/aspcf_cpp.cpp/ascat_io_fast.R/ascat_parallel.R）；wrapper 已接入 patch_ascat + aspcf_parallel + runAscat_parallel(plot=FALSE)；来源 `GSEdata/ASCAT/ASCAT_SNP6_优化报告.md` |

## 2. 已修复的真 bug（均已实现，附验证方式）

1. **CN 探针占位 BAF 污染 het 集**：birdseed 的 specialSNPs 含 CN 探针，
   其 BAF 列为占位值 2.0（无等位基因）；被标为 het 进入 BAF 分段后，
   段均值被拉向 2.0，`Tumor_BAF_segmented = 1 - bafPCFed` 产生 81% 负值。
   修复：prepare 脚本仅当 marker 的**所有样本 BAF ≤ 1** 才保留其调用；
   已生成数据就地修复（226 文件）。
   验证：清洗后 het 探针数从 ~89 万降到 ~21.6 万（真实 SNP het），分段
   值域回到 [0.235, 0.5]（镜像次等位基因，无负值）。
2. **缺 GC 校正**：ASCAT 3 对 SNP6 的标准步骤（`ascat.correctLogR`）。
   验证：校正后 GC 相关性自动选窗正常（25bp/5kb 等）。
3. **默认倍性网格上限 5.5 拒收高 WGD 队列**：`runAscat(max_ploidy=8)`
   手动验证 EA043 被接受（GoF 97.9）。修复：`--max-ploidy` 暴露为配置
   （`max_ploidy`，本队列用 8）。

## 3. 旧 L4 记录（已被 §0.2 的角点收敛分析取代与细化；过程留档）

**现象**：同一份宽表、同一组参数（gamma=0.55, penalty=70, seed=42,
max_ploidy=8），我手动 R 代码段跑 `ascat.runAscat` **接受**了 EA043
（purity 0.44 / ploidy 3.31 / GoF 97.9），而 wrapper 流水线跑
`ascat.runAscat` **拒绝**（"could not find an optimal ploidy and purity
value"）。

**已排除**：宽表内容（diff 逐字节一致）、`--max-ploidy` 传参（运行时
回显确认 max_ploidy=8）、探针名 rownames、样本列错配（原始文件与切分
产物逐值对照一致）、BAF 标度（直方图正常）。

**关键线索（未消化完）**：
- 手动成功用的 aspcf 对象（RDS）是**污染版 het 集**（CN 探针调用仍在，
  81% 分段值为负）产出的；干净版 aspcf（EA003: 215,911 het 探针，分段
  值域 [0.235, 0.5]）的距离面 EA043 上最优 GoF=97.9 @ (psi=5, rho=0.79)
  ——**高于 80 阈值但仍被拒**。
- 该最优格点处 percentzero=0.000（pass 1 的过滤条件之一）；perczeroAbb
  /percOddEven 也可能未达；且距离面上**严格局部极小**（7×7 邻域全更高）
  可能不存在（噪声面）——即 runASCAT 的接受逻辑组合拒绝了一个 GoF 97.9
  的解。
- 疑点：污染版为何反而"成功"——负值分段（b<0）进入距离矩阵后 GoF 仍
  97.9，说明负值段的贡献方式未被理解；其拟合结果 (0.44, 3.31) 因此
  **不可信**，不能作为管线正确的证据。

**下一步（未执行）**：对照污染版/干净版两套 aspcf 的 runASCAT 完整接受
路径（localmin/optima 逐格检查），锁定拒绝的确切条件。一次运行 ~30 分钟。

## 4. 已实现的回退（端到端已验证：EA043/EA098 fitted）

自动拟合失败的样本：逐失败样本用其自身分段算距离矩阵（`make_segments`
+ `ASCAT:::create_distance_matrix`），在**合理性窗口内**（重算倍性 ∈
[1.5, max_ploidy]，rho ∈ [0.2, 1]——排除退化角点）取最优格点 (psi*,
rho*)，以 `rho_manual/psi_manual` 逐样本重拟合（每失败样本一次，mclapply
并行；跳过画图）。note 列记 "manual-fit"。
依据：ASCAT 官方出口（文档明示 rho_manual/psi_manual 为专家参数）；
manual 分支无条件接受给定 (rho, psi)，绕过 percentzero/局部极小过滤。

实测（pilot 4 样本）：EA043/EA098 均经此路径拟合成功。风险与缓解：绕过
过滤器 → 每个手动拟合的 GoF 需记录，下游 ABS 统计（WGD/purity/ploidy
比较）只应使用 GoF > 80 的拟合；且 §0.2 的角点收敛问题意味着手动拟合
的合理性仍需领域审核（purity=1.0 的解在实体瘤中不合理）。

## 5. 运行成本（2026-10-13 实测修订；覆盖早前估计）

- 单样本完整管线（loadData+GC+aspcf+runAscat+回退）：**~3.5-7 分钟**
  （CN 清理后 het 集 21.6 万；污染时代曾 ~29 分钟）
- aspcf 内部：fastAspcf 随 het 探针数线性缩放；getMad/madWins 随总
  marker 数线性；exactPcf 随纯合 stretch 长度**平方**（长 hom 区样本的
  未知变量）
- 2 样本 chunk 并行（1 worker/样本）：~7 分钟墙钟；pilot 迭代 7-15 分钟
- 全队列（226 样本，jobs=4）：**~8-15 小时**；接 ASCAT 优化栈
  （`GSEdata/ASCAT/ASCAT_SNP6_优化报告.md`，13×，逐位等价已验证）后
  预计 ~1-2 小时——**未接入**，接入步骤见该报告 §4 方式 A/C
- 排查经验：RDS 缓存（`--cache-rds`）使 runAscat 参数实验降到 ~1-2
  分钟；对照实验务必一次设计

## 6. 待决定

1. ~~手动回退路径的端到端验证~~ ✅ 已完成（EA043/EA098 fitted）
2. **角点收敛的裁决**：EA043 sunrise 距离热图（一次 ~10 分钟实验）——
   判定 (rho=1.0, psi≈1.45) 是退化解还是数据真实形态
3. 是否接受"自动 + 手动回退 + GoF 阈值过滤"作为正式流程（写入 stats.md）
4. 更根本的：FGA 中位 0.077 的队列，ABS 的 WGD 频率比较只建立在拟合
   子集上——科学上是否可接受，需领域判断
5. 备选工具链（若 ASCAT 路线放弃）：HAPSEG+ABSOLUTE（原论文链，工具与
   数据类型原生匹配，安装门槛高）；FACETS 不支持阵列数据（源码已确认）
