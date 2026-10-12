# ABS 真实数据管线状态（Phase 3 留档）

> 状态：**管线机制已通，自动拟合被 ASCAT 过滤器拒绝，手动回退已实现但
> 未端到端验证**。本文档记录 2026-10 Phase 3 真实数据排查的全部结论，
> 供后续接手。本文只述真实数据（CUE/CHB SNP6）；机制/合成数据部分见
> `plan_absolute_cn.md` 与 CHANGELOG。

---

## 1. 状态总结

| 环节 | 状态 |
|---|---|
| 数据切分（`scripts/prepare_abs_inputs.py`） | ✅ 完成：CUE 136 / CHB 90 肿瘤样本，与 GISTIC 队列规模精确一致；探针网格对齐；样本名与 GISTIC 峰表一致（纯度校正逐峰模型可直接按名 join） |
| 数据质量 | ✅ 正常：het 调用位点 BAF 峰 0.507、hom 双峰 0/1；LRR 极值形态为未压缩真 log2 |
| GC 校正（ASCAT 3 标准 SNP6 步骤） | ✅ 已接入（Zenodo 14008443 官方 SNP6 参考文件） |
| 自动拟合（runAscat 默认过滤器） | ❌ pilot 全部样本被拒（含 FGA 0.71 的高畸变样本） |
| 手动回退（网格最优 rho/psi 重拟合） | 已实现；EA043 手动验证可接受（purity 0.44 / ploidy 3.31 / GoF 97.9），**端到端未验证** |
| 队列可拟合比例 | 待定：FGA(\|mean\|>0.15) 分布 p25=0.002 / med=0.077 / p75=0.257——约四分之一样本几乎零畸变（ASCAT 无法拟合是正确行为） |

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

## 3. L4 未解之谜（下一步排查的起点）

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

## 4. 已实现的回退（wrapper 内，待端到端验证）

自动拟合失败的样本：逐失败样本用其自身分段算距离矩阵（`make_segments`
+ `ASCAT:::create_distance_matrix`），取全局最优格点 (psi*, rho*)，以
`rho_manual/psi_manual`（逐样本向量，混合模式：自动成功的样本保持自动
结果）重跑 `ascat.runAscat`。note 列记 "manual-fit"。
依据：ASCAT 官方出口（文档明示 rho_manual/psi_manual 为专家参数）；
manual 分支无条件接受给定 (rho, psi)，绕过 percentzero/局部极小过滤。

风险与缓解：绕过过滤器 → 每个手动拟合的 GoF 需记录，下游 ABS 统计
（WGD/purity/ploidy 比较）只应使用 GoF > 80 的拟合。

## 5. 运行成本（规划全队列时用）

- aspcf 对 176 万 marker：~20 分钟/样本（瓶颈）；runAscat 网格：~3 分钟
- 226 样本、`jobs: 4` 并行 ≈ **14-28 小时**（约四分之一样本 non-aberrant
  会快速退出，实际更短）
- 每次 pilot 迭代 30-60 分钟——**排查时务必一次设计对照实验**，勿逐假设
  推进

## 6. 待决定

1. 手动回退路径的端到端验证（一次 ~10 分钟实验：4 个高 FGA 样本）
2. 是否接受"自动 + 手动回退 + GoF 阈值过滤"作为正式流程（写入 stats.md）
3. 更根本的：FGA 中位 0.077 的队列，ABS 的 WGD 频率比较只建立在拟合
   子集上——科学上是否可接受，需领域判断
4. 备选工具链（若 ASCAT 路线放弃）：HAPSEG+ABSOLUTE（原论文链，工具与
   数据类型原生匹配，安装门槛高）；FACETS 不支持阵列数据（源码已确认）
