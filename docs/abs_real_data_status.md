# ABS 真实数据管线状态（Phase 3 留档）

> 状态：**✅ 已解决并端到端验证（2026-10-13）**。根因 = gg het 掩码按
> 清单顺序位置构建，而 ascat.loadData 会按基因组位置重排探针——掩码
> 错位使 BAF 分段输入变成 0/1 垃圾。修复后 4 样本 pilot **4/4 自动拟合
> 成功**（purity 0.44-0.69 / ploidy 1.65-4.57，无一需要手动回退），
> 单样本全管线 **72 秒**。此前所有 purity=1.0 的拟合作废。官方 SNP6
> 规程的两处遗漏（CN/SNP 探针分类中心化、BAF=2→NA）亦已采纳。全文记录
> 排查过程供后续接手；机制部分见 `plan_absolute_cn.md` 与 CHANGELOG。

---

## 0. 实测耗时（修正此前所有估计；样本质量差异导致波动很大）

| 任务 | 实测 | 备注 |
|---|---|---|
| 数据切分（226 样本，7.6GB 流式） | ~45 分钟 | 一次成形 |
| **单样本完整管线**（loadData+GC+aspcf+runAscat+输出） | **~72 秒**（gg 修复后；分阶段：loadData 7 s / GC 9 s / aspcf 25 s / runAscat 2 s）| 掩码错位时代曾 28-40 分钟 |
| 2 样本 chunk 并行（1 worker/样本） | ~7 分钟墙钟 | 含回退第二次 pass |
| 手动回退（每失败样本） | +1-2 分钟 | 网格扫描秒级 + 单样本手动拟合 |
| 每次 pilot 迭代（4 样本，2 队列） | **~4 分钟**墙钟 | 修复后实测；此前 7-40 分钟 |
| 全队列实测（226 样本） | **1 h 31 m**（20 批 × 12 样本、4 并发、ncores=1；内存约束下的安全配置——4 worker×57 样本的粗分块会触发内存压缩风暴，CPU 掉到 20-30%） | 拟合率 226/226 = 100%，零 failedarray；55 个 non-aberrant（CHB 23 / CUE 33，多为 purity=1.0 的平坦基因组——官方 "Likely normal" 形态）|

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
10. **Bug #2（独立复核报告发现，2026-10-13）：`ascat.aspcf_parallel`
   不按样本切分 ascat.gg**——gg 经 `...` 原样透传给单样本子对象，内部
   `gg[, sample]` 恒取第 1 列 → 多样本 worker 进程里除首样本外全部用错
   germline 掩码。生产路径 round-robin 分块（jobs=4 → 4×~57 样本）会
   大面积触发；pilot 因每 worker 恰好 1 样本而侥幸未触发。已修（显式
   ascat.gg 参数 + 按样本切片 + `.realign_gg` 名字键守卫），同进程双样本
   验证 EA098 恢复 0.48/4.54。
11. **`resolve()` 未向量化**——`--cache-rds` 对 ≥2 样本 manifest 直接
   报错（R≥4.3 将 length>1 的 is.na 升级为 error）；pilot 侥幸未开缓存。
   已修（vapply 向量化）。
12. **`ascat.*_sequential` 回退函数被引用但从未定义**——快栈加载失败时
   fallback 会 "could not find function"。已补定义。

> 独立复核报告（GSEdata/ASCAT/ABS_真实数据诊断.md）的核对结论：除
> "修复已提交到 bench/" 一句位置有误（实际在 ~/project/ASCAT/diag/，
> bench 未动）外，其余结论全部核实成立，含对本留档 §3 误判的纠正
> （0.44/3.31 是正确答案，被污染的是角点解）。

## 0.2 角点收敛问题——已解决（根因：gg 掩码错位）

"100% 纯度近单倍体"拟合的根因不是 ASCAT，也不是数据，而是**我们的
het 掩码与 ASCAT 内部探针顺序错位**：

- `ascat.loadData`（stock 与 -fast 版均）将全部矩阵按基因组位置重排
  （官方行为；correctLogR 之后还会按 GC 参考文件取交集再排一次）；
- wrapper 曾把 gg 掩码按**清单/网格顺序**按位置传入 aspcf——网格并非
  严格位置排序（SNP/CN 探针交错，CN_473963@chr1:61735 在网格第 87 万行
  而非第 1 行）——het 标志随机落到同源/未调用探针上；
- BAF PCF 的输入因此变成 0/1 为主的噪声（真实 het 输入中位 0.49），
  分段输出 b≈0.105 处处 → 距离面退化 → 自动拟合全败 → 手动回退在垃圾
  上拟合出 purity=1.0。

**ASCAT 的拒绝自始至终是正确的。** 修复：掩码按键（chr+pos → 探针名）
映射到加载后的对象行序；并顺带修复"缓存命中时 gg 块被跳过、静默切换
predictGermlineGenotypes"的连带 bug。修复后（同一宽表、同一参数）：

- 自动搜索直接通过（runAscat 1.8 s，无警告）；
- aspcf 从 1406 s 降到 **24.6 s**（损坏输入曾使 PCF 触发段长重试与
  二次方行为——"神秘 38 分钟"同源）；
- EA043 = purity 0.44 / ploidy 3.31 / WGD+ 0.57，chr12 大段扩增恢复为
  CN 11(9+2)/CN 10(8+2)，与原始 logR(+0.67)/BAF(0.281) 完全吻合。

## 0.3 全队列结果（2026-10-13，CUE 136 vs CHB 90）

| 指标 | CHB | CUE | 组间检验 | 结论 |
|---|---|---|---|---|
| WGD+ 频率 | 40% (36/90) | 46% (62/136) | Fisher p=0.41, OR=0.80 [0.46-1.37] | **无显著差异** |
| WGD+（仅 purity≥0.30 敏感性） | 32% (23/71) | 37% (37/101) | p=0.63 | 仍无差异 |
| purity 中位 | 0.445 | 0.420 | MW p=0.59 | 无差异 |
| 绝对倍性中位 | 2.04 | 2.21 | MW p=0.58 | 无差异 |
| LOH 占比 | 0.148 | 0.132 | p=0.63 | 无差异 |
| CN-LOH 占比 | 0.013 | 0.028 | p=0.20 | 无差异（CUE 偏高） |
| CN 三级负荷 focal/arm/chrom | — | — | p=0.31 / 0.61 / 0.089 | 无显著（chrom 级 CHB 偏高，r=0.13，不显著） |

两队列在绝对拷贝数层面**全面无显著差异**。55 个 non-aberrant（多为
purity=1.0 的平坦基因组，官方 "Likely normal" 形态）保留在主分析中
（WGD− 处理），并另做了 purity≥0.30 的敏感性分析——结论一致。

## 0.4 HRD 疤痕评分（scarHRD 参考实现，2026-10-13）

vendored scarHRD（MIT，Sztupinszki 2018）对全部 226 样本计算
HRD-LOH / TAI / LST / HRD-sum（7 秒，零错误；表
`ABS_hrd_scores.tsv`）。**结论：组间全阴性**（四列 MW p=0.43-0.91；
HRD+率(≥42) 21% vs 20%，Fisher p=1.0）。**强混杂警告**：HRD-sum 与
分割段数 rho=0.90、与纯度 rho=-0.60——低纯度样本的过分割噪声直接灌分
（HRD_adj 的 ploidy 校正修不了这个维度）；35 个 HRD+(≥42) 样本中 24
个 purity<0.40。解读必须限定 informative + purity≥0.30 层。WGD 阶梯与
HRD_sum 无关（rho=0.07）——混杂源是分割质量而非倍增。

## 1. 状态总结

| 环节 | 状态 |
|---|---|
| 数据切分（`scripts/prepare_abs_inputs.py`） | ✅ 完成：CUE 136 / CHB 90 肿瘤样本，与 GISTIC 队列规模精确一致；探针网格对齐；样本名与 GISTIC 峰表一致（纯度校正逐峰模型可直接按名 join） |
| 数据质量 | ✅ 正常：het 调用位点 BAF 峰 0.507、hom 双峰 0/1；LRR 极值形态为未压缩真 log2 |
| GC 校正（ASCAT 3 标准 SNP6 步骤） | ✅ 已接入（Zenodo 14008443 官方 SNP6 参考文件） |
| 自动拟合（runAscat 默认过滤器） | ❌ pilot 全部样本被拒（含 FGA 0.71 的高畸变样本） |
| 自动拟合 | ✅ **4/4 pilot 样本全部自动通过**（purity 0.44-0.69 / ploidy 1.65-4.57）；修复前 0/4 |
| 手动回退 | 保留为安全网（小样本队列可能仍有个别样本触发），当前 pilot 未触发 |
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
