# Changelog

All notable changes to this project are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Fixed
- **ABS — `ascat.aspcf_parallel` 不按样本切分 ascat.gg（独立复核发现，
  全队列阻塞项）**：gg 经 `...` 原样透传给单样本子对象，内部 `gg[, sample]`
  恒取第 1 列——多样本 worker 进程中除首样本外全部用错 germline 掩码
  （生产路径 jobs=4 → 4×~57 样本会大面积触发；4 样本 pilot 因每 worker
  恰好 1 样本而未暴露）。修复：显式 `ascat.gg` 参数 + 按样本切片 +
  `.realign_gg` 名字键防御守卫（正确输入时 no-op）。同进程双样本验证：
  EA098 恢复 purity 0.48 / ploidy 4.54。独立复核：
  `GSEdata/ASCAT/ABS_真实数据诊断.md`（结论全部核实成立）。
- ABS：`resolve()` 向量化（此前 `--cache-rds` 对 ≥2 样本 manifest 因
  R≥4.3 的 length>1 is.na error 直接崩溃）；补定义
  `ascat.aspcf_sequential` / `ascat.runAscat_sequential`（快栈不可用时
  的 stock 回退此前被引用但未定义，会直接崩）。
- **ABS — 退化拟合根因修复（gg het 掩码与 loadData 排序错位）**：
  `ascat.loadData`（stock 与 fast 版均）把矩阵按基因组位置重排，
  `ascat.correctLogR` 再按 GC 参考文件取交集（本队列丢失 ~4.3k 探针），
  而 wrapper 的 het 掩码按清单/网格顺序按位置传入——het 标志随机错位，
  BAF 分段输入变成 0/1 噪声（分段 b≈0.105 处处，真实中位 0.49），导致
  自动拟合全败、(rho=1, psi=1.45) 角点解与 purity=1.0 的退化拟合。
  掩码现按键（chr+pos → 探针名）映射到加载后对象的行序，并对 GC 交集
  取并；同时修复缓存命中时 gg 块被跳过、静默切换
  predictGermlineGenotypes 的连带 bug。修复后 4 样本高 FGA pilot
  **4/4 自动拟合**（purity 0.44-0.69 / ploidy 1.65-4.57，无需手动回退；
  EA043 chr12 扩增恢复为 CN 11(9+2)/10(8+2)，与原始 logR/BAF 吻合），
  单样本墙钟 28-40 分钟 → **72 秒**（aspcf 1406 s → 24.6 s）。此前所有
  purity=1.0 的 ABS 结果作废重算。详见 `docs/abs_real_data_status.md`。
- ABS：pilot 规模（每队列 <3 拟合样本）时 post-pipeline 汇总对空
  metrics 表崩溃（KeyError: comparison）——加守卫跳过。

### Changed
- **Modules P4/P2 merged with their duplicated counterparts** — 23 module IDs
  → 21 (two fewer to run and read):
  - **P5 is no longer a module**: it became a *layer of P4*
    (`pairwise_regional.run_marker_layer`). Running `--modules P4` runs the
    1 Mb-window pipeline on post-CBS `tumors.seg` (the `P4_*` tables) **and**,
    when every cohort provides `D.cap1.5.mat`, on GISTIC's own marker matrix
    (the `P5_*` tables, same statistics, `analyze_windows` + `write_p4_like`
    reused directly). The layer is switchable with
    `analysis.regional.marker_matrix: false` and reports, but does not fail,
    when the marker files are absent.  `pairwise_marker.py` is gone and the
    P5 figures are drawn by `figures_p4` from `res["marker"]`.
  - **C and P2 are one implementation** (`sample_calls.run_calls`, the old
    `pairwise_calls.py` deleted): two workflow forms of the same call
    concordance, differing only in the comparison set and table prefix.  P2
    no longer skips unpaired pairs outright — it now writes the
    matched-pair **frequency** layer (plus `P2_frequency_scatter`), leaving
    the paired-only agreement/McNemar layers empty; and both forms now share
    the same `summary` columns (`kind`, `n_paired_samples` included).
    Two numeric/robustness side-effects of unifying: `a_freq`/`b_freq` are
    now computed on finite calls only (C silently counted NaN calls as
    0-weighted), and the per-comparison FDR loop is NaN-safe.
- **CLI**: `--workflow` help lists `pairwise: P1-P4, P6-P11` (P5 folded into
  P4); `--skip-reasons` no longer lists P5.
- **ABS — 采纳官方 SNP6 预处理规程**（ASCAT ExampleData/
  ASCAT_fromCELfiles.R）：logR 按 CN 探针/SNP 探针**分类中心化**（本队列
  实测类间偏移 +0.01~+0.11 log2）；**BAF=2 占位符转 NA**（CN 探针无
  基因型，NA 为官方预期状态；NA 卫生检查相应只按 logR 丢行）。
- ABS：新增 `--ncores`/`--max-ploidy`/`--cache-rds`（两级 RDS 缓存：
  GC 后 / ASPCF 后，内容键控）与分阶段耗时日志（`[abs][t]`）。

### Added
- **ABS — HRD genomic-scar scores via vendored scarHRD** (MIT reference
  implementation, `R/scarHRD/`): per-sample HRD-LOH / TAI / LST /
  HRD-sum / ploidy-adjusted score in `ABS_hrd_scores.tsv`
  (`scripts/run_hrd_scores.py`), plus four distribution figures
  (density + ECDF per metric).  Confounding documented: HRD-sum ~
  segmentation count rho = 0.90 and ~ purity rho = -0.60 — interpret on
  the informative tier with a purity floor; between-cohort comparison
  null (MW p = 0.43-0.91).
- **Modules P9/P11 (pairwise) — CNV landscape figures**:
  - **P9 arm / chromosome event-frequency figures** (`P9_arm_events`,
    `P9_chrom_events`): per cohort, the fraction of samples carrying each
    arm / chromosome broad event, Amplifications above and Deletions below
    the baseline, the two cohorts side by side (dark = cohort A, light =
    cohort B), alternating chromosome shading, and a green marker for a
    cohort difference (two-sided Fisher per unit, BH-FDR within each level
    and direction).  New table `P9_event_freq.tsv`.
  - **Module P11 — genome-wide focal q-value landscape**
    (`P11_focal_landscape`): each cohort's `scores.gistic` ``-log10(q)``
    track over the whole genome, Amplifications above / Deletions below a
    `q` threshold line, with a **non-linear** (symlog, linear to q=1e-5)
    axis so a single extreme deletion cannot squash the near-threshold
    peaks, equal-height panels, peak genes annotated (protein-coding
    preferred), and the lighter cohort filled.  Descriptive overlay only —
    q-values are within-cohort, so no between-cohort test.  New table
    `P11_focal_peaks.tsv`.
- **Module ABS (pairwise) — absolute copy-number comparisons (Phase 1,
  table mode)**: reads harmonized FACETS-output tables per cohort
  (`absolute_cn_<cid>.tsv`, `absolute_seg_<cid>.tsv`) and runs the full
  between-cohort suite: WGD frequency (Fisher + odds ratio), purity /
  ploidy / LOH / CN-LOH (MW + rank-biserial), P9's three-level burden
  transplanted into copy-number space, and purity-adjusted per-peak
  logistic models.  Low-purity QC with a high-purity-only sensitivity
  rerun.  The R/FACETS pipeline mode ships in Phase 2 (see
  `docs/plan_absolute_cn.md`); until then the tables must exist per
  cohort (synthetic generator in `tests/synth_absolute.py`).  Outputs:
  `ABS_wgd.tsv`, `ABS_purity_ploidy.tsv`, `ABS_cn_scores.tsv`,
  `ABS_cn_tests.tsv`, `ABS_purity_adjusted_peaks.tsv`, `ABS_summary.tsv`.
- **Module ABS Phase 2 — ASCAT pipeline mode**: if the harmonized tables
  are missing but `analysis.absolute_cn.{raw_dir,manifest}` exists and
  Rscript + the ASCAT package are available, the R wrapper
  `R/absolute_cn.R` runs first (per cohort; per-sample failure
  tolerance; logged to `ABS_r_pipeline.log`) and emits the harmonized
  tables consumed by table mode.  Tool decision corrected: **FACETS
  takes sequencing counts only** (no SNP-array logR/BAF input — the
  original plan §3 premise was wrong), the primary tool is now **ASCAT
  3.x** (native SNP-array ancestry; Birdseed genotypes give the het
  SNPs; tumor-only mode without germline columns).  WGD heuristics A
  (minor >= 1 & total >= 3 fraction, binary status) and B (ploidy
  >= 3.7) with a `wgd_disagreement` column are unchanged and
  tool-independent.  Synthetic ASCAT inputs (known-truth marker-level
  logR/BAF in ASCAT's gamma-convention) drive the end-to-end test;
  R tests skip cleanly when Rscript/ASCAT is unavailable.
- **Module ABS — tumour-only het handling**: ASCAT 3.2 crashes without
  germline input (the het genotype matrix must be provided explicitly);
  the wrapper now supports a `het_calls` manifest column — the samples'
  own Birdseed genotype calls (AA/AB/BB, or 0/1/2) passed to
  `ascat.aspcf(ascat.gg=...)` — with a fallback to ASCAT's built-in
  `predictGermlineGenotypes`.  Mode provenance is recorded in the
  per-sample table (`germline_mode`, `het_source`) and in
  `ABS_summary.tsv` (`germline_mode_a/b`, `het_source_a/b`).  Tumour-only
  makes `loh_frac`/`cnloh_frac` lower-bound estimates (somatic-LOH loci
  drop out of the het set); total-CN/WGD signals are unaffected.  On the
  synthetic truth the tumour-only run recovers purity/ploidy/WGD on par
  with the matched run.  Also fixed: `_cn_scores` crashed on
  zero-segment tables (all-samples-failed cohorts); leftover "FACETS"
  strings replaced.
- **Module PD (pairwise) — baseline/batch diagnostic**: per-sample
  baselines (median of neutral segments), genome-wide shift-structure
  tracks, and the centering-sensitivity experiment — P4-style window
  tests re-run on within-sample centered segments, with lambda compared
  raw-vs-centered and every window assigned a trust layer (robust /
  baseline_only / centered_only / ns).  Answers *why* PQ flags inflation
  and which loci survive baseline removal.  Outputs: `PD_lambda.tsv`,
  `PD_shift_regions.tsv`, `PD_shift_track.tsv`, `PD_baselines.tsv`,
  `PD_summary.tsv`; figures: baseline distributions, shift track, lambda
  bars.
- **P9 parameter sensitivity** (`analysis.burden.sensitivity: true`):
  81-combo grid over `arm_min_frac` x `arm_match_diff` x `focal_thr` x
  `focal_broad_overlap_max` with a per-metric sign-consistency verdict —
  the robustness statement for the burden conclusions.  Outputs:
  `P9_sensitivity.tsv`, `P9_sensitivity_summary.tsv`; figure:
  `P9_sensitivity`.
- **Module P10 (pairwise) — burden-controlled permutation on peaks**
  (adopted from Yuan et al., Cancer Cell 2018): union-peak x sample binary
  matrix (direction-matched `- CN values` at +/-0.25) with a margin-
  preserving permutation (row sums = per-sample peak counts, column sums =
  per-peak frequencies; vectorised checkerboard swaps) so only peak-
  specific preferential hitting is testable. Every peak carries its naive
  p alongside — the naive-vs-permutation gap quantifies burden-induced
  inflation (`burden_flag`). Peak QC: invariable, centromere/telomere and
  high-segment samples/peaks dropped. Outputs: `P10_peak_perm.tsv`,
  `P10_summary.tsv`; figures: volcano, naive-vs-permutation scatter, Z
  histogram.
- **Module PQ (pairwise) — p-value QC**: reads the run's saved P3/P4/P5/
  GS/P7/P10 tables and reports, per p-value set, the genomic inflation
  factor lambda = median(z^2)/0.4549 plus a Q-Q figure; P10's naive and
  permutation p are kept as separate sets (their lambda gap = the removed
  burden inflation). lambda > 1.1 is flagged as systematic bias.
  Outputs: `PQ_inflation.tsv`; figures: `PQ_qq_*`.
- **Module P9 (pairwise) — three-level SCNA burden + wGII** (adopted from
  Yuan et al., Cancer Cell 2018): per-sample copy-number events counted
  separately at focal, arm and chromosome scales (amplitude-weighted |C|,
  focal events from the GISTIC marker matrix with broad-overlap exclusion),
  rank-normalized within the pooled two-cohort set, plus the
  GISTIC-independent wGII (ploidy-adjusted). Cohort effects via OLS beta +
  Mann-Whitney. Outputs: `P9_burden_scores.tsv`, `P9_burden_tests.tsv`,
  `P9_crosscheck.tsv`, `P9_summary.tsv`; figures: burden distributions,
  effect forest, focal-vs-broad decomposition.
- **Module P6 (pairwise) — arm-level broad landscape**: per-arm Amp/Del
  frequency differences (Fisher + BH) and per-sample aneuploidy metrics
  (broad-amp/del counts, fraction of arms altered; Mann-Whitney +
  rank-biserial). Uses GISTIC's `broad_values_by_arm.txt` when every cohort
  provides it, otherwise derives arm calls from `tumors.seg` (thresholds
  `analysis.broad.*`; the source is recorded in every table). Arm q-values
  from `amp_qv`/`del_qv` are attached display-only when present.
  Outputs: `P6_arm_diff.tsv`, `P6_aneuploidy.tsv`, `P6_aneuploidy_tests.tsv`,
  `P6_summary.tsv`; figures: `P6_arm_landscape`, `P6_aneuploidy`.
- **Module P7 (pairwise) — amplitude among altered samples**: splits each
  locus's between-cohort signal into a frequency axis (Fisher) and an
  amplitude axis (Mann-Whitney + rank-biserial **within the altered
  samples**), plus a high-level fraction test among altered samples and a
  freq×amplitude quadrant classification (freq_driven / amp_driven / both /
  neither). Three unit layers: 1 Mb windows and cytobands (GISTIC's
  `D.cap1.5.mat`, seg-derived fallback) and union peaks (the previously
  discarded `- CN values` rows of all_lesions, now parsed).
  Outputs: `P7_window_amplitude.tsv`, `P7_band_amplitude.tsv`,
  `P7_peak_amplitude.tsv`, `P7_summary.tsv`; figures: quadrants, volcanoes,
  per-sample amplitude strips of the strongest units.
- **Module P8 (pairwise) — multivariate cohort separation**: PERMANOVA
  (pseudo-F, R², permutation p; within-pair restricted permutation for
  paired designs), PCA coordinates, out-of-fold cross-validated AUC of an
  L2/elastic-net logistic classifier with a permutation null, and bootstrap
  stability selection. Feature sets: 1 Mb grid (default, peak-call
  independent), cytoband, union peaks (documented circularity). Requires
  scikit-learn (new dependency). Outputs: `P8_separation.tsv`,
  `P8_pca_coords.tsv`, `P8_stability.tsv`, `P8_summary.tsv`; figures: PCA,
  AUC permutation null, stability bars.
- `io_gistic.load_all_lesions(..., return_cn=True)` now exposes the
  per-sample continuous CN rows (`peak_cn`), plus new parsers
  `load_broad_by_arm` (0/1/2 and 0/1/−1 conventions auto-detected) and
  `load_arm_qv`.
- Statistical caveats for P6/P7/P8 in `docs/stats.md`; pairwise workflow
  README tables now document P1–P8 + CF + GS.

### Changed
- `scikit-learn>=1.3` added to dependencies (module P8); sklearn's
  post-1.8 `penalty` deprecation is handled (`l1_ratio` API).
- Synthetic test bundle now also writes `broad_values_by_arm.txt`,
  `amp_qv`/`del_qv` tables, differentiates all_lesions `- CN values` rows
  by direction/alteredness, and carries a broad chr2q signal; test suite
  grows to 65 tests.

### Fixed
- **Module P10 — two permutation-validity defects** (found by new
  calibration tests against exact enumeration of the margin class):
  1. `_perturb`'s batch collision filter dropped both sides of any
     colliding swap pair; with a batch far larger than the matrix nearly
     every valid swap was discarded, the Markov chain barely moved (2
     applied swaps per 4800 attempts on a 60x8 matrix) and `p_perm`
     collapsed to 1/(B+1) or 1.  Batch size is now scaled to the matrix
     (~10% of attempts apply).
  2. `Z_perm` was stored as float32: ties at the observed state rounded
     to a slightly smaller |Z| and the `>=` comparison counted them as
     non-exceedances, pushing `p_perm` to the floor for null peaks.
     Stored as float64 now.
- **Silent data-misalignment fixes**: module F G-scores are attached to
  peaks by `peak_id` (Del peaks used to receive Amp peaks' G-scores);
  module PD pairs raw/centered window statistics by window id (trust
  layers were mislabelled whenever the two variants skipped different
  windows); module P10's binary matrix selects values by sample label
  (shared-sample drops used to shift values onto the wrong samples) and
  subsets to metadata-complete samples correctly (KeyError on partial
  metadata).
- **Module P9**: per-comparison rank normalization is computed within each
  comparison's own pooled sample set and returned in a per-comparison long
  table (`P9_burden_scores.tsv` gains a `comparison` column); overlapping
  comparisons used to overwrite each other's percentiles.
- **Module PQ**: exact-zero (underflowed) p-values are kept for lambda and
  the reported `n` (they were dropped, biasing lambda downward); corrupt
  p > 1 can no longer enter the Q-Q.
- **P4/P5 null calibration** now uses pooled-sample label permutations at
  the comparison's own group sizes with the same binary pipeline as the
  main analysis (Fisher + MW, BH per direction); the earlier within-cohort
  half-split z-test null calibrated neither the sample size nor the test.
- P10 volcano ring legend no longer inverts the meaning of `burden_flag`;
  module F whiskers are drawn only on the count panel (they were Wilson
  CIs of a different quantity on the weighted panels).
- Paired comparisons: P8 no longer crashes on partially overlapping
  sample sets (incomplete pairs are dropped with a notice), and skip
  messages state explicitly that only unpaired tests are implemented.

## [0.1.0] - 2026-10-07

First open-source release. Reorganized from the internal
`GISTIC_compare` project into a distributable, reusable package.

### Added
- **ABS — HRD genomic-scar scores via vendored scarHRD** (MIT reference
  implementation, `R/scarHRD/`): per-sample HRD-LOH / TAI / LST /
  HRD-sum / ploidy-adjusted score in `ABS_hrd_scores.tsv`
  (`scripts/run_hrd_scores.py`), plus four distribution figures
  (density + ECDF per metric).  Confounding documented: HRD-sum ~
  segmentation count rho = 0.90 and ~ purity rho = -0.60 — interpret on
  the informative tier with a purity floor; between-cohort comparison
  null (MW p = 0.43-0.91).
- `gistic-compare` console command (and `python -m gistic_compare`) replacing
  the ad-hoc `run_analysis.py`, with `--outdir`, `--list-modules`, and
  module pre-checks (G/P5/GS skip with a reason instead of failing).
- **Run manifest**: every run writes `run_manifest.json` (resolved config
  snapshot, tool/dependency versions, git commit, per-module status/timings).
- **MDE / power table** (`P4_power.tsv` / `P5_power.tsv`): minimum detectable
  frequency difference per comparison × direction × baseline at
  alpha=0.05/power=0.8, plus `max_abs_diff_*` and
  `n_windows_diff_gt_mde_*` columns in the summary tables.
- **Union-peak OncoPrint** (`P1_onco_calls.tsv` + `P1_oncoprint[_<comp>]`
  figures): per-sample 0/1/2 calls across cohorts with a Δfreq strip;
  grey encodes "no significant peak called in that cohort".
- **Scheme-A labelling** for crowded scatters (`output.figure.annotate_top_k`,
  default 14): pre-specified loci/bands render as black rings instead of
  text labels so they cannot crowd out the strongest generic signals, which
  get collision-safe labels. Applied to the CF volcano / required-effect
  scatter, GS volcano and genome Manhattan panels, P1 unique-peak genome,
  and the P2 McNemar volcano.
- Preset system: paper-/study-specific values moved out of code into YAML
  (`presets/reference_panels/*.yaml` for module G, `presets/loci_sets/*.yaml`
  for P4/P5/GS/CF locus annotations). Unconfigured = fully generic runs.
- Bundled hg19 cytoband + refgene TSVs (`src/gistic_compare/data/refs/`) with
  `scripts/build_refs.py` to regenerate them from GISTIC's `hg19.mat`;
  user-supplied `.mat`/TSV paths still take precedence.
- Test suite (37 tests): parsers, config topology, interval matching, stats
  helpers, and both workflows end-to-end on a synthetic GISTIC bundle
  (`tests/synth.py`, incl. a MATLAB v7.3 `D.cap1.5.mat`).
- Packaging metadata, MIT license, bilingual README (EN/中文), CONTRIBUTING,
  CITATION.cff, changelog, GitHub Actions CI.
- `plotting` split from one 3,000-line module into a per-analysis-module
  package.

### Changed
- Results are written relative to the config file (or `--outdir`), never
  inside the installed package; `output.tables`/`output.figures` are honored.
- Cohort colors are fully config-driven (semantic defaults removed).
- Reference tables resolve relative to the config file; bundled tables are
  the fallback instead of machine-specific absolute paths.
- hg19 chromosome constants moved to `genomes.py` (single source of truth);
  unsupported genome builds raise a clear error.
