# Changelog

All notable changes to this project are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added
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
