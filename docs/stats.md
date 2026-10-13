# Statistical notes and caveats

Things to keep in mind when quoting numbers from `gistic-compare`.

## Peak matching (A, P1)

Tiers (same chromosome + same type + overlap ≥ 50 bp):

| tier | rule | role |
|---|---|---|
| **moderate** | one-sided reciprocal ≥ 50% (shorter peak ≥ 50% covered) | **primary** |
| **strict** | both-direction reciprocal ≥ 50% | sensitivity (same-scale peaks) |
| **loose** | any overlap ≥ 50 bp | lower bound |

The permutation null preserves the length multiset and per-chromosome peak
counts; p = (k+1)/(n+1).

## q-value comparability (A)

q-values are compared as within-cohort rank-normalized −log10(q) (Blom z),
with Pearson r + bootstrap 95% CI. Raw q is not comparable across GISTIC runs
with different cohort sizes.

## What each tier of P4/P5 can support

- **Per-window FDR** (`*_regional_diff` `*_fdr`) — exploratory only. 1 Mb
  windows are highly correlated; BH across ~2,800 windows is far too lenient.
- **Cytoband-level tests** (`*_cytoband_diff`) — pre-specified units, no
  data-driven selection; the more defensible genome scan. MW (continuous) is
  sensitive to small systematic shifts (can reflect batch/processing
  differences); Fisher (binary frequency) is more conservative.
- **`*_paper_loci`** (loci from your preset) — pre-registered hypothesis
  tests; suitable for confirmatory citation.
- **`*_differential_regions`** — windows are *selected* on the same data that
  is then tested; p-values are optimistic. Candidates only.

## Power / minimum detectable difference (P4/P5)

`*_power.tsv` answers the question a negative result always raises: **how
large a difference could this comparison have detected?**  For each
comparison × direction (Amp/Del) × baseline frequency (5/10/20/30/50%),
it reports the minimum detectable absolute difference (MDE) at
`alpha=0.05, power=0.8` (arcsine / Cohen's-h approximation, configurable via
`analysis.regional.power_alpha` / `power_target`).

Interpretation rules:

- An observed |Δfreq| **below** the MDE means "not detectable at this sample
  size" — it does **not** mean "no difference".
- For small cohorts (n ≲ 20) the arcsine MDE is **optimistic** relative to
  exact tests (Fisher is discrete and conservative at these n); treat it as a
  lower bound on the detectable difference.
- `*_summary.tsv` additionally reports `max_abs_diff_{amp,del}` (largest
  observed window difference) and `n_windows_diff_gt_mde_{amp,del}` — windows
  whose observed difference exceeds the per-window MDE, i.e. findings that
  were **empowered to be found**.

## Null calibration (P4/P5)

Random label permutations of the pooled two-cohort samples at the
comparison's own group sizes (n_a / n_b), run through the **same** binary
pipeline as the main analysis (thresholded amp/del counts, two-sided
Fisher + Mann-Whitney, BH per direction), quantify the false-positive
background of the window pipeline.  They do **not** simulate batch effects
— a clean null does not exclude batch/processing differences as the source
of the signal.

## Module CF — bounding, not adjusting

CF only uses per-cohort **marginal** distributions (sex/age/smoking/
mutations). It computes: universal E-values (RR and OR forms), the minimum
confounder effect δ* = D/Δp needed to fully explain an observed difference,
and bias bounds under assumed effect grids. It **cannot** adjust for
confounders — that requires sample-level metadata (module P3). Binary
confounders with δ* above `plausible_max_binary` (default 0.20) are "cannot
explain" (robust). Mutations used as confounders may be mediators or
colliders — treat CF results as sensitivity analysis, not proof.

## Sample-level concordance (C, P2)

Exact McNemar (binomial) on discordant pairs, weighted Cohen's κ for 0/1/2
calls, BH across peaks. Only meaningful for **paired** comparisons
(same samples in both cohorts); unpaired comparisons are skipped
automatically.

## P3 modelling

Unit = sample × union-peak. Crude OR from `altered ~ cohort`; each covariate
is added one at a time (single-covariate adjusted OR). Sample-level
correlation within paired designs is ignored (the McNemar in P2 covers that
setting).

## Module P6 — arm-level broad landscape

When GISTIC's own `broad_values_by_arm.txt` is present in **every** cohort it
is used; otherwise arm calls are **derived from `tumors.seg`** (arm called
altered when ≥ `min_frac` of its length is covered by segments beyond
±`amp_thr`/`del_thr`, defaults 0.3/−0.3). The source is recorded in the
`source` column (`gistic` vs `derived_seg`) — the two definitions are related
but not identical, and one run always uses a single source. Arm q-values
(`amp_qv`/`del_qv`) are display-only: like GS, the difference of two
cohort-specific q-values has no null and is not tested. BH-FDR across ~39
arms per direction is far more defensible than per-window FDR, but arms are
still correlated (whole-arm events co-occur with aneuploidy); treat the
per-arm FDR as exploratory and the per-sample aneuploidy-score test as the
primary cohort-level readout.

## Module P7 — amplitude among altered samples

Alteredness at a unit is defined per sample from that sample's own value
(unit amplitude > `amp_thr`), so the altered subset does not depend on
cohort labels; the Mann-Whitney on amplitudes **within the altered subset**
is therefore an exact two-sample test — no extra permutation is needed.
Caveats: (i) with very few altered samples on one side the amplitude test is
underpowered and the table reports the counts (`n_alt_*`) so this is visible;
(ii) `p_high` compares high-level fractions **among altered samples only**
(Fisher), which is a different question from the frequency test — read the
two together, never alone; (iii) the peak layer reuses union peaks, so
loci called as peaks in only one cohort have no between-cohort statistics
(`present_*` = a_only/b_only); (iv) the freq×amplitude quadrant classes are
built from the same FDR thresholds — they are a descriptive taxonomy, not
independent tests. GISTIC's own amplitude thresholds (0.1 / 0.9, as printed
in the `Amplitude Threshold` column of all_lesions) are the defaults.

## Module P8 — multivariate separability

Three statistics per comparison: PERMANOVA pseudo-F with label permutation
(paired designs permute **within pairs** — restricted permutation), an
out-of-fold cross-validated AUC of an L2 (or elastic-net) logistic
classifier with its own permutation null, and bootstrap stability selection.
Non-negotiable caveats:

- **Separability is not biology.** Batch, platform, sample provenance and
  processing pipelines separate cohorts exactly like real CN differences do.
  The permutation p only says the observed separation exceeds what random
  labels produce *within this dataset*; it says nothing about the cause.
  Read P8 with module CF and, ideally, balance the design.
- **Feature choice matters.** The default `grid` features (1 Mb windows) are
  independent of peak calling. `peaks` features are **circular** (peaks were
  selected from these same cohorts, so separability is expected even under
  the null) — use them descriptively only.
- The OOF-AUC permutation null uses the same CV folds; with n ≲ 30 per
  cohort the AUC confidence is wide — report `auc_perm_p` alongside, never
  the bare AUC.
- Stability selection frequencies describe the fitted model, not effect
  sizes; cross-check the top features against P4/P7/GS per-locus tables.

## Module P9 — three-level SCNA burden + wGII

Per-sample burden counted at three scales (focal / arm / chromosome), each
event weighted by amplitude |C| (2 at |c|≥1, 1 at |c|≥0.25), rank-normalized
on the **pooled** two-cohort sample set, plus wGII (ploidy-adjusted,
copy-number units).

- **Read P4–P8 against this panel.** A cohort with globally higher burden
  shows higher frequency at *every* locus; a significant per-locus result
  with a burden asymmetry may just reflect the background, not the locus.
- **Rank normalization constrains the means**: by construction each
  percentile metric averages ~0.5 in the pooled set, so the OLS beta is a
  *relative positioning* effect, not an absolute-load difference. The raw
  medians (`median_a`/`median_b`) are the absolute-scale summary; report
  both.
- **The S/wGII/FGA cross-checks** (`P9_crosscheck.tsv`) should show
  moderate-to-strong positive correlation (ρ ≳ 0.4). A weak ρ means the
  definitions disagree and something needs investigation.
- The focal definition (thresholded GISTIC marker matrix, merged adjacent
  equal-C intervals, broad-overlapping intervals dropped) is a
  reproduction-specific parameterisation (Yuan et al. 2018 gives no
  executable definition) — sensitivity-analyse `focal_thr` and
  `focal_broad_overlap_max` before quoting numbers.
- With metadata, covariates are added as a single block (not one at a
  time); without metadata the comparison is unadjusted — clinical
  confounding can only be bounded (module CF).

## Module PD — baseline/batch diagnostic

PD explains *why* module PQ flags inflation and what survives it.

- **Baseline check**: each sample's baseline is the median `seg_mean` over
  neutral segments (|c| < `neutral_band`, default 0.25). A between-cohort
  baseline difference is the classic batch signature.
- **Shift structure**: the per-window median difference along the genome.
  A roughly uniform offset = pure batch/plane effect; regional structure
  (e.g. whole arms) is more likely biology or a real CN difference.
- **Centering sensitivity**: the P4-style window tests are re-run on
  within-sample centered segments (baseline subtracted). Windows
  significant in both variants are **robust** (quotable); windows
  significant only before centering are **baseline-shift candidates** and
  must be downgraded.
- Centering also removes true whole-sample scale differences (real
  aneuploidy shifts) — always read it together with module P9, which
  establishes whether event-count burden differs.

## Module P9 — parameter sensitivity

Setting `analysis.burden.sensitivity: true` runs the full 81-combo grid
over the reproduction-specific choices (`arm_min_frac` x `arm_match_diff`
x `focal_thr` x `focal_broad_overlap_max`).  The verdict per metric —
whether the beta sign is consistent in every combination — is the
robustness statement to quote; inconsistent signs mean the conclusion is
an artefact of the parameterisation, not a finding.

## Module P10 — burden-controlled permutation on peaks

The permutation randomises which sample-peak pairs are altered while
holding **both margins** of the binary matrix fixed: row sums (each
sample's total altered-peak count = global disruption) and column sums
(each peak's overall frequency). Only peak-specific *preferential
hitting* remains testable. Interpretation notes:

- Read every peak against its **naive p** (`p_naive`) and the summary
  lambda pair. A peak that is naive-significant but **not** permutation-robust
  (`burden_flag`) is what a global burden asymmetry would manufacture.
- With **no covariates** the default `stat: zprop` (two-proportion z) is
  used — asymptotically the logistic Wald Z and fully vectorised. With
  metadata, `stat: wald` fits the covariate-adjusted logistic per
  permutation (use a smaller B).
- Margin-fixed Z statistics have a slightly *larger* null variance than
  the textbook normal (hypergeometric variance correction), so
  `lambda_perm` can sit a little above 1 without indicating bias.
- The permutation controls **global burden only**. It is not a batch
  correction: platform-level differences between the cohorts need module
  PQ and, ideally, a co-normalised re-run.

## Module PQ — p-value QC

For every p-value set the module reports the genomic inflation factor
lambda = median(z^2)/0.4549 (z = Phi^-1(1 - p/2)) and a Q-Q plot.

- lambda ~ 1 is healthy; **lambda > 1.1 means systematic bias** — batch
  effects, stratification or an unadjusted confounder. The affected
  sets' "significant" units must be re-read (as hypotheses, not
  findings) until the bias is handled.
- A Q-Q lift driven by a *few* extreme points can be real signal; a lift
  of the *whole curve* (most p-values smaller than expected) is bias.
- P10's `p_naive` and `p_perm` sets are reported separately on purpose:
  their lambda gap is the burden inflation that the permutation removed.
- Small sets (n < ~50, e.g. peak layers) give noisy lambda — read them
  together with the Q-Q shape, not by lambda alone.


## Module ABS — absolute copy-number comparisons (optional extension)

Runs on ASCAT output (`R/absolute_cn.R` over marker-level logR/BAF from
the SNP6 Birdseed chain; plan `docs/plan_absolute_cn.md`).  All
comparisons are unpaired between-cohort tests on per-sample absolute-CN
summaries.

- **WGD frequency**: 2x2 Fisher (WGD+ vs WGD-) with odds ratio + CI.
  `wgd_status` is **binary by construction** — a single ASCAT snapshot
  cannot count multiple whole-genome-doubling rounds reliably (R05 §4.6);
  the continuous `wgd_frac` carries the graded information.  A table
  `mean_status` column is descriptive only, never tested.
- **WGD heuristics (declared, sensitivity-checked)**: A (primary,
  allelic): autosomal length fraction with minor CN >= 1 AND total CN >=
  `cn_thr_wgd` (3) >= `loh_frac_wgd` (0.30; sensitivity 0.2-0.5).  B
  (reference, classical ploidy rule): ploidy >= `ploidy_wgd` (3.7).
  `wgd_disagreement` marks samples where A and B differ — read those
  calls as tool-dependent (R05 §4.6 warns different tools call WGD
  differently).  Sanity range for a solid-tumour cohort: 10-60% WGD+.
- **purity / ploidy / LOH / CN-LOH**: Mann-Whitney + rank-biserial.
  Purity and ploidy are simultaneously **confounders of every
  logR-based comparison** in P1-P10 — this module is where their
  between-cohort balance is checked.
- **CN-space burden**: P9's three-level counting transplanted into
  absolute copy-number space (neutral band +/- `cn_focal_band` copies;
  broad arm event when the length-weighted mean CN exceeds
  `cn_arm_gain`/below `cn_arm_del` over > `cn_min_frac` of the arm;
  pooled rank normalisation + MW).  Thresholds are parameterised, not
  canonically fixed.
- **Purity-adjusted per-peak models**: `altered ~ cohort` vs
  `altered ~ cohort + purity` logistic per union peak (P3 machinery).
  Purity attenuates logR signal, so a cohort effect that survives
  purity adjustment is stronger evidence than the crude OR — the
  covariate only the absolute-CN stage can provide.
- **WGD multiplicity classes (strict signature)**: per segment the
  parsimonious whole-genome-doubling evidence is `nMinor >= 2` (both
  parental copies >= 2 — a single-allele gain cannot produce it);
  `wgd_class` = 0x (strict signature < `wgd_class_frac` of the
  genome) / 1x / >=2x (minor >= 3 or total >= 8 fraction >=
  `wgd_multi_frac`).  The loose heuristic (A above) is reported
  alongside: its WGD+ rate is inflated by (2,1)-type single-allele
  gains — in this cohort the two criteria even disagree on DIRECTION
  (loose: CUE 59% > CHB 50%; strict: CHB 23% > CUE 18% among
  informative samples).  `qc_tier` marks purity >= 0.99 fits as
  non-informative (purity/ploidy unidentifiable: flat "likely normal"
  genomes or balanced-elevated degenerate fits); they are excluded
  from all class statistics.
- **HRD genomic-scar scores** (scarHRD, vendored MIT reference
  implementation, Sztupinszki 2018): HRD-LOH (>=15 Mb LOH segments,
  whole-chromosome losses excluded), TAI (telomeric allelic
  imbalance), LST (large-scale transitions, >=10 Mb neighbours <3 Mb
  apart) and their sum (clinical cutpoint 42, Telli 2016 — derived in
  HGSOC, not validated in LUAD; use tertiles for stratification).
  **Purity/segmentation confounding is severe**: HRD-sum correlates
  with the segmentation count (rho = 0.90) and negatively with purity
  (rho = -0.60) — low-purity samples' segmentation noise feeds
  directly into LST/TAI counts.  Compute and interpret on the
  informative tier with a purity floor (>=0.30) only; the
  ploidy-adjusted score (`-15.5 x ploidy`) corrects the ploidy axis
  but not the segmentation-noise axis.  Full-cohort comparison CUE vs
  CHB: null at every component (MW p = 0.43-0.91).
- **Low-purity QC** (R05 §4.5): samples with purity < `min_purity`
  (0.30) are flagged; the headline WGD test is re-run without them and
  the delta is recorded in `ABS_summary.tsv`.
- **Tumour-only caveat** (`germline_mode = tumor_only`): without a
  matched/pooled germline, het SNPs come from the samples' own Birdseed
  genotype calls (`het_source = birdseed_calls`; a SNP6 array measures
  the germline genotype at most loci) or from ASCAT's built-in
  inference (`inferred`).  Somatic-LOH loci look homozygous in tumour
  calls and drop out of the het set, so **`loh_frac`/`cnloh_frac` are
  lower-bound estimates** and CN-LOH segments may be missed entirely;
  total-CN (logR) signals are unaffected.  `ABS_summary.tsv` records
  `germline_mode_a/b` + `het_source_a/b` for provenance — quote them
  with any WGD frequency.
