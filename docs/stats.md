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

Random within-cohort half-splits quantify the false-positive background of
the window pipeline. They do **not** simulate batch effects — a clean null
does not exclude batch/processing differences as the source of the signal.

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
