# Figure conventions

All figures are drawn through `src/gistic_compare/theme.py`; the per-module
code in `src/gistic_compare/plotting/` only supplies content.

## Canvas & export

- Double-column width **7.0 in** is a hard cap (`theme.capped_width`); more
  panels wrap via `theme.grid_shape` instead of growing the canvas.
- Fonts: Arial (fallback Helvetica / DejaVu Sans). Axis labels 8 pt, ticks
  7 pt, legend 6.5 pt, annotations 6 pt.
- Vector fonts embedded: `pdf.fonttype = 42` — PDFs carry no Type 3 fonts and
  text stays editable in Illustrator.
- Every figure exports as **PDF (vector) + PNG (600 dpi)**.
- Multi-panel lettering (a/b/c) is off by default; enable with
  `output.figure.panel_labels: true`.

## Colour semantics

- Qualitative palette: **Okabe-Ito** (colour-blind safe).
- `Amp = #D55E00` (vermillion), `Del = #0072B2` (blue) everywhere.
- Cohort colours come from `output.figure.cohort_colors` in your config; any
  cohort without an explicit colour gets the tab10 entry **farthest from every
  colour already in use** (never tab10 orange next to Okabe orange).
- Δ / direction always uses **cohort colours**, never the Amp/Del colours —
  "red = amplification" and "red = cohort A higher" must not co-exist.
- Correlation heatmaps switch to `RdBu_r` (−1..1) automatically when any
  negative correlation is present.

## Statistics on figures

- Annotation format: `r = 0.49 (95% CI 0.10–0.75, n = 13), p = 0.012`;
  p < 0.001 renders as `p < 0.001`.
- Multiple-testing annotations use BH-adjusted q (`theme.bh_adjust`).
- Proportions carry **Wilson 95% CIs** (`theme.wilson_ci` /
  `theme.bar_with_ci`), not normal approximations.
- Small-n distributions overlay raw points (`theme.strip_overlay`).
- Sample sizes are annotated (`theme.annot_n`) from counts recorded at load
  time (`cfg["_n_samples"]`).
- Paired/unpaired status of a comparison is printed in the figure footnote.

## Point labels on crowded scatters (scheme A)

Label placement follows one rule everywhere:

- **Black rings** mark the pre-specified loci/bands from your `loci_set`
  preset (`ring_max_per_label` rings only the strongest row per label so
  thick clouds don't turn into rings on rings).
- **Text labels** go to the strongest *generic* points (ranked by the panel's
  evidence axis, e.g. −log10 FDR or |Δfreq|), capped by the budget
  `output.figure.annotate_top_k` (default 14), placed with screen-space
  collision avoidance.

Applied to: CF E-value volcano + required-effect scatter, GS volcano (band /
window), GS genome Manhattan panels, P1 unique-peak genome, P2 McNemar
volcano, P7 quadrants/volcanoes (labels only on units significant on at
least one axis), P10 volcano, P8 stability panels, and the peak-gene labels
of `P11_focal_landscape` (variant — see the landscape section below).
Module ABS adds
`ABS_wgd_freq` (WGD+ fraction bars with Wilson CIs),
`ABS_metric_dists` (purity/ploidy/LOH/CN-LOH per cohort), and
`ABS_purity_adjusted` (crude vs purity-adjusted per-peak OR).

- No descriptive titles inside figures — keep captions in your manuscript.

## OncoPrint (P1)

Rows are union peak keys (moderate-match chains), columns are samples grouped
by cohort. Cells: 0 = white, 1/2 = mid/deep amplitude (vermillion = Amp, blue
= Del), **grey = the cohort has no member peak for that key** (GISTIC did not
call a peak — distinct from "called, no alteration"). The right-hand strip
shows Δfreq (A − B) in cohort colours, treating a missing peak as frequency 0
(stated in the figure footnote). Cohort labels carry n.

## CNV-landscape panels (P9 event frequencies, P11 focal landscape)

These two figure families share one convention that is the **reverse of the
default colour rule**: cohort is encoded as a **shade of the Amp/Del hue**,
not by the config cohort colour.  Rationale — both families separate the
directions *positionally* (up/down bars, mirrored panels), so the panel
already says which side is Amp; that frees the colour for the direction hue
(vermillion/blue) and the cohort for its lightness (`_core._blend`: full hue
= cohort A, ~half-blended = cohort B).  Two distinct cohort
colours would clash with the Amp/Del hue inside the same panel.  Applies to
`P9_arm_events`, `P9_chrom_events` and `P11_focal_landscape`.

### P9 event frequencies (`P9_arm_events`, `P9_chrom_events`)

One chart per comparison.  For every arm (or chromosome) an up (Amp) and a
down (Del) pair of bars — one bar per cohort — sits on a single mirrored
frequency scale.

- **Event definition** follows `analysis.burden` exactly (`arm_min_frac`,
  `arm_thr`, `arm_match_diff`), i.e. the same arm and chromosome events the
  burden scores are built from, so the figure and `P9_burden_dists` never
  disagree.  Arm frequencies count every arm event *before* chromosome
  merging (merging is a burden-accounting step only).
- **Cohort-difference marker**: a green dot on the baseline, upper = Amp,
  lower = Del.  Two-sided Fisher per unit; BH-FDR is computed **within each
  (comparison, level, direction) cell**, the same convention as P4/P6, not
  across the whole table.
- **No Wilson CIs, on purpose.**  44 arms × 2 cohorts × 2 directions would
  turn them into ink; the counts live in `P9_event_freq.tsv`.
- The difference test is the **unpaired Fisher**; on a `kind=paired`
  comparison it must be read as conservative (it does not exploit pairing) —
  the same stance as P6/P7/GS, while paired inference is on the roadmap
  (`docs/stats.md` §E1).

### P11 focal landscape (`P11_focal_landscape`)

Each cohort's `scores.gistic` `-log10(q)` track drawn genome-wide,
Amplifications above / Deletions below a `q` threshold line from
`analysis.focal.q_threshold` (default: the first entry of
`analysis.q_thresholds`).

- **Non-linear y axis (deliberate).**  Linear in `-log10(q)` down to
  q = 10⁻⁵, symlog-compressed beyond (`plotting/pairwise_p11.py`,
  `_SYMLOG_LINTHRESH`).  A single extreme deletion (10⁻²⁵ or deeper) must not
  hide the 10⁻²…10⁻⁸ peaks that carry the biology; the compression point is a
  hardcoded constant, so a dataset whose signal sits deeper than 10⁻⁵ is
  compressed harder.
- Equal-height panels with **independent y ranges**; both start at the
  baseline `q = 1` (the tick labels show the q magnitude, the threshold line
  is the only significance reference).
- The **lighter (last) cohort is drawn as line + fill, the darker as line
  only**, so overlapping tracks stay readable.
- **Descriptive overlay only.**  q-values are computed within each cohort and
  GISTIC's q has no between-cohort null: the figure carries no cohort-difference
  marker or test (same caveat as P6's arm q-values and module GS).  Cohort
  sizes also differ, and q is sample-size sensitive — read it next to P1/P7/GS
  /P10 for any actual comparison.
- **Peak-gene labels** pick the strongest union peaks of *both* cohorts with
  q below the threshold, one gene per peak, **protein-coding names preferred**
  (miRNA / `LOC*` / `linc*` are skipped unless the peak contains nothing
  else), capped by `annotate_top_k` with the same screen-space collision
  avoidance as scheme A.  Source is `P11_focal_peaks.tsv`.
