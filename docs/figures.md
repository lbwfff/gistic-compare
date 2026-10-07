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
- No descriptive titles inside figures — keep captions in your manuscript.

## OncoPrint (P1)

Rows are union peak keys (moderate-match chains), columns are samples grouped
by cohort. Cells: 0 = white, 1/2 = mid/deep amplitude (vermillion = Amp, blue
= Del), **grey = the cohort has no member peak for that key** (GISTIC did not
call a peak — distinct from "called, no alteration"). The right-hand strip
shows Δfreq (A − B) in cohort colours, treating a missing peak as frequency 0
(stated in the figure footnote). Cohort labels carry n.
