# Changelog

All notable changes to this project are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

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
