# gistic-compare

Compare **GISTIC 2.0** results across cohorts and processing variants —
peak-level overlap, sample-level concordance, regional between-cohort
differences, and publication-ready figures. One CLI, two workflows, fully
config-driven.

[中文文档](README.zh-CN.md)

## Why

GISTIC runs are easy to produce and hard to compare: every run reports its own
peak list, its own q-values, and its own per-sample calls. `gistic-compare`
answers the questions you actually have when you put two (or more) GISTIC
outputs side by side:

- *Do the two runs find the same peaks?* (A, P1)
- *Do the same samples get the same calls?* (C, P2)
- *Is a region's copy-number distribution actually different between cohorts?* (P4/P5)
- *How much of the observed difference could confounders explain?* (CF)

## Installation

```bash
pip install gistic-compare          # from PyPI (once published)
# or, from a clone of this repository:
pip install .
```

Requires Python ≥ 3.10 and GISTIC 2.0 outputs (hg19) per cohort.

Verify:

```bash
gistic-compare --list-modules
gistic-compare --version
```

## Quick start

1. **Point a config at your GISTIC output directories** (`examples/config.pairwise.yaml`
   is a commented template):

   ```yaml
   workflow: pairwise            # or: benchmark
   data_root: /path/to/data      # dir containing your cohort folders
   run_id: my_first_run
   comparison_mode: experiment_pairs
   cohorts:
     exp1: {label: "Exp 1", path: exp1_dir, prefix: exp1}
     exp2: {label: "Exp 2", path: exp2_dir, prefix: exp2}
   analysis: {q_thresholds: [0.25, 0.05]}
   ```

2. **Run it:**

   ```bash
   gistic-compare --config my_config.yaml --all
   ```

3. **Read the outputs** in `results/my_first_run/{tables,figures}` next to
   your config (or wherever `--outdir` points).

Each cohort folder must contain the GISTIC files named `<prefix>.<suffix>`:

| file | used by |
|---|---|
| `<prefix>.all_lesions.conf_90.txt` | A, C, P1–P3 |
| `<prefix>.amp_genes.conf_90.txt` / `.del_genes.conf_90.txt` | D, G, P3 |
| `<prefix>.scores.gistic` | B |
| `<prefix>.tumors.seg` | E, P4 |
| `<prefix>.sample_seg_counts.txt` | E |
| `<prefix>.regions_track.conf_90.bed` | A (cross-check) |
| `<prefix>.D.cap1.5.mat` | P5, GS (optional) |

## The two workflows

### benchmark (modules A–G) — one control, many experiments

Star topology: every experiment is compared against a shared control cohort
(e.g. TCGA). Answers *"how does my processed data compare to a reference?"*
and *"which signals are specific to my experiment?"*.

| ID | content |
|---|---|
| A | peak overlap: Jaccard, q/frequency correlations + bootstrap CI, permutation null, shared-peak membership |
| B | scores.gistic → 1 Mb grid resampling, genome-track correlations |
| C | per-sample 0/1/2 call concordance / κ / McNemar (paired only) + frequency correlations |
| D | gene / cytoband set Jaccard, driver-gene presence |
| E | per-sample FGA and segment QC |
| F | experiment-specific signal fraction (vs control), 4 weightings |
| G | benchmark against a **reference panel** (numbers externalized to `presets/reference_panels/*.yaml`; disabled unless configured) |

### pairwise (modules P1–P5, CF, GS) — cohort vs cohort, no control

| ID | question | basis |
|---|---|---|
| P1 | do the peak lists agree? (matched pairs, Δq/Δfreq, unique peaks, tiered Jaccard, UpSet, per-sample **OncoPrint** of union peaks) | GISTIC peaks |
| P2 | do the same samples get the same calls? (auto-skips when unpaired) | per-peak calls |
| P3 | per-peak logistic `altered ~ cohort + covariates` (needs a sample-metadata table) | per-peak calls |
| **P4** | is the CN distribution different in this region? (1 Mb windows, Fisher + MW, null calibration, **minimum-detectable-difference table**) | `tumors.seg` |
| P5 | same as P4 on GISTIC's own marker matrix | `D.cap1.5.mat` |
| CF | how large would a confounder have to be to explain the P4/P5 findings? (bounds, not adjustment) | P4/P5 tables |
| GS | which windows/bands differ in **amplitude** (G-score) rather than frequency? | `D.cap1.5.mat` |

**P1–P3 answer detection/consistency** — a locus must first be a significant
peak in some cohort to enter them. **P4/P5 answer group difference** and do not
depend on peak calling at all. See `docs/` for statistical caveats (per-window
FDR is exploratory; the pre-specified cytoband tests and preset loci are the
confirmatory tier).

## Study-specific values live in presets, not code

Anything tied to one paper or cohort set is a YAML preset:

```yaml
analysis:
  regional:
    loci_set: presets/loci_sets/we_ea.yaml   # named loci + bookmarked bands
paper_benchmark:
  enabled: true
  panel: presets/reference_panels/we_lung.yaml   # module G reference values
```

Copy an example from `presets/`, edit the values for your study, and point your
config at it. With no preset configured the package runs fully generically:
no locus annotations, no paper figures, no hard-coded assumptions.

## Reference tables

Cytoband and refgene tables (hg19) are bundled at
`src/gistic_compare/data/refs/*.tsv` and used automatically when your config
does not supply `analysis.regional.cytoband_file` / `refgene_file`. To
regenerate them from GISTIC's own `hg19.mat`:

```bash
python scripts/build_refs.py --mat /path/to/hg19.mat
```

You can also pass plain TSVs with the same columns, or GISTIC `.mat` files,
via the config paths. P5/GS additionally need each cohort's
`<prefix>.D.cap1.5.mat` (produced by your GISTIC run) — P5/GS are skipped with
a notice when any cohort lacks it.

## Output

```
results/<run_id>/
├── tables/    # one TSV per module, 00_cohort_overview.tsv always written
└── figures/   # PDF (vector, Type-42 fonts) + PNG (600 dpi) per table
```

Every run also writes **`run_manifest.json`** next to the tables: the resolved
config snapshot, package/dependency versions, git commit (when run from a
checkout), and per-module status/timings — any result table can be traced back
to exactly what produced it.

Figures follow a publication theme (Okabe-Ito palette, Arial, capped
double-column width, Wilson CIs, BH-FDR annotations) — see `docs/figures.md`.
Crowded scatters follow one labelling rule: black rings mark the
pre-specified loci/bands from your `loci_set` preset, text labels go to the
strongest generic points, capped by `output.figure.annotate_top_k` (default 14).
The P1 OncoPrint shows per-sample 0/1/2 calls on union peaks (grey = that
cohort did not call a peak there) with the cohort-coloured Δ-frequency strip
alongside.

## Command line

```
gistic-compare --config CONFIG [--workflow {benchmark,pairwise}]
               [--modules A,B,...] [--outdir DIR] [--skip-figures]
               [--fail-fast] [--list-modules] [--version]
```

Modules run independently; one failing module does not abort the rest
(unless `--fail-fast`). Exit code is 0 when nothing failed.

## Development

```bash
git clone https://github.com/lbwfff/gistic-compare
cd gistic-compare
pip install -e ".[dev]"
pytest          # 27 tests, incl. two end-to-end runs on synthetic GISTIC data
ruff check src tests
```

The test suite generates a tiny synthetic GISTIC bundle (peaks, calls, seg
files, even a MATLAB v7.3 `D.cap1.5.mat`) and runs both workflows end-to-end
through the CLI — no external data needed.

See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## Scope & limitations

- hg19 only (chromosome set and lengths are fixed; see `genomes.py`).
- Inputs are GISTIC 2.0 outputs (conf_90 tables); the package does not run
  GISTIC itself.
- CF bounds confounding from marginal distributions only — it cannot adjust
  for confounders without sample-level metadata (that is P3).
- `results/` and any patient-level data should never be committed to this
  repository.

## License

MIT — see [LICENSE](LICENSE). If you use gistic-compare in your research,
please cite it (see [CITATION.cff](CITATION.cff)).
