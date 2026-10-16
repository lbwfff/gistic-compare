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
- *Is a region's copy-number distribution actually different between cohorts?* (P4)
- *Is one cohort globally more rearranged — focal or broad?* (P9)
- *Which differences survive burden control and quality checks?* (P10, PQ, PD)
- *How much of the observed difference could confounders explain?* (CF)
- *Do WGD frequency, purity or absolute ploidy differ between cohorts — and do
  peak effects survive purity adjustment?* (ABS — optional, needs raw SNP6 data)

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
| `<prefix>.all_lesions.conf_90.txt` | A, C, P1–P3, P7/P10 (peak CN layer) |
| `<prefix>.amp_genes.conf_90.txt` / `.del_genes.conf_90.txt` | D, G, P3 |
| `<prefix>.scores.gistic` | B |
| `<prefix>.tumors.seg` | E, P4, P6, PD (arm calls / baseline), P7–P9 (fallback) |
| `<prefix>.sample_seg_counts.txt` | E, P10 (sample QC) |
| `<prefix>.regions_track.conf_90.bed` | A (cross-check) |
| `<prefix>.D.cap1.5.mat` | P4's marker layer, GS, P7–P9 (optional) |
| `<prefix>.broad_values_by_arm.txt` | P6 (optional; else derived from seg) |
| `<prefix>.amp_qv.conf_90.txt` / `.del_qv.conf_90.txt` | P6 (optional, display-only) |

Module **ABS** is the one optional extension that needs data *beyond* GISTIC
outputs: the two harmonized absolute-CN tables per cohort.  ASCAT itself runs
**outside** this package — normally on a server — from this repository's
`R/absolute_cn.R` (marker-level logR/BAF in, `absolute_cn_<cid>.tsv` +
`absolute_seg_<cid>.tsv` out).  Run book and frozen table schema:
`docs/ascat_handoff.md`.  Without the tables ABS skips with a reason; every
other module runs on the tables above.

## The two workflows

### benchmark (modules A–G) — one control, many experiments

Star topology: every experiment is compared against a shared control cohort
(e.g. TCGA). Answers *"how does my processed data compare to a reference?"*
and *"which signals are specific to my experiment?"*.

| ID | content |
|---|---|
| A | peak overlap: Jaccard, q/frequency correlations + bootstrap CI, permutation null, shared-peak membership |
| B | scores.gistic → 1 Mb grid resampling, genome-track correlations |
| C | per-sample 0/1/2 call concordance / κ / McNemar (paired comparisons) + matched-pair frequency correlations — **the benchmark form of P2** |
| D | gene / cytoband set Jaccard, driver-gene presence |
| E | per-sample FGA and segment QC |
| F | experiment-specific signal fraction (vs control), 4 weightings |
| G | benchmark against a **reference panel** (numbers externalized to `presets/reference_panels/*.yaml`; disabled unless configured) |

### pairwise (modules P1–P4, P6–P11, PQ, PD, CF, GS, ABS) — cohort vs cohort, no control

| ID | question | basis |
|---|---|---|
| P1 | do the peak lists agree? (matched pairs, Δq/Δfreq, unique peaks, tiered Jaccard, UpSet, per-sample **OncoPrint** of union peaks) | GISTIC peaks |
| **P2** | do the same samples get the same calls? (κ / exact McNemar on **paired** pairs; matched-pair frequencies for every pair) | per-peak calls |
| P3 | per-peak logistic `altered ~ cohort + covariates` (needs a sample-metadata table) | per-peak calls |
| **P4** | is the CN distribution different in this region? (1 Mb windows, Fisher + MW, null calibration, **minimum-detectable-difference table**) | `tumors.seg` |
| P4 (marker layer) | the **same statistics on GISTIC's own marker matrix** (`P5_*` tables), run as a second layer of the module when every cohort provides `D.cap1.5.mat`; disable with `analysis.regional.marker_matrix: false` | `D.cap1.5.mat` |
| **P6** | does the **arm-level (broad) landscape** differ — the aneuploidy background? (per-arm Fisher, per-sample aneuploidy score) | `broad_values_by_arm.txt` or derived from `tumors.seg` |
| **P7** | same frequency but **stronger alteration**? (amplitude compared *among altered samples*, high-level fraction, freq×amplitude quadrants) | `D.cap1.5.mat` + all_lesions `- CN values` rows |
| **P8** | are the cohorts separable **as a whole**, and which loci drive it? (PERMANOVA + permutation, CV-AUC elastic-net, stability selection) | window/peak feature matrix |
| **P9** | is one cohort globally **more rearranged** — and at which scale (focal vs broad)? (three-level SCNA burden score + wGII) | `tumors.seg` + `D.cap1.5.mat` (or seg fallback) |
| **P10** | which peaks survive once the **global burden** is controlled? (margin-preserving permutation, naive-vs-permutation inflation) | all_lesions `- CN values` rows |
| **P11** | what does the **focal landscape** look like genome-wide? (per-cohort `-log10(q)` from `scores.gistic`, Amp up / Del down, peak genes annotated; descriptive overlay, no between-cohort test) | `scores.gistic` + `amp/del_genes.conf_90.txt` |
| **PQ** | are the p-value sets healthy? (Q-Q plots + genomic inflation λ for every scan) | saved tables of the same run |
| **PD** | **why** is lambda inflated — and what survives baseline removal? (per-sample baseline check, shift-structure track, centering sensitivity with trust layers) | `tumors.seg` |
| **ABS** | **absolute-copy-number extension** (optional): WGD frequency (Fisher+OR), purity/absolute-ploidy/LOH/CN-LOH (MW; LOH also within WGD strata), LOH x GISTIC-Del class cross (del+LOH / del-only / LOH-only / neither, copy-neutral split), P9's three-level burden in CN space, purity-adjusted per-peak logistic | `absolute_cn_<cid>.tsv` + `absolute_seg_<cid>.tsv` (ASCAT runs on a server from `R/absolute_cn.R` — `docs/ascat_handoff.md`; skips with a reason without the tables) |
| CF | how large would a confounder have to be to explain the P4/P5 findings? (bounds, not adjustment) | P4/P5 tables |
| GS | which windows/bands differ in **amplitude** (G-score) rather than frequency? | `D.cap1.5.mat` |

**P1–P3 answer detection/consistency** — a locus must first be a significant
peak in some cohort to enter them. **P4–P7 answer group difference** and do not
depend on peak calling at all (P7's peak layer is the exception). P8 is a global
multivariate summary, **P9 the burden background** every per-locus result
should be read against, **P10 the burden-controlled peak test**, and **PQ/PD
the quality gates** (inflation, batch/baseline diagnostics); **P11** is the
descriptive scan-level **focal q-value landscape** (amplifications and
deletions mirrored genome-wide) — an overlay, not a test. See `docs/` for
statistical caveats (per-window FDR is exploratory; the pre-specified cytoband
tests and preset loci are the confirmatory tier; P8 separability is not
evidence of biology).

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
via the config paths. P4's marker layer, GS, P7–P9 additionally use each
cohort's `<prefix>.D.cap1.5.mat` (produced by your GISTIC run) — the marker
layer and GS are skipped with a notice when any cohort lacks it, while
P7–P9 fall back to seg-derived window values.

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
pytest          # 100+ tests: parsing/config/intervals/stats + end-to-end runs
                # on synthetic GISTIC data + per-module tests; the ABS
                # R/ASCAT pipeline tests (auto_run) skip without Rscript+ASCAT
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
- P6 derives arm-level broad calls from `tumors.seg` when GISTIC's own
  `broad_values_by_arm.txt` is absent — the two definitions are related but
  not identical, and the source is recorded in every P6 table.
- P8 separability is a descriptive global statistic: batch, platform and
  processing differences separate cohorts just as well as biology. Read it
  together with CF.
- P10 controls the **global burden** confounder only; it is not a batch
  correction. PQ/PD exist to surface (and, for baselines, remove) exactly
  that class of artifact — run them before quoting scan results.
- An optional absolute-copy-number extension (WGD/purity/ploidy/LOH via
  ASCAT, R-backed) is implemented; ASCAT runs on a server and only its two
  harmonized tables enter this package — see `docs/ascat_handoff.md` and
  `docs/plan_absolute_cn.md`.
- CF bounds confounding from marginal distributions only — it cannot adjust
  for confounders without sample-level metadata (that is P3).
- `results/` and any patient-level data should never be committed to this
  repository.

## License

MIT — see [LICENSE](LICENSE). If you use gistic-compare in your research,
please cite it (see [CITATION.cff](CITATION.cff)).
