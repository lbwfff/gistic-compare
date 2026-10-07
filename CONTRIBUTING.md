# Contributing to gistic-compare

Thanks for considering a contribution!

## Development setup

```bash
git clone https://github.com/lbwfff/gistic-compare
cd gistic-compare
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
pre-commit install          # optional (ruff)
```

## Running the checks

```bash
pytest                      # unit + end-to-end tests (synthetic data only)
ruff check src tests
```

Both must pass before a PR can be merged. The end-to-end tests generate a
synthetic GISTIC bundle (`tests/synth.py`) — **never commit real patient
data or large analysis outputs**.

## How the code is organized

```
src/gistic_compare/
├── cli.py               # argument parsing + module orchestration
├── config.py            # YAML loading/validation, output paths, presets wiring
├── genomes.py           # hg19 chromosome constants (single source of truth)
├── presets.py           # loci-set / reference-panel preset loading
├── io_gistic.py         # all GISTIC file parsers (tidy DataFrames out)
├── intervals.py         # interval matching / Jaccard / permutation null
├── theme.py             # publication figure style + small stat helpers
├── plotting/            # one file per analysis module (split by section)
├── <module>.py          # one file per analysis module (A..G, P1..P5, CF, GS)
└── data/refs/*.tsv      # bundled hg19 cytoband + refgene tables
```

Rules of thumb:

- **No study-specific constants in code.** If a value belongs to one paper or
  cohort set, it belongs in a preset YAML (`presets/`) or the run config.
- Every module takes `(cfg, data)` and returns a dict of DataFrames; figures
  take `(res, cfg)` and only read from that dict.
- All paths in configs resolve relative to the config file; results are never
  written inside the installed package.
- New analysis modules: add `run_module_x(cfg, data)` + a `figures_x` entry
  point, register them in `cli.py`, and cover them with an end-to-end test.

## Commit / PR style

- Small, focused commits; imperative subject lines ("Add P6 module", not "added").
- PRs should state what changed, why, and include test coverage.
- Bump the version in `pyproject.toml` + `CHANGELOG.md` for user-visible changes.
