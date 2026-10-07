"""Shared fixtures: a synthetic GISTIC project in a temp dir."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))

from synth import make_cohort_dir  # noqa: E402

REPO = Path(__file__).resolve().parents[1]

COMMON_ANALYSIS = """
analysis:
  q_thresholds: [0.25, 0.05]
  paired_frac: 0.8
  peak_limit_field: wide_peak
  overlap: {min_bp: 50, reciprocal_shorter: 0.5, reciprocal_both: 0.5,
            default_tier: moderate}
  permutation: {n: 100, seed: 42}
  bootstrap: {n: 100, seed: 42, ci: 0.95}
  score_track: {bin_size: 1000000}
  regional:
    bin_size: 1000000
    amp_thr: 0.1
    del_thr: -0.1
    excl_freq_min: 0.10
    excl_freq_max: 0.05
    fdr: 0.05
    null_splits: 3
    seed: 42
  gscore: {cap: 1.5, method: normal, perm: 100, seed: 42}
  fga: {log2_threshold: 0.2, log2_threshold_sens: 0.3}
  genes: {strip_brackets: true, keep_mir: true,
          max_wide_len_mb_amp: 30, max_wide_len_mb_del: 10}
  drivers: [EGFR, MYC, TP53, FHIT, CDKN2A]
"""

FIGURE = """
output:
  figure: {style: paper, palette: okabe_ito, width_in: 7.0, dpi: 100,
           rasterize_tracks: false}
"""


@pytest.fixture()
def synth_root(tmp_path: Path) -> Path:
    make_cohort_dir(tmp_path, "A", "ca")
    make_cohort_dir(tmp_path, "B", "cb")
    return tmp_path


@pytest.fixture()
def benchmark_config(synth_root: Path) -> Path:
    """Star topology: A and B experiments vs control C."""
    make_cohort_dir(synth_root, "C", "cc")
    cfg = synth_root / "config_benchmark.yaml"
    cfg.write_text(
        "workflow: benchmark\ngenome: hg19\n"
        f"data_root: {synth_root}\nrun_id: test_bench\n"
        "comparison_mode: control_vs_all\ncontrol: c\n"
        "cohorts:\n"
        "  a: {label: CohortA, path: A, prefix: ca}\n"
        "  b: {label: CohortB, path: B, prefix: cb}\n"
        "  c: {label: ControlC, path: C, prefix: cc}\n"
        + COMMON_ANALYSIS + FIGURE)
    return cfg


@pytest.fixture()
def pairwise_config(synth_root: Path) -> Path:
    cfg = synth_root / "config_pairwise.yaml"
    cfg.write_text(
        "workflow: pairwise\ngenome: hg19\n"
        f"data_root: {synth_root}\nrun_id: test_pair\n"
        "comparison_mode: experiment_pairs\n"
        "cohorts:\n"
        "  a: {label: CohortA, path: A, prefix: ca}\n"
        "  b: {label: CohortB, path: B, prefix: cb}\n"
        "confounding:\n"
        "  enabled: true\n"
        "  sources: [P4]\n"
        "  plausible_max_binary: 0.20\n"
        "  plausible_max_age: 0.01\n"
        "  cohorts:\n"
        "    a: {sex: {male: 5, female: 3}, age_median: 60,\n"
        "        smoking: {smoker: 4, nonsmoker: 4, unknown: 0}}\n"
        "    b: {sex: {male: 4, female: 4}, age_median: 65,\n"
        "        smoking: {smoker: 2, nonsmoker: 6, unknown: 0}}\n"
        + COMMON_ANALYSIS + FIGURE)
    return cfg
