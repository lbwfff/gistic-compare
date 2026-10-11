"""Tests for module PD (baseline/batch diagnostic) and P9 sensitivity."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from gistic_compare import pairwise_diag

# ---------------------------------------------------------------------------
# PD units
# ---------------------------------------------------------------------------

def test_pd_baseline_recovery():
    """Baseline = median seg_mean over the neutral band (|c| < 0.25)."""
    seg = pd.DataFrame(dict(
        sample=["s1"] * 5,
        chrom=["chr1"] * 5,
        start=[0, 10_000_000, 20_000_000, 30_000_000, 40_000_000],
        end=[10_000_000, 20_000_000, 30_000_000, 40_000_000, 50_000_000],
        seg_mean=[0.10, 0.12, 0.60, -0.80, 0.08]))
    b = pairwise_diag._baseline_per_sample(seg, 0.25)
    # neutral segments: 0.10, 0.12, 0.08 -> median 0.10
    assert b["baseline"].iloc[0] == pytest.approx(0.10)
    cen = pairwise_diag._center_seg(seg, 0.25)
    # after centering, the neutral segments sit at ~0 and the events at
    # +/- their original distance from the baseline
    assert cen["seg_mean"].iloc[2] == pytest.approx(0.50)
    assert cen["seg_mean"].iloc[3] == pytest.approx(-0.90)
    assert cen["seg_mean"].iloc[0] == pytest.approx(0.0, abs=1e-12)


def test_pd_window_stats_direction():
    rng = np.random.default_rng(9)
    grid = pd.DataFrame(dict(chrom=["chr1"], start=[0], end=[1_000_000]))
    va = rng.normal(0.5, 0.05, 20)     # A clearly higher
    vb = rng.normal(0.0, 0.05, 20)
    st = pairwise_diag._window_stats(va.reshape(1, -1), vb.reshape(1, -1),
                                     grid, dict(amp_thr=0.1, del_thr=-0.1))
    r = st.iloc[0]
    assert r["p_mw"] < 1e-6
    assert r["delta_median"] > 0.4
    assert r["freq_amp_a"] == pytest.approx(1.0)
    assert r["freq_amp_b"] < 0.2


# ---------------------------------------------------------------------------
# PD end-to-end on the synthetic bundle
# ---------------------------------------------------------------------------

def test_pd_end_to_end(pairwise_config, tmp_path):
    from gistic_compare.cli import main
    out = tmp_path / "out_pd"
    assert main(["--config", str(pairwise_config), "--modules", "PD",
                 "--outdir", str(out), "--skip-figures"]) == 0
    tables = out / "test_pair" / "tables"
    lam = pd.read_csv(tables / "PD_lambda.tsv", sep="\t")
    assert {"comparison", "variant", "set_col", "n",
            "lambda_infl"} <= set(lam.columns)
    assert set(lam["variant"]) == {"raw", "centered"}
    assert (lam["n"] > 0).all()
    layers = pd.read_csv(tables / "PD_shift_regions.tsv", sep="\t")
    assert set(layers["layer"]) <= {"robust", "baseline_only",
                                    "centered_only", "ns"}
    assert len(layers) > 0
    summary = pd.read_csv(tables / "PD_summary.tsv", sep="\t")
    assert {"n_robust", "n_baseline_only", "median_baseline_a",
            "median_baseline_b"} <= set(summary.columns)


# ---------------------------------------------------------------------------
# P9 sensitivity
# ---------------------------------------------------------------------------

def test_p9_sensitivity(pairwise_config, tmp_path):
    import yaml
    cfg = yaml.safe_load(open(pairwise_config))
    cfg.setdefault("analysis", {})["burden"] = {"sensitivity": True}
    path = pairwise_config.parent / "config_p9sens.yaml"
    path.write_text(yaml.safe_dump(cfg))
    from gistic_compare.cli import main
    out = tmp_path / "out_p9s"
    assert main(["--config", str(path), "--modules", "P9",
                 "--outdir", str(out), "--skip-figures"]) == 0
    tables = out / "test_pair" / "tables"
    sens = pd.read_csv(tables / "P9_sensitivity.tsv", sep="\t")
    # 3^4 = 81 combos x 5 metrics
    assert len(sens) == 81 * 5
    assert {"arm_min_frac", "arm_match_diff", "focal_thr",
            "focal_broad_overlap_max", "metric", "beta"} <= set(sens.columns)
    verdict = pd.read_csv(tables / "P9_sensitivity_summary.tsv", sep="\t")
    assert len(verdict) == 5
    # synth: A carries the chr2 chromosome-level event -> A has higher
    # burden for every parameterisation
    st = verdict[verdict["metric"] == "S_total"].iloc[0]
    # parse the TSV round-tripped value strictly: NaN/"False"/"False-y"
    # strings must not sneak through a truthiness check
    assert str(st["sign_consistent"]).strip().lower() in ("true", "1")
    assert st["n_combos"] == 81
