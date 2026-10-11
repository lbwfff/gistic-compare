"""Tests for the P6 (arm landscape), P7 (amplitude) and P8 (multivariate)
modules: parsers, unit statistics, and end-to-end runs on the synthetic
GISTIC bundle."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from gistic_compare import (
    io_gistic,
    pairwise_amplitude,
    pairwise_arm,
    pairwise_multivariate,
)

# ---------------------------------------------------------------------------
# parsers / io
# ---------------------------------------------------------------------------

def test_load_broad_by_arm(synth_root):
    b = io_gistic.load_broad_by_arm(
        synth_root / "A" / "ca.broad_values_by_arm.txt", "a")
    # A: 2q amp in 8 samples + 1q amp in 5 = 13 events
    assert len(b) == 13
    assert set(b["type"]) == {"Amp"}
    assert set(b["arm"]) == {"1q", "2q"}
    two_q = b[b["arm"] == "2q"]
    assert len(two_q) == 8


def test_load_broad_by_arm_sign_convention(tmp_path):
    """Negative values select the 1/-1 convention (1=amp, -1=del)."""
    p = tmp_path / "broad.tsv"
    p.write_text("sample\t1p\t1q\t2p\n"
                 "S1\t0\t1\t0\n"
                 "S2\t0\t-1\t0\n")
    b = io_gistic.load_broad_by_arm(p, "x")
    got = {(r["sample"], r["arm"]): r["type"] for _, r in b.iterrows()}
    assert got[("S1", "1q")] == "Amp"
    assert got[("S2", "1q")] == "Del"


def test_load_all_lesions_cn_rows(synth_root):
    peaks, calls, cn = io_gistic.load_all_lesions(
        synth_root / "A" / "ca.all_lesions.conf_90.txt", "a", return_cn=True)
    assert len(cn) > 0
    assert set(cn.columns) == {"peak_id", "sample", "cn"}
    assert set(cn["peak_id"]) <= set(peaks["peak_id"])
    # altered samples carry the direction-specific value (0.6 / -0.6)
    vals = set(np.round(cn["cn"].to_numpy(float), 3))
    assert vals <= {0.0, 0.6, -0.6}


def test_load_cohort_bundle_includes_p6_p7_inputs(synth_root,
                                                  benchmark_config):
    from gistic_compare.config import load_config
    cfg = load_config(benchmark_config)
    d = io_gistic.load_cohort(cfg, "a")
    assert d["broad"] is not None and len(d["broad"]) == 13
    assert set(d["arm_qv"]) == {"Amp", "Del"}
    assert "peak_cn" in d


# ---------------------------------------------------------------------------
# P6 — derived arm calls + end-to-end
# ---------------------------------------------------------------------------

def _arms_for_test():
    return pd.DataFrame([
        dict(arm="1p", chrom="chr1", start=0, end=120_000_000),
        dict(arm="1q", chrom="chr1", start=120_000_000, end=249_250_621),
        dict(arm="2p", chrom="chr2", start=0, end=90_000_000),
        dict(arm="2q", chrom="chr2", start=90_000_000, end=243_199_373),
    ])


def test_p6_broad_from_seg_detects_broad_amp(synth_root):
    seg = io_gistic.load_seg(synth_root / "A" / "ca.tumors.seg", "a")
    calls = pairwise_arm._broad_from_seg(seg, _arms_for_test(), "a",
                                         0.30, -0.30, 0.50)
    # whole chr2 at seg_mean 0.5 -> 2q broad amp for every A sample
    got = calls[(calls["arm"] == "2q") & (calls["type"] == "Amp")]
    assert len(got) == 8
    assert set(got["sample"]) == set(seg["sample"].unique())


def test_p6_end_to_end(pairwise_config, tmp_path):
    from gistic_compare.cli import main
    out = tmp_path / "out_p6"
    assert main(["--config", str(pairwise_config), "--modules", "P6",
                 "--outdir", str(out), "--skip-figures"]) == 0
    tables = out / "test_pair" / "tables"
    diff = pd.read_csv(tables / "P6_arm_diff.tsv", sep="\t")
    assert (diff["source"] == "gistic").all()
    # synth broad file: 2q amp in 8/8 A vs 0/8 B -> significant
    hit = diff[(diff["arm"] == "2q")]
    assert (hit["freq_amp_a"] == 1.0).all()
    assert (hit["freq_amp_b"] == 0.0).all()
    assert (hit["fdr_amp"] < 0.05).all()
    aneu = pd.read_csv(tables / "P6_aneuploidy.tsv", sep="\t")
    assert {"n_broad_amp", "n_broad_del", "broad_frac"} <= set(aneu.columns)
    assert aneu["sample"].nunique() == 16
    summary = pd.read_csv(tables / "P6_summary.tsv", sep="\t")
    assert (summary["n_sig_amp"] >= 1).all()


# ---------------------------------------------------------------------------
# P7 — statistics + end-to-end
# ---------------------------------------------------------------------------

def test_p7_direction_stats():
    va = np.array([0.0, 0.5, 0.6, 1.0, 0.2, 0.0])       # 4 altered (>=0.1)
    vb = np.array([0.0, 0.0, 0.15, 0.0, 0.0, 0.12])     # 2 altered, weak
    st = pairwise_amplitude._direction_stats(va, vb, 0.1, 0.9, True, "amp")
    assert st["n_alt_amp_a"] == 4 and st["n_alt_amp_b"] == 2
    assert st["freq_amp_a"] == pytest.approx(4 / 6)
    assert st["delta_amp_amp"] > 0      # A alters more strongly
    assert np.isfinite(st["p_amp_amp"]) and np.isfinite(st["r_rb_amp"])
    assert st["p_high_amp"] == pytest.approx(1.0)  # no high-level on either side


def test_p7_end_to_end(pairwise_config, tmp_path):
    from gistic_compare.cli import main
    out = tmp_path / "out_p7"
    assert main(["--config", str(pairwise_config), "--modules", "P7",
                 "--outdir", str(out), "--skip-figures"]) == 0
    tables = out / "test_pair" / "tables"
    win = pd.read_csv(tables / "P7_window_amplitude.tsv", sep="\t")
    # the synthetic chr1 5-10 Mb window: altered in A (0.4) never in B
    hit = win[(win["chrom"] == "chr1") & (win["start"] == 5_000_000)]
    assert len(hit) == 1
    h = hit.iloc[0]
    assert h["freq_amp_a"] == pytest.approx(1.0)
    assert h["freq_amp_b"] == pytest.approx(0.0)
    assert h["n_alt_amp_a"] == 8 and h["n_alt_amp_b"] == 0
    assert h["class_amp"] in ("freq_driven", "both", "neither")
    band = pd.read_csv(tables / "P7_band_amplitude.tsv", sep="\t")
    assert len(band) > 0
    peak = pd.read_csv(tables / "P7_peak_amplitude.tsv", sep="\t")
    assert {"peak_key", "comparison", "class_amp"} <= set(peak.columns)
    both = peak[peak["present_amp"] == "both"]
    assert len(both) > 0
    summary = pd.read_csv(tables / "P7_summary.tsv", sep="\t")
    assert set(summary["layer"]) == {"window", "band", "peak"}


# ---------------------------------------------------------------------------
# P8 — PERMANOVA unit test + end-to-end
# ---------------------------------------------------------------------------

def test_p8_permanova_null_and_signal():
    from scipy.spatial.distance import pdist, squareform
    rng = np.random.default_rng(7)
    # two well-separated clusters -> large F, minimal p
    a = rng.normal(0, 0.2, (20, 3))
    b = rng.normal(10, 0.2, (20, 3))
    X = np.vstack([a, b])
    labels = np.array([0] * 20 + [1] * 20)
    D2 = squareform(pdist(X)) ** 2
    f, r2, p = pairwise_multivariate._permanova(
        D2, labels, 199, np.random.default_rng(1), None)
    assert f > 50 and r2 > 0.8 and p <= 1 / 200
    # one homogeneous cloud, random labels -> small R2, large p
    X0 = rng.normal(0, 1.0, (40, 3))
    D20 = squareform(pdist(X0)) ** 2
    f0, r20, p0 = pairwise_multivariate._permanova(
        D20, labels, 199, np.random.default_rng(2), None)
    assert r20 < 0.15 and p0 > 0.05


def test_p8_end_to_end(pairwise_config, tmp_path):
    from gistic_compare.cli import main
    out = tmp_path / "out_p8"
    assert main(["--config", str(pairwise_config), "--modules", "P8",
                 "--outdir", str(out), "--skip-figures"]) == 0
    tables = out / "test_pair" / "tables"
    sep = pd.read_csv(tables / "P8_separation.tsv", sep="\t")
    assert (sep["status"] == "ok").all()
    assert (sep["n_a"] == 8).all() and (sep["n_b"] == 8).all()
    assert sep["permanova_p"].between(0, 1).all()
    assert not (tables / "P8_skipped.txt").exists()
    pca = pd.read_csv(tables / "P8_pca_coords.tsv", sep="\t")
    assert len(pca) == 16
    stab = pd.read_csv(tables / "P8_stability.tsv", sep="\t")
    assert {"feature", "selection_freq", "mean_coef"} <= set(stab.columns)


def test_p8_peaks_features_circular_flag(pairwise_config, tmp_path):
    import yaml
    cfg = yaml.safe_load(open(pairwise_config))
    cfg["analysis"]["multivariate"] = {"features": "peaks"}
    path = pairwise_config.parent / "config_p8peaks.yaml"
    path.write_text(yaml.safe_dump(cfg))
    from gistic_compare.cli import main
    out = tmp_path / "out_p8p"
    assert main(["--config", str(path), "--modules", "P8",
                 "--outdir", str(out), "--skip-figures"]) == 0
    sep = pd.read_csv(out / "test_pair" / "tables" / "P8_separation.tsv",
                      sep="\t")
    assert (sep["feature_set"] == "peaks").all()
