"""Tests for module P9 (three-level SCNA burden + wGII)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from gistic_compare import io_gistic, pairwise_burden


def test_threshold_c_boundaries():
    f = pairwise_burden._threshold_c
    vals = np.array([-1.5, -1.0, -0.5, -0.25, -0.2, 0.0, 0.2, 0.25, 0.5,
                     1.0, 1.5])
    expect = np.array([-2, -1, -1, -1, 0, 0, 0, 1, 1, 2, 2])
    assert np.array_equal(f(vals, 0.25, 1.0), expect.astype(float))


def test_p9_focal_merge_and_exclusion():
    """Adjacent touching same-C segments merge; gapped ones stay separate;
    broad-covered intervals drop."""
    segs = pd.DataFrame(dict(
        sample=["s"] * 4, chrom=["chr1"] * 4,
        start=[0, 1_000_000, 5_000_000, 9_000_000],
        end=[1_000_000, 2_000_000, 6_000_000, 10_000_000],
        value=[0.5, 0.6, 0.4, 0.3]))
    ev = pairwise_burden._focal_from_segments(segs, 0.25, 1.0)
    # segments 1+2 touch and share C=1 -> merged [0, 2Mb]; segments 3/4 are
    # separated by gaps -> their own events
    assert len(ev) == 3
    assert (ev[0]["start"], ev[0]["end"]) == (0, 2_000_000)
    # a broad event covering 0-8 Mb excludes the first two intervals
    # (overlap 100% > 0.5) but keeps the third (overlap 0%)
    kept = pairwise_burden._drop_broad_covered(
        ev, [dict(sample="s", chrom="chr1", direction="Amp", start=0,
                  end=8_000_000)], 0.5)
    assert len(kept) == 1 and kept[0]["start"] == 9_000_000


def test_p9_chrom_merge(synth_root):
    """Synth cohort A: whole chr2 at +0.5 -> both arms broad, same value ->
    merged into ONE chromosome-level event (not two arm events)."""
    seg = io_gistic.load_seg(synth_root / "A" / "ca.tumors.seg", "a")
    arms = pd.DataFrame([
        dict(arm="2p", chrom="chr2", start=0, end=90_000_000),
        dict(arm="2q", chrom="chr2", start=90_000_000, end=243_199_373),
    ])
    arm_ev = pairwise_burden._arm_events(seg, arms, 0.25, 0.70)
    assert set(arm_ev["arm"]) == {"2p", "2q"}       # both arms, 8 samples
    assert len(arm_ev) == 16
    arm_out, chrom_ev = pairwise_burden._merge_chromosome_events(
        arm_ev, 0.10)
    assert len(chrom_ev) == 8                       # one event per sample
    assert len(arm_out) == 0                        # arm events consumed
    assert set(chrom_ev["arms"]) == {"2p+2q"}


def test_p9_wgii(synth_root):
    seg = io_gistic.load_seg(synth_root / "A" / "ca.tumors.seg", "a")
    wg = pairwise_burden._wgii_fga(seg)
    assert wg["wGII"].between(0, 1).all()
    assert wg["fga"].between(0, 1).all()
    # ploidy near 2 (most of the genome at seg_mean ~ 0)
    assert wg["ploidy"].between(1.5, 2.5).all()
    # cohort A carries chr2 at CN = 2*2^0.5 -> |CN - 2| = 0.83 > 0.3 for the
    # whole chromosome: wGII must be at least 1/22 for every A sample
    assert (wg["wGII"] >= 1 / 22 - 1e-9).all()


def test_p9_end_to_end(pairwise_config, tmp_path):
    from gistic_compare.cli import main
    out = tmp_path / "out_p9"
    assert main(["--config", str(pairwise_config), "--modules", "P9",
                 "--outdir", str(out), "--skip-figures"]) == 0
    tables = out / "test_pair" / "tables"
    sc = pd.read_csv(tables / "P9_burden_scores.tsv", sep="\t")
    assert {"S_focal", "S_arm", "S_chrom", "S_focal_star", "S_arm_star",
            "S_chrom_star", "S_total", "wGII", "ploidy"} <= set(sc.columns)
    # rank-normalization properties (pooled over both cohorts, n=16)
    n = len(sc)
    assert sc["S_focal_star"].mean() == pytest.approx((n + 1) / (2 * n),
                                                      abs=1e-6)
    assert sc["S_total"].between(0, 3).all()
    # synth: A carries the chr2 chromosome-level event, B does not
    a = sc[sc["cohort"] == "a"]
    b = sc[sc["cohort"] == "b"]
    assert (a["n_chrom_events"] == 1).all()
    assert (b["n_chrom_events"] == 0).all()
    assert (a["S_chrom"] >= 1).all() and (b["S_chrom"] == 0).all()
    # focal signal on chr1 5-10 Mb present in both cohorts (D matrix)
    assert (sc["S_focal"] >= 1).all()
    # cohort effect points at a
    tests = pd.read_csv(tables / "P9_burden_tests.tsv", sep="\t")
    st = tests[tests["metric"] == "S_total"].iloc[0]
    assert st["cohort_higher"] == "a"
    assert st["beta"] < 0          # x = is_b; A higher -> negative beta
    assert np.isfinite(st["p_ols"]) and 0 <= st["p_ols"] <= 1
    # cross-check correlations are finite and in range (NaN allowed when an
    # input is constant in the tiny synthetic data)
    cross = pd.read_csv(tables / "P9_crosscheck.tsv", sep="\t")
    assert cross["spearman_rho"].dropna().between(-1, 1).all()


def test_p9_end_to_end_figures(pairwise_config, tmp_path):
    from gistic_compare.cli import main
    out = tmp_path / "out_p9f"
    assert main(["--config", str(pairwise_config), "--modules", "P9",
                 "--outdir", str(out)]) == 0
    figures = out / "test_pair" / "figures"
    for name in ("P9_burden_dists.png", "P9_effect_forest.png",
                 "P9_level_decomposition.png"):
        assert (figures / name).exists()
