"""Regression tests for P1 review findings.

Covers: H5 (P9 per-comparison rank normalization), H4 (P10 covariate path
alignment + metadata id normalisation), PQ p=0 handling, the P4
label-permutation null, and the paired-design mechanical fixes (P8 complete
pairs, comparison-kind assignment).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from gistic_compare import (
    config as gc_config,
)
from gistic_compare import (
    pairwise_burden as P9,
)
from gistic_compare import (
    pairwise_model,
    pairwise_regional,
)
from gistic_compare import (
    pairwise_multivariate as P8,
)
from gistic_compare import (
    pairwise_permutation as P10,
)
from gistic_compare import (
    pval_qc as PQ,
)

# ---------------------------------------------------------------------------
# H5: P9 rank normalization is per comparison, never cross-overwritten
# ---------------------------------------------------------------------------

def _burden_scores():
    """Three cohorts x 3 samples with distinct S_focal values.

    Within pool(a,b) sample a1 (value 10) is rank 1/6 -> star 1/6;
    within pool(a,c) it is rank 2/6 -> star 2/6.  A shared-column
    implementation cannot represent both at once.
    """
    return pd.DataFrame(dict(
        cohort=["a"] * 3 + ["b"] * 3 + ["c"] * 3,
        sample=["a1", "a2", "a3", "b1", "b2", "b3", "c1", "c2", "c3"],
        S_focal=[10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 5.0, 70.0, 80.0],
        S_arm=[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 0.5, 6.5, 7.0],
        S_chrom=[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 0.5, 6.5, 7.0],
        wGII=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.05, 0.65, 0.7],
        fga=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.05, 0.65, 0.7],
        ploidy=[2.0] * 9,
    ))


def test_p9_ranks_computed_per_comparison():
    cfg = {"comparisons": [{"a": "a", "b": "b"}, {"a": "a", "b": "c"}]}
    pooled, tests = P9._rank_and_test(cfg, _burden_scores(), meta_df=None)
    assert len(tests) == 10                     # 2 comparisons x 5 metrics
    assert len(pooled) == 12                      # 6 pooled samples x 2 comps
    star = pooled.set_index(["comparison", "sample"])["S_focal_star"]
    assert star[("a_vs_b", "a1")] == pytest.approx(1 / 6)   # 10 lowest of 6
    assert star[("a_vs_c", "a1")] == pytest.approx(2 / 6)   # 10 is 2nd of 6
    assert star[("a_vs_b", "b3")] == pytest.approx(1.0)
    assert star[("a_vs_c", "c3")] == pytest.approx(1.0)
    # S_total is the sum of that comparison's own three percentiles
    row = pooled[(pooled["comparison"] == "a_vs_b")
                 & (pooled["sample"] == "a1")].iloc[0]
    assert row["S_total"] == pytest.approx(
        row["S_focal_star"] + row["S_arm_star"] + row["S_chrom_star"])


# ---------------------------------------------------------------------------
# H4: P10 covariate path — subset X by sample label, keep alignment
# ---------------------------------------------------------------------------

def test_p10_covariate_design_aligns_by_sample():
    rng = np.random.default_rng(0)
    samples = [f"s{i}" for i in range(12)]
    X = rng.integers(0, 2, size=(12, 3))
    is_b = np.array([0.0] * 6 + [1.0] * 6)
    # metadata for only 10 of the 12 samples (s0, s7 missing)
    meta = pd.DataFrame(dict(
        sample=[s for s in samples if s not in ("s0", "s7")],
        sex=["m", "f"] * 5))
    X2, s2, ib2, cov, note, dropped = P10._covariate_design(
        X, samples, is_b, meta)
    assert dropped == 2
    assert set(samples) - set(s2) == {"s0", "s7"}
    assert len(s2) == len(X2) == len(ib2) == cov.shape[0] == 10
    # each kept row carries its own original values (label lookup, not a
    # positional slice — the old code shifted rows past the first gap)
    for row, s in enumerate(s2):
        assert X2[row].tolist() == X[samples.index(s)].tolist()
        assert ib2[row] == is_b[samples.index(s)]


def test_p10_covariate_design_unusable_metadata_returns_unchanged():
    samples = ["a", "b", "c"]
    X = np.zeros((3, 2), dtype=np.int8)
    is_b = np.array([0.0, 0.0, 1.0])
    meta = pd.DataFrame(dict(sample=samples, sex=["m", "f", "m"]))
    X2, s2, ib2, cov, note, dropped = P10._covariate_design(
        X, samples, is_b, meta)
    assert cov is None and dropped == 0
    assert s2 == samples and X2.shape == (3, 2)


def test_metadata_id_column_normalised(tmp_path):
    p = tmp_path / "meta.csv"
    p.write_text("sample_id,sex\ns1,m\ns2,f\ns3,m\n")
    cfg = {"metadata": {"_abs": str(p), "sample_id": "sample_id",
                        "covariates": ["sex"]}}
    df = pairwise_model._load_metadata(cfg)
    assert df is not None
    assert list(df.columns) == ["sample", "sex"]   # id renamed, not a covariate


# ---------------------------------------------------------------------------
# PQ: exact-zero p-values kept in lambda; corrupt p > 1 excluded everywhere
# ---------------------------------------------------------------------------

def test_pq_lambda_keeps_underflowed_zeros():
    p = np.array([0.0] * 11 + [0.5])
    lam = PQ._lambda(p)
    assert np.isfinite(lam)          # zeros used to be dropped -> len<10 -> NaN
    assert lam > 100.0               # eleven extreme z^2 values dominate
    # lambda and the reported set share one cleaned set (n consistency)
    cleaned = PQ._clean_p(p)
    assert len(cleaned) == 12


def test_pq_clean_p_filters_corrupt_values():
    out = PQ._clean_p(np.array([0.5, 1.5, -0.1, np.nan, 0.0, 1.0]))
    assert out.tolist() == [0.5, 0.0, 1.0]


def test_pq_qq_rejects_p_above_one():
    exp, obs = PQ._qq_arrays(np.array([0.01, 0.5, 1.5]))
    assert len(exp) == len(obs) == 2   # 1.5 never enters; 0.01, 0.5 plotted
    assert (obs >= 0).all()            # no -log10 of p>1 (negative y)


# ---------------------------------------------------------------------------
# P4 null: label permutation at real group sizes, same pipeline
# ---------------------------------------------------------------------------

def _rcfg(**kw):
    base = dict(amp_thr=0.1, del_thr=-0.1, fdr=0.05, null_splits=20,
                seed=42, bin_size=1_000_000)
    base.update(kw)
    return base


def test_p4_null_label_permutation_calibrates():
    grid = pd.DataFrame(dict(chrom=["chr1"], start=[0], end=[1_000_000]))
    vals_a = np.ones((1, 8))     # every sample amplified
    vals_b = np.zeros((1, 8))    # nobody
    out = pairwise_regional._null_label_perm(vals_a, vals_b, grid, _rcfg())
    assert len(out) == 1
    r = out.iloc[0]
    assert (r["n_a"], r["n_b"]) == (8, 8)      # real group sizes preserved
    assert r["n_windows"] == 1
    # under random labels the perfect separation is essentially never
    # recovered; the old half-split z null reported inflated counts here
    assert r["sig_windows_mean"] < 1.0
    assert 0.0 <= r["sig_windows_p95"] <= r["n_windows"]


def test_p4_null_needs_two_samples_per_group():
    grid = pd.DataFrame(dict(chrom=["chr1"], start=[0], end=[1_000_000]))
    out = pairwise_regional._null_label_perm(
        np.ones((1, 1)), np.zeros((1, 5)), grid, _rcfg())
    assert out.empty


# ---------------------------------------------------------------------------
# paired mechanics: complete-pair extraction + comparison-kind assignment
# ---------------------------------------------------------------------------

def test_complete_pairs_drops_singletons():
    keep, pair_ids, pairs_idx, dropped = P8._complete_pairs(
        ["s1", "s2", "s3", "s1", "s2"])       # a: s1 s2 s3 ; b: s1 s2
    assert keep.tolist() == [True, True, False, True, True]
    assert dropped == ["s3"]
    assert pair_ids.tolist() == [0, 1, 0, 1]
    assert pairs_idx.tolist() == [[0, 2], [1, 3]]


def test_complete_pairs_empty_when_no_overlap():
    keep, pair_ids, pairs_idx, dropped = P8._complete_pairs(["a1", "a2"])
    assert not keep.any()
    assert pairs_idx.shape == (0, 2)
    assert sorted(dropped) == ["a1", "a2"]


def test_assign_comparison_kinds_boundary():
    def mk(sa, sb):
        return {"a": {"calls": pd.DataFrame(dict(sample=sa)),
                      }, "b": {"calls": pd.DataFrame(dict(sample=sb))}}
    cfg = {"comparisons": [{"a": "a", "b": "b", "kind": "unknown"}]}
    gc_config.assign_comparison_kinds(cfg, mk([f"s{i}" for i in range(5)],
                                              [f"s{i}" for i in range(5)]),
                                      paired_frac=0.8)
    assert cfg["comparisons"][0]["kind"] == "paired"      # 5/5 overlap
    cfg = {"comparisons": [{"a": "a", "b": "b", "kind": "unknown"}]}
    gc_config.assign_comparison_kinds(cfg, mk(["s0", "s1", "s2", "s3", "x"],
                                              ["s0", "s1", "s2", "s3", "y"]),
                                      paired_frac=0.8)
    assert cfg["comparisons"][0]["kind"] == "paired"      # 4/5 = 0.8
    cfg = {"comparisons": [{"a": "a", "b": "b", "kind": "unknown"}]}
    gc_config.assign_comparison_kinds(cfg, mk(["s0", "s1", "s2", "x", "y"],
                                              ["s0", "s1", "s2", "u", "v"]),
                                      paired_frac=0.8)
    assert cfg["comparisons"][0]["kind"] == "unpaired"    # 3/5 = 0.6
