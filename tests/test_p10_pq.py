"""Tests for module P10 (burden-controlled permutation) and PQ (p-value
QC)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from gistic_compare import pairwise_permutation as P10
from gistic_compare import pval_qc as PQ

# ---------------------------------------------------------------------------
# permutation machinery
# ---------------------------------------------------------------------------

def test_perturb_preserves_margins_exactly():
    """The core invariant: row and column sums must be element-wise
    identical after perturbation."""
    rng = np.random.default_rng(3)
    X = (rng.random((60, 25)) < 0.25).astype(np.int8)
    Xw = X.copy()
    swaps = P10._perturb(Xw, rng, n_attempt=200_000)
    assert swaps > 1000                    # the perturbation actually moved
    assert np.array_equal(Xw.sum(axis=0), X.sum(axis=0))
    assert np.array_equal(Xw.sum(axis=1), X.sum(axis=1))


def test_perturb_changes_matrix():
    rng = np.random.default_rng(4)
    X = (rng.random((80, 40)) < 0.3).astype(np.int8)
    Xw = X.copy()
    P10._perturb(Xw, rng, n_attempt=100_000)
    assert not np.array_equal(Xw, X)       # and it is a real reshuffling


def test_zprop_matches_formula():
    n_a, n_b = 10, 14
    col = np.zeros(n_a + n_b)
    col[:2] = 1                    # a: 2/10 altered
    col[n_a:n_a + 5] = 1           # b: 5/14 altered
    is_b = np.zeros(n_a + n_b)
    is_b[n_a:] = 1
    z = P10._zprop(col.reshape(-1, 1), is_b)[0]
    ka, kb, na, nb = 2, 5, n_a, n_b
    pooled = (ka + kb) / (na + nb)
    se = np.sqrt(pooled * (1 - pooled) * (1 / na + 1 / nb))
    assert z == pytest.approx((kb / nb - ka / na) / se, abs=1e-12)


def test_lambda_detects_inflation():
    rng = np.random.default_rng(5)
    p_null = rng.uniform(size=5000)
    p_infl = np.concatenate([rng.uniform(size=4000),
                             10 ** (-rng.uniform(0, 8, size=1000))])
    lam_null = PQ._lambda(p_null)
    lam_infl = PQ._lambda(p_infl)
    assert lam_null == pytest.approx(1.0, abs=0.1)
    assert lam_infl > 1.5


# ---------------------------------------------------------------------------
# end-to-end on the synthetic bundle
# ---------------------------------------------------------------------------

def test_p10_end_to_end(pairwise_config, tmp_path):
    from gistic_compare.cli import main
    out = tmp_path / "out_p10"
    assert main(["--config", str(pairwise_config), "--modules", "P10",
                 "--outdir", str(out), "--skip-figures"]) == 0
    tables = out / "test_pair" / "tables"
    assert not (tables / "P10_skipped.txt").exists()
    peaks = pd.read_csv(tables / "P10_peak_perm.tsv", sep="\t")
    assert {"peak_key", "p_naive", "p_perm", "q_perm",
            "burden_flag"} <= set(peaks.columns)
    assert (peaks["p_perm"] > 0).all() and (peaks["p_perm"] <= 1).all()
    assert (peaks["q_perm"] >= peaks["p_perm"] - 1e-9).all()  # BH >= raw
    # only both-side peaks survive
    assert len(peaks) >= 1
    summary = pd.read_csv(tables / "P10_summary.tsv", sep="\t")
    r = summary.iloc[0]
    assert r["n_peaks"] == len(peaks)
    assert r["n_perm"] == 10_000
    assert {"lambda_naive", "lambda_perm"} <= set(summary.columns)


def test_pq_after_p10_and_p4(pairwise_config, tmp_path):
    from gistic_compare.cli import main
    out = tmp_path / "out_pq"
    assert main(["--config", str(pairwise_config),
                 "--modules", "P4,P10,PQ", "--outdir", str(out),
                 "--skip-figures"]) == 0
    tables = out / "test_pair" / "tables"
    inf = pd.read_csv(tables / "PQ_inflation.tsv", sep="\t")
    assert {"set_id", "comparison", "n", "n_p_lt_005",
            "lambda_infl"} <= set(inf.columns)
    set_ids = set(inf["set_id"])
    assert any(s.startswith("P4_regional_diff") for s in set_ids)
    assert any(s.startswith("P10_peak_perm__p_naive") for s in set_ids)
    assert any(s.startswith("P10_peak_perm__p_perm") for s in set_ids)
    lam = inf["lambda_infl"].dropna()
    assert np.isfinite(lam).all() and (lam >= 0).all()
    assert not (tables / "PQ_skipped.txt").exists()


def test_pq_alone_skips(pairwise_config, tmp_path):
    from gistic_compare.cli import main
    out = tmp_path / "out_pq_only"
    assert main(["--config", str(pairwise_config), "--modules", "PQ",
                 "--outdir", str(out), "--skip-figures"]) == 0
    assert (out / "test_pair" / "tables" / "PQ_skipped.txt").exists()
