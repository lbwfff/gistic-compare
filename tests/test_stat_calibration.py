"""Statistical validity tests (review item A3).

Two previously untested areas:
  * module CF — E-value / RR / OR / required-effect arithmetic checked
    against hand calculations (the module had zero numeric coverage);
  * module P10 — the vectorised two-proportion z cross-checked against
    statsmodels, and the margin-preserving permutation p checked for
    calibration under the null (uniform) and power under a real signal.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from gistic_compare import confounding
from gistic_compare import pairwise_permutation as P10

# ---------------------------------------------------------------------------
# CF: hand-checked arithmetic
# ---------------------------------------------------------------------------

def test_cf_evalue_vanderweele_ding():
    # E(RR) = RR + sqrt(RR * (RR - 1)); symmetric in RR <-> 1/RR
    assert confounding._evalue(1.0) == pytest.approx(1.0)
    assert confounding._evalue(2.0) == pytest.approx(2.0 + math.sqrt(2.0))
    assert confounding._evalue(0.5) == pytest.approx(2.0 + math.sqrt(2.0))
    assert confounding._evalue(4.0) == pytest.approx(4.0 + math.sqrt(12.0))
    assert np.isnan(confounding._evalue(0.0))
    assert np.isnan(confounding._evalue(np.nan))


def test_cf_rr_or_cells_hand_values():
    # fa=0.5, fb=0.25, na=nb=8 -> Haldane cells 4.5/9 and 2.5/9
    rr, or_ = confounding._rr_or_cells(0.5, 0.25, 8, 8)
    assert rr == pytest.approx((4.5 / 9) / (2.5 / 9))          # 1.8
    assert or_ == pytest.approx((4.5 * 6.5) / (4.5 * 2.5))     # 2.6
    # zero cells are regularised, never inf/nan
    rr0, or0 = confounding._rr_or_cells(1.0, 0.0, 8, 8)
    assert np.isfinite(rr0) and np.isfinite(or0) and rr0 > 1 and or0 > 1


def _cf_cfg(tmp_path):
    """Config mirroring tests/conftest.py's confounding block."""
    return {
        "output": {"_tables": tmp_path},
        "comparisons": [{"a": "a", "b": "b"}],
        "confounding": {
            "sources": ["P4"],
            "plausible_max_binary": 0.20,
            "plausible_max_age": 0.01,
            "delta_grid": [0.02, 0.05, 0.10, 0.20],
            "age_grid": [0.001, 0.002, 0.005, 0.010],
            "cohorts": {
                "a": {"sex": {"male": 5, "female": 3}, "age_median": 60,
                      "smoking": {"smoker": 4, "nonsmoker": 4}},
                "b": {"sex": {"male": 4, "female": 4}, "age_median": 65,
                      "smoking": {"smoker": 2, "nonsmoker": 6}},
            },
        },
    }


def test_cf_envelope_exposure_gaps(tmp_path):
    cfg = _cf_cfg(tmp_path)
    env = confounding.build_envelope(cfg).set_index("confounder")
    assert env.loc["sex", "delta_p"] == pytest.approx(5 / 8 - 4 / 8)   # 0.125
    assert env.loc["age", "delta_p"] == pytest.approx(-5.0)            # years
    assert env.loc["smoking", "delta_p"] == pytest.approx(0.5 - 0.25)  # 0.25


def test_cf_required_effect_and_robustness(tmp_path):
    """delta* = D / dp exactly, and the robust flag follows its own rule."""
    cfg = _cf_cfg(tmp_path)
    env = confounding.build_envelope(cfg)
    # one strong window (diff 0.5) and one weak-but-significant (diff 0.02)
    rows = [dict(comparison="a_vs_b", level="window", chrom="chr1",
                 start=0, end=1_000_000, band="strong", n_a=8, n_b=8,
                 freq_amp_a=0.75, freq_amp_b=0.25, freq_del_a=0.0,
                 freq_del_b=0.0, diff_amp=0.5, diff_del=0.0,
                 p_amp=0.001, p_del=1.0, p_amp_fdr=0.01, p_del_fdr=1.0,
                 p_mw_fdr=0.01),
            dict(comparison="a_vs_b", level="window", chrom="chr1",
                 start=2_000_000, end=3_000_000, band="weak", n_a=8, n_b=8,
                 freq_amp_a=0.27, freq_amp_b=0.25, freq_del_a=0.0,
                 freq_del_b=0.0, diff_amp=0.02, diff_del=0.0,
                 p_amp=0.01, p_del=1.0, p_amp_fdr=0.02, p_del_fdr=1.0,
                 p_mw_fdr=0.02)]
    pd.DataFrame(rows).to_csv(tmp_path / "P4_regional_diff.tsv",
                              sep="\t", index=False)
    sens, ev = confounding._sensitivity_rows(cfg, "P4", env)
    # ev rows are emitted for both directions (Del is gated only by the
    # shared MW FDR): 2 findings x 2 directions
    assert len(ev) == 4
    for _, r in sens.iterrows():
        assert r["required_effect"] == pytest.approx(
            r["diff"] / r["delta_p"])
        per_year = r["required_unit"] == "per_year"
        max_plausible = 0.01 if per_year else 0.20
        g = 0.001 if per_year else 0.02      # first value of each effect grid
        assert bool(r["robust"]) == (abs(r["required_effect"]) > max_plausible)
        assert r[f"bias_e{g:g}"] == pytest.approx(abs(r["delta_p"]) * g)
        assert r[f"residual_e{g:g}"] == pytest.approx(
            abs(r["diff"]) - abs(r["delta_p"]) * g)
    # windows carry their paper-locus label in "band", so key on coordinates
    strong_sex = sens[(sens["start"] == 0)
                      & (sens["confounder"] == "sex")].iloc[0]
    assert strong_sex["required_effect"] == pytest.approx(0.5 / 0.125)  # 4.0
    assert bool(strong_sex["robust"]) is True
    weak_smoke = sens[(sens["start"] == 2_000_000)
                      & (sens["confounder"] == "smoking")].iloc[0]
    assert weak_smoke["required_effect"] == pytest.approx(0.02 / 0.25)  # 0.08
    assert bool(weak_smoke["robust"]) is False
    # E-values: RR form on rr, OR form on sqrt(OR) (rare-outcome conversion)
    for _, r in ev.iterrows():
        assert r["evalue_rr"] == pytest.approx(confounding._evalue(r["rr"]))
        assert r["evalue_or"] == pytest.approx(
            confounding._evalue(math.sqrt(r["or_"])))


# ---------------------------------------------------------------------------
# P10: statistic and permutation-p validity
# ---------------------------------------------------------------------------

def test_p10_zprop_matches_statsmodels():
    """Independent cross-check of the vectorised two-proportion z."""
    sm_prop = pytest.importorskip("statsmodels.stats.proportion")
    rng = np.random.default_rng(3)
    X = (rng.random((30, 6)) < 0.4).astype(np.int8)
    is_b = np.array([False] * 15 + [True] * 15)
    z = P10._zprop(X, is_b)
    for j in range(6):
        kb, ka = int(X[is_b, j].sum()), int(X[~is_b, j].sum())
        z_ref, _ = sm_prop.proportions_ztest([kb, ka], [15, 15])
        assert z[j] == pytest.approx(float(z_ref))


def _perm_p(X: np.ndarray, is_b: np.ndarray, seed: int, B: int = 99):
    """p_perm through the module's own permutation sampler."""
    z_obs = P10._zprop(X, is_b)
    Z, _ = P10._perm_z_sample(X.copy(), is_b, B,
                              np.random.default_rng(seed))
    with np.errstate(invalid="ignore"):
        exceed = (np.abs(Z) >= np.abs(z_obs)[None, :]).sum(axis=0)
    return (exceed + 1) / (B + 1), z_obs


def test_p10_perm_p_counts_ties_at_observed_state():
    """Regression for the float32 Z_perm storage bug: ties at the observed
    |Z| must count as exceedances.

    4x2 margins (rows [1,1,1,1], cols [2,2]): |Z| takes {0, 2} with
    P(|Z| = 2) = 1/3, and the extreme states tie with the observed one —
    those ties must count as exceedances (float32 storage rounded them
    just below the float64 |Z_obs| for non-representable z values, and p
    collapsed to 1/(B+1)).
    """
    X = np.array([[1, 0], [1, 0], [0, 1], [0, 1]], dtype=np.int8)
    is_b = np.array([False, False, True, True])
    p, z = _perm_p(X, is_b, seed=3, B=199)
    assert abs(z[0]) == pytest.approx(2.0)
    assert 0.15 <= p[0] <= 0.6          # exact = 1/3; lost ties -> 0.005
    Z, _ = P10._perm_z_sample(X.copy(), is_b, 5, np.random.default_rng(1))
    assert Z.dtype == np.float64        # float32 storage loses ties


def test_perturb_mixes_at_reasonable_rate():
    """Guard against the frozen-chain regression: with a batch far larger
    than the matrix, the collision filter dropped nearly every valid swap
    (2 applied out of 4800 attempts) and the permutation p-values
    collapsed to 1/(B+1) or 1."""
    rng = np.random.default_rng(2)
    X = (rng.random((60, 8)) < 0.4).astype(np.int8)
    attempts = 4 * X.size
    done = P10._perturb(X, rng, n_attempt=attempts)
    assert done >= 0.03 * attempts        # observed ~9%; frozen chain ~0.1%


def test_p10_perm_p_valid_under_null():
    """Null X (same alteration rate in both groups): p_perm must be
    *super-uniform* — the exact conditional test is valid but conservative,
    because the margin-fixed Z is discrete (the modal allocation gives
    Z = 0 exactly, and ``>=`` counts it as an exceedance).  What must never
    happen is anti-conservative behaviour: more small p than uniform.
    """
    rng = np.random.default_rng(11)
    n_a = n_b = 30
    ps = []
    for rep in range(30):
        X = (rng.random((n_a + n_b, 8)) < 0.4).astype(np.int8)
        is_b = np.array([False] * n_a + [True] * n_b)
        p, _ = _perm_p(X, is_b, seed=1000 + rep)
        ps.extend(p[np.isfinite(p)].tolist())
    ps = np.asarray(ps)
    assert len(ps) > 150
    # validity: P(p <= alpha) never exceeds the uniform rate by much
    for alpha, slack in ((0.01, 0.03), (0.05, 0.08), (0.10, 0.12)):
        assert (ps <= alpha).mean() <= alpha + slack
    # not systematically deflated either (would hide real signal)
    assert ps.mean() >= 0.4


def test_p10_perm_p_detects_preferential_hitting():
    """Column sums held fixed, a preferentially hit column must light up.

    Every row carries exactly two altered peaks so the global-burden
    margins cannot explain the column-allocation asymmetry — the exact
    quantity P10 conditions on.
    """
    n_a = n_b = 20
    rng = np.random.default_rng(21)
    X = np.zeros((n_a + n_b, 6), dtype=np.int8)
    X[:2, 0] = 1            # a:  2/20 hit col 0
    X[20:35, 0] = 1         # b: 15/20 hit col 0 — column sum 17
    for i in range(n_a + n_b):           # fill to exactly 2 per row
        n_fill = 1 if X[i, 0] else 2
        X[i, rng.choice(np.arange(1, 6), size=n_fill, replace=False)] = 1
    assert (X.sum(axis=1) == 2).all()    # equal row margins
    is_b = np.array([False] * n_a + [True] * n_b)
    p, z = _perm_p(X, is_b, seed=5, B=199)
    assert z[0] > 3.0
    assert p[0] <= 0.01                  # strong signal survives margin fix
    assert np.median(p[1:]) > 0.1        # untouched columns stay null-ish
