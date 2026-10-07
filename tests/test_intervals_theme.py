"""Interval utilities and stats helpers."""

from __future__ import annotations

import numpy as np
import pandas as pd

from gistic_compare import intervals as iv
from gistic_compare.theme import bh_adjust, wilson_ci


def test_overlap_len():
    ov = iv.overlap_len(100, 200, [150, 250], [300, 300])
    assert ov[0] == 50
    assert ov[1] == 0


def _peaks(cohort, rows):
    return pd.DataFrame(
        rows, columns=["chrom", "start", "end", "type", "q"],
    ).assign(cohort=cohort,
             peak_id=[f"{cohort}:{i}" for i in range(len(rows))])


def test_match_peaks_and_jaccard():
    a = _peaks("a", [("chr1", 1000, 2000, "Amp", 0.01),
                     ("chr2", 1000, 2000, "Del", 0.3)])
    b = _peaks("b", [("chr1", 1100, 2100, "Amp", 0.02),
                     ("chr3", 500, 900, "Del", 0.01)])
    res = iv.match_peaks(a, b, min_bp=50, reciprocal=0.5, mode="either")
    assert isinstance(res, iv.MatchResult)
    assert len(res.pairs) == 1          # only chr1 amp matches
    assert len(res.only_a) == 1 and len(res.only_b) == 1
    j = iv.jaccard(len(a), len(b), len(res.pairs))
    assert 0.0 < j < 1.0


def test_match_peaks_no_cross_chromosome():
    a = _peaks("a", [("chr1", 1000, 2000, "Amp", 0.01)])
    b = _peaks("b", [("chr2", 1000, 2000, "Amp", 0.01)])
    res = iv.match_peaks(a, b, min_bp=50, reciprocal=0.5, mode="either")
    assert res.pairs.empty
    # matching is type-agnostic at this layer (callers split by type first)
    a2 = _peaks("a", [("chr1", 1000, 2000, "Amp", 0.01)])
    b2 = _peaks("b", [("chr1", 1000, 2000, "Del", 0.01)])
    res2 = iv.match_peaks(a2, b2, min_bp=50, reciprocal=0.5, mode="either")
    assert len(res2.pairs) == 1


def test_bh_adjust_monotone():
    p = np.array([0.001, 0.01, 0.02, 0.5, 0.9])
    q = bh_adjust(p)
    assert np.all(np.diff(q) >= -1e-12)
    assert q[0] <= 0.005 + 1e-9
    assert np.all(q >= p - 1e-12)


def test_bh_adjust_nan_safe():
    p = np.array([np.nan, 0.01, 0.02])
    q = bh_adjust(p)
    assert np.isnan(q[0])
    assert np.isfinite(q[1:]).all()


def test_wilson_ci_contains_estimate():
    lo, hi = wilson_ci(7, 10)
    assert 0.0 < lo < 0.7 < hi < 1.0
    lo0, hi0 = wilson_ci(0, 10)
    assert lo0 == 0.0 and 0 < hi0 < 0.4
