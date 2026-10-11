"""Regression tests for silent data-misalignment bugs (review findings H1-H3).

Each test constructs a case where the *positional* pairing used to be wrong
and asserts values are attached by identity (peak_id / window id / sample
label) instead.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from gistic_compare import novelty, pairwise_diag
from gistic_compare import pairwise_permutation as P10

# ---------------------------------------------------------------------------
# H1: novelty module F — G-scores attached to peaks by peak_id
# ---------------------------------------------------------------------------

def _fake_scores_peaks():
    """Full peaks table in sort_peaks order (Amp block first, then Del).

    Each peak has a distinctive G-score so a positional mix-up is visible:
    Amp peaks carry 10/20/30, Del peaks carry 0.1/0.2.
    """
    peaks = pd.DataFrame(dict(
        peak_id=["c:Amp:0", "c:Amp:1", "c:Amp:2", "c:Del:3", "c:Del:4"],
        type=["Amp"] * 3 + ["Del"] * 2,
        chrom=["chr1"] * 5,
        start=[0, 100, 200, 0, 100],
        end=[50, 150, 250, 50, 150]))
    rows = []
    for pid, g in zip(peaks["peak_id"], [10.0, 20.0, 30.0, 0.1, 0.2]):
        r = peaks.loc[peaks["peak_id"] == pid].iloc[0]
        rows.append(dict(type=r["type"], chrom=r["chrom"], start=r["start"],
                         end=r["end"], g_score=g))
    return pd.DataFrame(rows), peaks


def test_aggregate_gscore_indexed_by_peak_id():
    scores, peaks = _fake_scores_peaks()
    gs = novelty._aggregate_gscore(scores, peaks)
    assert (gs.index == peaks["peak_id"]).all()
    assert gs["c:Del:3"] == pytest.approx(0.1)
    assert gs["c:Amp:2"] == pytest.approx(30.0)


def test_attach_gscore_filtered_subset_by_id():
    """A q-filtered Del-only subset must keep its own G-scores.

    Positional reindexing hands the Del peaks the first m rows of the full
    table — i.e. the Amp peaks' G-scores (10/20 instead of 0.1/0.2).
    """
    scores, peaks = _fake_scores_peaks()
    gs = novelty._aggregate_gscore(scores, peaks)
    fe = peaks[peaks["type"] == "Del"].reset_index(drop=True)  # as filter_peaks
    got = novelty._attach_gscore(fe, gs)
    np.testing.assert_allclose(got, [0.1, 0.2])   # not [10.0, 20.0]
    # reordered frame must follow the frame, not the source order
    rev = fe.iloc[::-1].reset_index(drop=True)
    np.testing.assert_allclose(novelty._attach_gscore(rev, gs), [0.2, 0.1])


# ---------------------------------------------------------------------------
# H2: PD — raw vs centered window tables paired by window id
# ---------------------------------------------------------------------------

def _win_stats(wins, deltas, p_mw):
    return pd.DataFrame(dict(
        win=wins,
        chrom=["chr1"] * len(wins),
        start=[w * 1_000_000 for w in wins],
        end=[(w + 1) * 1_000_000 for w in wins],
        delta_median=[float(d) for d in deltas],
        p_mw=[float(p) for p in p_mw],
        p_mw_fdr=[float(p) for p in p_mw],
        p_amp_fdr=[1.0] * len(wins)))


def test_trust_layer_pairs_by_window_id():
    """Centered variant skips win=1 (e.g. NaN-baseline sample); rows after
    the gap must not shift onto the wrong raw window."""
    raw = _win_stats([0, 1, 2, 3], [10.0, 20.0, 30.0, 40.0],
                     [0.001, 0.001, 0.5, 0.5])
    cen = _win_stats([0, 2, 3], [100.0, 300.0, 400.0],
                     [0.001, 0.001, 0.5])
    raw_map, cen_map = {0: 10.0, 1: 20.0, 2: 30.0, 3: 40.0}, \
        {0: 100.0, 2: 300.0, 3: 400.0}
    layer = pairwise_diag._trust_layer(raw, cen, fdr=0.05)
    # inner join on win: window 1 (missing from centered) is dropped
    assert sorted(layer["win"]) == [0, 2, 3]
    for _, r in layer.iterrows():
        w = int(r["win"])
        assert r["delta_median_raw"] == raw_map[w]
        assert r["delta_median_cen"] == cen_map[w]   # never the shifted row
    by_win = layer.set_index("win")["layer"].to_dict()
    assert by_win[0] == "robust"          # sig in both
    assert by_win[2] == "centered_only"   # sig only after centering
    assert by_win[3] == "ns"


# ---------------------------------------------------------------------------
# H3: P10 — binary matrix values follow sample labels when shared ids drop
# ---------------------------------------------------------------------------

def test_binary_matrix_shared_sample_alignment():
    """With a shared sample id dropped from an unpaired comparison, the
    remaining samples must keep their own CN values (the old positional
    slice shifted every value after the dropped one)."""
    amp_key, del_key = "Amp:chr1:100-200", "Del:chr1:300-400"
    keys = pd.DataFrame([
        dict(cohort="a", peak_key=amp_key, ptype="Amp", peak_id="p1",
             cytoband="1q21", chrom="chr1", start=100, end=200),
        dict(cohort="b", peak_key=amp_key, ptype="Amp", peak_id="q1",
             cytoband="1q21", chrom="chr1", start=100, end=200),
        dict(cohort="a", peak_key=del_key, ptype="Del", peak_id="p2",
             cytoband="1p36", chrom="chr1", start=300, end=400),
        dict(cohort="b", peak_key=del_key, ptype="Del", peak_id="q2",
             cytoband="1p36", chrom="chr1", start=300, end=400),
    ])
    cn_a = pd.DataFrame([  # samples A1, S (shared)
        dict(peak_id="p1", sample="A1", cn=0.5),
        dict(peak_id="p1", sample="S", cn=0.5),
        dict(peak_id="p2", sample="A1", cn=-0.6),
        dict(peak_id="p2", sample="S", cn=-0.6),
    ])
    cn_b = pd.DataFrame([  # samples S (shared), B1
        dict(peak_id="q1", sample="S", cn=0.5),
        dict(peak_id="q1", sample="B1", cn=0.0),
        dict(peak_id="q2", sample="S", cn=-0.6),
        dict(peak_id="q2", sample="B1", cn=0.0),
    ])
    data = {"a": {"peak_cn": cn_a}, "b": {"peak_cn": cn_b}}
    built = P10._binary_matrix(None, data, "a", "b", keys, cn_thr=0.25)
    assert built is not None
    X, samples, meta = built
    assert samples == ["A1", "B1"]          # shared "S" dropped
    assert list(meta["peak_key"]) == [amp_key, del_key]
    # A1 altered on both peaks; B1 unaltered on both — and must stay so
    assert X[:, 0].tolist() == [1, 0]
    assert X[:, 1].tolist() == [1, 0]
