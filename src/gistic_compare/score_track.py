"""Module B: genome-wide score track resampling and correlation."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from .config import CHROMS, HG19_CHROM_SIZES, out_table


def make_grid(bin_size: int) -> pd.DataFrame:
    rows = []
    for c in CHROMS:
        L = HG19_CHROM_SIZES[c]
        starts = np.arange(0, L, bin_size, dtype=np.int64)
        ends = np.minimum(starts + bin_size, L)
        for s, e in zip(starts, ends):
            rows.append((c, int(s), int(e)))
    return pd.DataFrame(rows, columns=["chrom", "start", "end"])


def resample_track(scores: pd.DataFrame, grid: pd.DataFrame,
                   value_col: str, ptype: str) -> pd.Series:
    """Length-weighted mean of `value_col` over grid bins for one Type.

    Returns Series indexed like grid (same length), NaN where no coverage.
    """
    sub = scores[scores["type"] == ptype]
    vals = np.full(len(grid), np.nan)
    # pre-split grid by chrom
    grid_idx: dict[str, np.ndarray] = {
        c: np.flatnonzero((grid["chrom"] == c).to_numpy()) for c in CHROMS
    }
    weight_sum = np.zeros(len(grid))
    value_sum = np.zeros(len(grid))

    for c in CHROMS:
        gidx = grid_idx[c]
        if len(gidx) == 0:
            continue
        g_start = grid["start"].to_numpy()[gidx]
        g_end = grid["end"].to_numpy()[gidx]
        srows = sub[sub["chrom"] == c]
        if srows.empty:
            continue
        s_start = srows["start"].to_numpy()
        s_end = srows["end"].to_numpy()
        s_val = srows[value_col].to_numpy(dtype=float)
        # for each segment, find overlapping bins
        for i in range(len(srows)):
            ov_start = np.maximum(g_start, s_start[i])
            ov_end = np.minimum(g_end, s_end[i])
            ov = ov_end - ov_start
            m = ov > 0
            if not m.any():
                continue
            w = ov[m].astype(float)
            idx = gidx[m]
            # ignore NaN values
            v = s_val[i]
            if np.isnan(v):
                continue
            value_sum[idx] += w * v
            weight_sum[idx] += w
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(weight_sum > 0, value_sum / weight_sum, np.nan)
    return pd.Series(out, index=grid.index)


def run_module_b(cfg: dict[str, Any], data: dict[str, dict]) -> dict[str, Any]:
    ana = cfg["analysis"]
    bin_size = ana["score_track"]["bin_size"]
    grid = make_grid(bin_size)
    metrics = ["neg_log10_q", "g_score", "frequency"]

    # build per-cohort × type × metric vectors aligned to grid
    tracks: dict[tuple[str, str, str], np.ndarray] = {}
    for cid, d in data.items():
        sc = d["scores"]
        for ptype in ("Amp", "Del"):
            for m in metrics:
                tracks[(cid, ptype, m)] = resample_track(sc, grid, m, ptype).to_numpy()

    # pairwise correlations
    corr_rows = []
    for cmp_ in cfg["comparisons"]:
        ca, cb = cmp_["a"], cmp_["b"]
        for ptype in ("Amp", "Del"):
            for m in metrics:
                xa = tracks[(ca, ptype, m)]
                xb = tracks[(cb, ptype, m)]
                mask = np.isfinite(xa) & np.isfinite(xb)
                n = int(mask.sum())
                if n >= 3:
                    r_p, p_p = stats.pearsonr(xa[mask], xb[mask])
                    r_s, p_s = stats.spearmanr(xa[mask], xb[mask])
                else:
                    r_p = p_p = r_s = p_s = np.nan
                corr_rows.append(dict(
                    comparison=f"{ca}_vs_{cb}",
                    kind=cmp_.get("kind", "unpaired"),
                    ptype=ptype, metric=m, n_bins=n,
                    pearson_r=r_p, pearson_p=p_p,
                    spearman_r=r_s, spearman_p=p_s,
                ))

    # top differing bins (using neg_log10_q, Amp+Del combined via max)
    diff_rows = []
    for cmp_ in cfg["comparisons"]:
        ca, cb = cmp_["a"], cmp_["b"]
        for ptype in ("Amp", "Del"):
            xa = tracks[(ca, ptype, "neg_log10_q")]
            xb = tracks[(cb, ptype, "neg_log10_q")]
            delta = np.nan_to_num(xa, nan=0.0) - np.nan_to_num(xb, nan=0.0)
            order = np.argsort(-np.abs(delta))[:50]
            for i in order:
                if delta[i] == 0:
                    continue
                diff_rows.append(dict(
                    comparison=f"{ca}_vs_{cb}", ptype=ptype,
                    chrom=grid["chrom"].iloc[i],
                    start=int(grid["start"].iloc[i]),
                    end=int(grid["end"].iloc[i]),
                    delta_neg_log10_q=float(delta[i]),
                    a_val=float(xa[i]) if np.isfinite(xa[i]) else np.nan,
                    b_val=float(xb[i]) if np.isfinite(xb[i]) else np.nan,
                ))

    # save long grid table for significant bins (q>1.3 ~ q<0.05)
    grid_out = grid.copy()
    for (cid, ptype, m), vec in tracks.items():
        if m != "neg_log10_q":
            continue
        col = f"{cid}_{ptype}"
        grid_out[col] = vec

    corr_df = pd.DataFrame(corr_rows)
    diff_df = pd.DataFrame(diff_rows)
    grid_out.to_csv(out_table(cfg, "B_score_grid_1mb.tsv"), sep="\t", index=False)
    corr_df.to_csv(out_table(cfg, "B_track_correlations.tsv"), sep="\t", index=False)
    diff_df.to_csv(out_table(cfg, "B_top_diff_bins.tsv"), sep="\t", index=False)
    grid_indexed = grid_out.set_index(["chrom", "start", "end"])
    return dict(grid=grid_indexed, tracks=tracks, correlations=corr_df, top_diff=diff_df)
