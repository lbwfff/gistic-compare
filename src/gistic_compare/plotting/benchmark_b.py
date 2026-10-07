"""Figures for module B (score tracks).

(Split from the original single-module ``plotting``.)
"""

from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ._core import *  # noqa: F401,F403

# ---------------------------------------------------------------------------
# B: score tracks
# ---------------------------------------------------------------------------

def plot_track_correlations(corr: pd.DataFrame, cfg: dict[str, Any],
                            method: str = "spearman", *,
                            prefix: str = "B") -> None:
    """Heatmap of Pearson or Spearman r: rows = comparison x ptype, cols = metric."""
    if corr.empty:
        return
    if method not in ("pearson", "spearman"):
        raise ValueError(f"unknown method: {method!r}")
    value_col = f"{method}_r"
    if value_col not in corr.columns:
        return
    style(cfg)
    metric_order = ["neg_log10_q", "g_score", "frequency"]
    metric_labels = ["-log10(q)", "G-score", "frequency"]
    sub = corr[corr["metric"].isin(metric_order)].copy()
    if sub.empty:
        return
    sub["row"] = sub["comparison"] + " | " + sub["ptype"]
    mat = sub.pivot_table(index="row", columns="metric",
                          values=value_col, aggfunc="first")
    mat = mat.reindex(columns=metric_order)
    row_order = []
    for comp in sub["comparison"].unique():
        for pt in ("Amp", "Del"):
            r = f"{comp} | {pt}"
            if r in mat.index:
                row_order.append(r)
    mat = mat.reindex(index=[r for r in row_order if r in mat.index] +
                            [r for r in mat.index if r not in row_order])

    cbar_label = "Pearson r" if method == "pearson" else "Spearman r"
    vals = mat.to_numpy(dtype=float)
    # a correlation matrix must be able to show negative r -> diverging scale
    has_neg = bool(np.isfinite(vals).any() and np.nanmin(vals) < 0)
    cmap = "RdBu_r" if has_neg else "Blues"
    vmin, vmax = (-1.0, 1.0) if has_neg else (0.0, 1.0)
    n_bins = np.nan
    if "n_bins" in sub.columns and len(sub):
        n_bins = float(np.nanmax(sub["n_bins"].to_numpy(dtype=float)))
    fig, ax = plt.subplots(
        figsize=(theme.COL_WIDTH + 1.1, 0.34 * len(mat) + 1.2))
    sns.heatmap(mat, annot=True, fmt=".2f", cmap=cmap,
                vmin=vmin, vmax=vmax, linewidths=0.6, linecolor="white",
                annot_kws={"fontsize": 6.5},
                cbar_kws={"label": cbar_label, "shrink": 0.8}, ax=ax)
    ax.set_xticklabels(metric_labels, rotation=30, ha="right")
    ax.set_xlabel(f"n = {int(n_bins)} bins" if np.isfinite(n_bins) else "")
    ax.set_ylabel("")
    theme.finalize(fig, "tight")
    save(fig, cfg, f"{prefix}_track_correlations_{method}")


def plot_genome_tracks(grid: pd.DataFrame, cfg: dict[str, Any],
                       ptype: str = "Amp") -> None:
    """Stacked genome-wide -log10(q) tracks for one ptype (Amp or Del)."""
    if grid is None or grid.empty:
        return
    style(cfg)
    reset = (grid.reset_index()
             if isinstance(grid.index, pd.MultiIndex)
             or grid.index.name in ("chrom", "start", "end")
             else grid.copy())
    cohorts = [c for c in cfg["cohorts"]
               if f"{c}_{ptype}" in reset.columns]
    if not cohorts:
        return
    n = len(cohorts)
    fig, axes = plt.subplots(
        n, 1, figsize=theme.figsize(cfg, nrows=n, panel_h=1.75),
        sharex=True, squeeze=False)
    reset = reset.reset_index(drop=True)
    x = np.arange(len(reset))
    chrom_bounds = []
    prev = None
    for i, c in enumerate(reset["chrom"]):
        if c != prev:
            chrom_bounds.append((i, c))
            prev = c
    pal = _palette(cfg)
    color = _pc(ptype)
    for i, cid in enumerate(cohorts):
        ax = axes[i][0]
        col = f"{cid}_{ptype}"
        vals = reset[col].to_numpy(dtype=float)
        ax.fill_between(x, 0, np.nan_to_num(vals), alpha=0.85,
                        color=pal.get(cid, color), linewidth=0)
        ax.set_ylabel(cfg["cohorts"][cid].get("label", cid), fontsize=7)
        ax.set_xlim(0, len(reset) - 1)
        ax.set_ylim(bottom=0)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
        ax.tick_params(axis="y", pad=2)
        # alternating chromosome shading + boundaries between chromosomes
        for k, (pos, _) in enumerate(chrom_bounds):
            if k % 2 == 1:
                end = chrom_bounds[k + 1][0] if k + 1 < len(chrom_bounds) else len(reset)
                ax.axvspan(pos, end, color="#F2F2F2", zorder=0)
            if k > 0:
                ax.axvline(pos, color="#CCCCCC", lw=0.5, zorder=1)
    axes[-1][0].set_xticks([p for p, _ in chrom_bounds])
    axes[-1][0].set_xticklabels([c.replace("chr", "") for _, c in chrom_bounds],
                                fontsize=6.5)
    axes[-1][0].set_xlabel("chromosome")
    axes[-1][0].tick_params(axis="x", pad=2)
    fig.supylabel("-log10 q", fontsize=8)
    theme.finalize(fig, "constrained", hspace=0.35)
    save(fig, cfg, f"B_genome_tracks_{ptype.lower()}")


