"""Figures for module PD (baseline/batch diagnostic)."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ._core import *  # noqa: F401,F403

_LAYER_COLORS = {
    "robust": "#009E73",          # okabe bluish green
    "baseline_only": "#E69F00",   # orange
    "centered_only": "#CC79A7",   # reddish purple
    "ns": "#D9D9D9",
}


def plot_pd_baseline_dists(seg_values: dict, baselines: pd.DataFrame, cfg
                           ) -> None:
    """Pooled seg_mean distributions + per-sample baseline, per cohort."""
    if baselines is None or baselines.empty or not seg_values:
        return
    style(cfg)
    pal = _palette(cfg)
    cids = list(dict.fromkeys(baselines["cohort"]))
    fig, axes = plt.subplots(1, 2, figsize=theme.figsize(
        cfg, panel_h=2.6, total_w=7.0))
    ax = axes[0]
    bins = np.linspace(-1.5, 1.5, 90)
    for cid in cids:
        v = np.asarray(seg_values.get(cid, []), dtype=float)
        if len(v) == 0:
            continue
        ax.hist(v, bins=bins, density=True, histtype="step",
                lw=1.1, color=pal.get(cid, "#888888"),
                label=theme.cohort_label(cfg, cid))
    ax.axvline(0, color="#999999", lw=0.6, ls=":")
    ax.set_xlabel("Segment_Mean (log2)", fontsize=7)
    ax.set_ylabel("density", fontsize=7)
    ax.set_title("pooled segment distribution", fontsize=7)
    ax.legend(fontsize=5.5)
    ax.tick_params(labelsize=6.5)
    ax.grid(alpha=0.2)

    ax = axes[1]
    order = cids
    data_list = [baselines.loc[baselines["cohort"] == c, "baseline"]
                 .dropna().to_numpy(float) for c in order]
    bp = ax.boxplot(data_list, tick_labels=order, patch_artist=True,
                    widths=0.5, medianprops=dict(color="#333333"))
    for patch, c in zip(bp["boxes"], order):
        patch.set_facecolor(pal.get(c, "#888888"))
        patch.set_alpha(0.45)
    rng = np.random.default_rng(1)
    for xi, (c, vals) in enumerate(zip(order, data_list), start=1):
        ax.scatter(rng.normal(xi, 0.07, len(vals)), vals, s=6,
                   color=pal.get(c, "#888888"), alpha=0.6, edgecolor="none",
                   zorder=4)
    ax.axhline(0, color="#999999", lw=0.6, ls=":")
    ax.set_ylabel("per-sample baseline\n(median of neutral segments)",
                  fontsize=7)
    ax.set_title("baseline per sample", fontsize=7)
    ax.tick_params(labelsize=6.5)
    ax.grid(axis="y", alpha=0.2)
    fig.suptitle("Module PD \u2014 baseline check", fontsize=8)
    theme.finalize(fig, "tight")
    save(fig, cfg, "PD_baseline_dists")


def plot_pd_shift_track(track: pd.DataFrame, cfg) -> None:
    """Genome-wide median difference (CUE - CHB), raw vs centered."""
    if track is None or track.empty:
        return
    style(cfg)
    comps = list(dict.fromkeys(track["comparison"]))
    nrows, _ = theme.grid_shape(len(comps), max_cols=1)
    fig, axes = plt.subplots(nrows, 1, figsize=theme.figsize(
        cfg, panel_h=2.4, total_w=7.0), squeeze=False)
    for row, comp in enumerate(comps):
        ax = axes[row][0]
        d = track[track["comparison"] == comp].copy()
        d["x"] = _genome_x(d)
        la, lb, ca, cb = _sides(cfg, comp)
        ax.scatter(d["x"], d["delta_median_raw"], s=3, color="#D0D0D0",
                   edgecolor="none", zorder=2)
        ax.plot(d["x"], d["raw_smooth"], color=theme.OKABE["vermillion"],
                lw=1.1, label=f"raw (\u0394median, {la} \u2212 {lb})",
                zorder=4)
        ax.plot(d["x"], d["cen_smooth"], color=theme.OKABE["blue"], lw=1.1,
                label="after baseline removal", zorder=4)
        ax.axhline(0, color="#555555", lw=0.6)
        bnds = _chrom_bounds(d)
        ax.set_xticks([b for b, _ in bnds])
        ax.set_xticklabels([c.replace("chr", "") for _, c in bnds],
                           fontsize=5)
        ax.set_ylabel("\u0394median seg_mean", fontsize=7)
        ax.set_title(f"{comp}: is the shift uniform (batch) or regional?",
                     fontsize=7)
        ax.legend(fontsize=5.5, loc="upper right")
        ax.tick_params(labelsize=6)
        ax.grid(alpha=0.2)
    theme.finalize(fig, "tight")
    save(fig, cfg, "PD_shift_track")


def plot_pd_lambda_bars(lambda_df: pd.DataFrame, cfg) -> None:
    """Lambda raw vs centered per p-value set (healthy lambda ~ 1)."""
    if lambda_df is None or lambda_df.empty:
        return
    style(cfg)
    comps = list(dict.fromkeys(lambda_df["comparison"]))
    nrows, _ = theme.grid_shape(len(comps), max_cols=1)
    fig, axes = plt.subplots(nrows, 1, figsize=theme.figsize(
        cfg, panel_h=2.3, total_w=5.0), squeeze=False)
    for row, comp in enumerate(comps):
        ax = axes[row][0]
        d = lambda_df[lambda_df["comparison"] == comp]
        sets = list(dict.fromkeys(d["set_col"]))
        x = np.arange(len(sets))
        w = 0.35
        for off, (variant, color) in enumerate(
                (("raw", "#888888"), ("centered", "#0072B2"))):
            sub = d[d["variant"] == variant].set_index("set_col")[
                "lambda_infl"].reindex(sets)
            ax.bar(x + (off - 0.5) * w, sub, width=w, color=color,
                   edgecolor="white", linewidth=0.3,
                   label=f"{variant}")
        ax.axhline(1.0, color="#B22222", lw=0.8, ls="--")
        ax.text(len(sets) - 0.5, 1.0, " \u03bb = 1", fontsize=5.5,
                color="#B22222", va="bottom", ha="right")
        ax.set_xticks(x)
        ax.set_xticklabels(sets, fontsize=6.5)
        ax.set_ylabel("\u03bb", fontsize=7)
        ax.set_title(f"{comp}", fontsize=7)
        ax.legend(fontsize=5.5)
        ax.tick_params(labelsize=6.5)
        ax.grid(axis="y", alpha=0.2)
    theme.finalize(fig, "tight")
    save(fig, cfg, "PD_lambda_bars")


def figures_pd(res: dict, cfg) -> None:
    plot_pd_baseline_dists(res.get("seg_values"), res.get("baselines"), cfg)
    plot_pd_shift_track(res.get("track"), cfg)
    plot_pd_lambda_bars(res.get("lambda_df"), cfg)
