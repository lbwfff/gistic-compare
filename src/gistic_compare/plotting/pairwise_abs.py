"""Figures for module ABS (absolute copy-number comparisons)."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ._core import *  # noqa: F401,F403


def plot_abs_wgd(wgd: pd.DataFrame, cfg) -> None:
    """WGD+ fraction per cohort with Wilson CIs and the Fisher p."""
    if wgd is None or wgd.empty:
        return
    style(cfg)
    main = wgd[~wgd["comparison"].str.contains("high-purity", na=False)]
    comps = list(dict.fromkeys(main["comparison"]))
    nrows, _ = theme.grid_shape(len(comps), max_cols=2)
    fig, axes = plt.subplots(nrows, 2, figsize=theme.figsize(
        cfg, panel_h=2.6, total_w=6.2), squeeze=False)
    axes = axes.ravel()
    for k, comp in enumerate(comps):
        ax = axes[k]
        r = main[main["comparison"] == comp].iloc[0]
        la, lb, _, _ = _sides(cfg, comp)
        heights = [r["wgd_frac_a"], r["wgd_frac_b"]]
        lo, hi = theme.wilson_ci(
            [r["wgd_n_a"], r["wgd_n_b"]], [r["n_a"], r["n_b"]])
        theme.bar_with_ci(ax, [0, 1], heights, lo, hi,
                          color=[theme.palette(cfg).get(c, "#888888")
                                 for c in (comp.split("_vs_")[0],
                                           comp.split("_vs_")[1])])
        ax.set_xticks([0, 1])
        ax.set_xticklabels([theme.cohort_label(cfg, comp.split("_vs_")[0]),
                            theme.cohort_label(cfg, comp.split("_vs_")[1])],
                           fontsize=7)
        theme.annot_stat(ax, f"Fisher p = {r['p_fisher']:.3g}, "
                             f"OR = {r['odds_ratio']}", loc="upper left")
        ax.set_ylabel("WGD+ fraction", fontsize=7)
        ax.set_ylim(0, max(0.1, float(np.nanmax(hi)) * 1.25))
        ax.set_title(comp, fontsize=7)
        ax.tick_params(labelsize=6.5)
        ax.grid(axis="y", alpha=0.2)
    for k in range(len(comps), len(axes)):
        axes[k].set_visible(False)
    theme.finalize(fig, "tight")
    save(fig, cfg, "ABS_wgd_freq")


def plot_abs_metric_dists(metrics: pd.DataFrame, scores: pd.DataFrame, cfg
                          ) -> None:
    """purity / ploidy / LOH / CN-LOH distributions per cohort."""
    if metrics is None or metrics.empty or scores is None or scores.empty:
        return
    style(cfg)
    pal = _palette(cfg)
    comps = list(dict.fromkeys(metrics["comparison"]))
    nrows, _ = theme.grid_shape(len(comps), max_cols=1)
    fig, axes = plt.subplots(nrows, 4, figsize=theme.figsize(
        cfg, panel_h=2.4, total_w=7.0), squeeze=False)
    panels = [("purity", "purity"), ("ploidy_abs", "absolute ploidy"),
              ("loh_frac", "LOH fraction"),
              ("cnloh_frac", "copy-neutral LOH")]
    for row, comp in enumerate(comps):
        d = scores[scores["cohort"].isin(comp.split("_vs_"))]
        t = metrics[metrics["comparison"] == comp]
        order = list(dict.fromkeys(d["cohort"]))
        for col, (metric, label) in enumerate(panels):
            ax = axes[row][col]
            data = [d.loc[d["cohort"] == c, metric].dropna().to_numpy(float)
                    for c in order]
            if all(len(x) == 0 for x in data):
                ax.set_visible(False)
                continue
            bp = ax.boxplot(data, tick_labels=order, patch_artist=True,
                            widths=0.5, medianprops=dict(color="#333333"))
            for patch, c in zip(bp["boxes"], order):
                patch.set_facecolor(pal.get(c, "#888888"))
                patch.set_alpha(0.45)
            rng = np.random.default_rng(row * 4 + col)
            for xi, vals in enumerate(data, start=1):
                if len(vals) == 0:
                    continue
                ax.scatter(rng.normal(xi, 0.07, len(vals)), vals, s=5,
                           color=pal.get(order[xi - 1], "#888888"),
                           alpha=0.55, edgecolor="none", zorder=4)
            if len(t) and "p_mw" in t.columns:
                r = t[t["metric"] == metric]
                if len(r):
                    pv = r["p_mw"].iloc[0]
                    theme.annot_stat(
                        ax, theme.fmt_p(float(pv)
                                        if np.isfinite(pv) else None),
                        loc="upper left", fontsize=6.0)
            ax.set_title(label, fontsize=6.8)
            ax.tick_params(labelsize=6)
            ax.grid(axis="y", alpha=0.2)
    fig.suptitle("Module ABS \u2014 absolute-CN metrics per sample",
                 fontsize=8)
    theme.finalize(fig, "tight")
    save(fig, cfg, "ABS_metric_dists")


def plot_abs_purity_adjusted(peaks: pd.DataFrame, cfg) -> None:
    """Crude vs purity-adjusted OR per peak: does the cohort effect survive
    purity adjustment?"""
    if peaks is None or peaks.empty:
        return
    d = peaks[np.isfinite(peaks["crude_or"])
              & np.isfinite(peaks["adj_or_purity"])]
    if d.empty:
        return
    style(cfg)
    comp = d["comparison"].iloc[0]
    fig, ax = plt.subplots(figsize=theme.figsize(cfg, panel_h=2.9,
                                                 total_w=3.5))
    sig = d["adj_q_purity"] < 0.05 if "adj_q_purity" in d.columns \
        else pd.Series(False, index=d.index)
    ax.scatter(d.loc[~sig, "crude_or"], d.loc[~sig, "adj_or_purity"],
               s=10, color="#8FA8C8", alpha=0.8, edgecolor="none")
    ax.scatter(d.loc[sig, "crude_or"], d.loc[sig, "adj_or_purity"], s=14,
               color=theme.OKABE["vermillion"], edgecolor="none",
               label="adj q < 0.05")
    lim = max(0.5, float(np.nanmax([d["crude_or"].max(),
                                    d["adj_or_purity"].max()])) * 1.05)
    ax.plot([0, lim], [0, lim], color="#999999", lw=0.7, ls="--")
    ax.axvline(1, color="#777777", lw=0.6, ls=":")
    ax.axhline(1, color="#777777", lw=0.6, ls=":")
    ax.set_xlabel("crude OR (cohort)", fontsize=7)
    ax.set_ylabel("purity-adjusted OR", fontsize=7)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_title(f"{comp}: does the effect survive purity adjustment?",
                 fontsize=7)
    ax.legend(fontsize=5.5, loc="upper left")
    ax.grid(alpha=0.2)
    theme.finalize(fig, "tight")
    save(fig, cfg, "ABS_purity_adjusted")


def figures_abs(res: dict, cfg) -> None:
    plot_abs_wgd(res.get("wgd"), cfg)
    plot_abs_metric_dists(res.get("metrics"), res.get("scores"), cfg)
    plot_abs_purity_adjusted(res.get("purity_peaks"), cfg)
