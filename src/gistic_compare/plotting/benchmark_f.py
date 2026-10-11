"""Figures for module F (novelty vs control).

(Split from the original single-module ``plotting``.)
"""

from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ._core import *  # noqa: F401,F403

# ---------------------------------------------------------------------------
# F: novelty (experiment-specific signal vs control)
# ---------------------------------------------------------------------------

_NOV_WEIGHTS = [
    ("frac_count", "count"),
    ("frac_freq", "frequency"),
    ("frac_neglogq", "-log10(q)"),
    ("frac_gscore", "G-score"),
]


def plot_novelty_facet(summary: pd.DataFrame, cfg: dict[str, Any]) -> None:
    """One panel per weight; x=experiment, hue=Amp/Del, y=fraction.

    Whiskers are drawn **on the count panel only**: the 95% Wilson interval
    of ``n_unique / n_exp`` matches the bar quantity there and only there.
    The weighted fractions (frequency / -log10(q) / G-score) are ratios of
    summed weights with no closed-form binomial CI, so they are drawn as
    bare bars rather than wearing a CI that belongs to another quantity.
    """
    if summary.empty:
        return
    style(cfg)
    q_main = cfg["analysis"]["q_thresholds"][0]
    sub = summary[summary["q_thresh"] == q_main].copy()
    if sub.empty:
        return
    weights = [(c, lab) for c, lab in _NOV_WEIGHTS if c in sub.columns]
    if not weights:
        return
    exp_order = [e for e in cfg["cohorts"] if e in set(sub["experiment"])]
    ptypes = [p for p in ("Amp", "Del") if p in set(sub["ptype"])]
    nrows, ncols = theme.grid_shape(len(weights), max_cols=2)
    fig, axes = plt.subplots(nrows, ncols, squeeze=False, sharex=True,
                             sharey=True,
                             figsize=(theme.width(cfg), 2.6 * nrows))
    flat = axes.ravel()
    w = 0.34
    for k, (col, lab) in enumerate(weights):
        ax = flat[k]
        has_ci = col == "frac_count"
        for pi, pt in enumerate(ptypes):
            xs, ys, lo, hi = [], [], [], []
            for ei, e in enumerate(exp_order):
                r = sub[(sub["experiment"] == e) & (sub["ptype"] == pt)]
                if r.empty:
                    continue
                r = r.iloc[0]
                if has_ci:
                    cl, ch = theme.wilson_ci([r["n_unique"]], [r["n_exp"]])
                    clo, chi = float(cl[0]), float(ch[0])
                else:
                    clo, chi = np.nan, np.nan
                xs.append(ei + (pi - 0.5) * (w + 0.04))
                ys.append(float(r[col]))
                lo.append(clo)
                hi.append(chi)
            if not xs:
                continue
            bars = theme.bar_with_ci(ax, xs, ys, lo, hi, width=w,
                                     color=_pc(pt), label=pt,
                                     edgecolor="white", linewidth=0.4)
            for xv, yv, hv in zip(xs, ys, hi):
                top = yv if not np.isfinite(hv) else max(yv, hv)
                ax.text(xv, top, f"{yv:.2f}", ha="center", va="bottom",
                        fontsize=5, color="#444444")
        ax.set_title(lab, fontsize=8)
        ax.set_ylim(0, 1.20)
        ax.set_xticks(range(len(exp_order)))
        ax.set_xticklabels([theme.cohort_label(cfg, e) for e in exp_order],
                           rotation=30, ha="right", fontsize=6.5)
        theme.panel_label(ax, theme.panel_letter(k), cfg)
    for ax in flat[len(weights):]:
        ax.set_visible(False)
    flat[0].set_ylabel("fraction unique to experiment")
    if ptypes:
        flat[0].legend(frameon=False, fontsize=6, loc="upper left", ncol=2)
    # figure-level footnote: inside the last panel it lands on the bars
    fig.text(0.5, 0.004,
             "whiskers (count panel): 95% Wilson CI of n_unique/n_exp; "
             "weighted fractions have no closed-form CI",
             fontsize=5.5, color="#555555", ha="center", va="bottom")
    fig.tight_layout(rect=(0.0, 0.028, 1.0, 1.0))
    save(fig, cfg, "F_novelty_facet")


def plot_novelty_bars(summary: pd.DataFrame, cfg: dict[str, Any]) -> None:
    """Single grouped bar: x=weight, hue=experiment, facet=ptype.

    Whiskers appear on the ``count`` weight only (95% Wilson CI of
    n_unique/n_exp); the weighted fractions have no closed-form CI."""
    if summary.empty:
        return
    style(cfg)
    q_main = cfg["analysis"]["q_thresholds"][0]
    sub = summary[summary["q_thresh"] == q_main].copy()
    if sub.empty:
        return
    weights = [(c, lab) for c, lab in _NOV_WEIGHTS if c in sub.columns]
    if not weights:
        return
    exp_order = [e for e in cfg["cohorts"] if e in set(sub["experiment"])]
    ptypes = [p for p in ("Amp", "Del") if p in set(sub["ptype"])]
    pal = _palette(cfg)
    nrows, ncols = theme.grid_shape(len(ptypes), max_cols=2)
    fig, axes = plt.subplots(nrows, ncols, squeeze=False, sharey=True,
                             figsize=(theme.width(cfg), 3.4 * nrows))
    flat = axes.ravel()
    x = np.arange(len(weights))
    bw = 0.8 / max(len(exp_order), 1)
    for k, pt in enumerate(ptypes):
        ax = flat[k]
        for ei, e in enumerate(exp_order):
            xs, ys, lo, hi = [], [], [], []
            for wi, (col, _lab) in enumerate(weights):
                r = sub[(sub["experiment"] == e) & (sub["ptype"] == pt)]
                if r.empty:
                    continue
                r = r.iloc[0]
                # CI only where it matches the bar quantity (count fraction)
                if col == "frac_count":
                    cl, ch = theme.wilson_ci([r["n_unique"]], [r["n_exp"]])
                    clo, chi = float(cl[0]), float(ch[0])
                else:
                    clo, chi = np.nan, np.nan
                xs.append(wi + (ei - (len(exp_order) - 1) / 2) * bw)
                ys.append(float(r[col]))
                lo.append(clo)
                hi.append(chi)
            if not xs:
                continue
            bars = theme.bar_with_ci(
                ax, xs, ys, lo, hi, width=bw * 0.9,
                color=pal.get(e, "#333333"), label=theme.cohort_label(cfg, e),
                edgecolor="white", linewidth=0.4)
            for xv, yv, hv in zip(xs, ys, hi):
                top = yv if not np.isfinite(hv) else max(yv, hv)
                ax.text(xv, top, f"{yv:.2f}", ha="center", va="bottom",
                        fontsize=4.5, color="#444444")
        ax.set_title(pt, fontsize=8)
        ax.set_ylim(0, 1.20)
        ax.set_xticks(x)
        ax.set_xticklabels([lab for _c, lab in weights], rotation=20,
                           ha="right", fontsize=6.5)
        theme.panel_label(ax, theme.panel_letter(k), cfg)
    for ax in flat[len(ptypes):]:
        ax.set_visible(False)
    flat[0].set_ylabel("fraction unique to experiment")
    flat[0].set_xlabel("weight")
    h, lab = flat[0].get_legend_handles_labels()
    if h:
        fig.legend(h, lab, loc="lower center", ncol=len(exp_order),
                   bbox_to_anchor=(0.5, 1.0), frameon=False, fontsize=6)
    theme.finalize(fig, "tight")
    save(fig, cfg, "F_novelty_bars")


