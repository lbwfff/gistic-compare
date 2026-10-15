"""Figures for module C (sample calls).

(Split from the original single-module ``plotting``.)
"""

from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ._core import *  # noqa: F401,F403

# ---------------------------------------------------------------------------
# C: sample calls
# ---------------------------------------------------------------------------

def plot_agreement(peak_agreement: pd.DataFrame, cfg: dict[str, Any],
                   prefix: str = "C") -> None:
    if peak_agreement.empty:
        return
    style(cfg)
    order = [p for p in ("Amp", "Del")
             if p in set(peak_agreement["ptype"])]
    fig, axes = plt.subplots(1, 2, figsize=theme.figsize(cfg, panel_h=3.0))
    metrics = [("exact_agreement", "exact 0/1/2 agreement"),
               ("weighted_kappa", "quadratic-weighted kappa")]
    for k, (ax, (metric, ylab)) in enumerate(zip(axes, metrics)):
        sns.violinplot(data=peak_agreement, x="ptype", y=metric,
                       order=order, hue="ptype", palette=theme.PTYPE_COLORS,
                       inner="box", cut=0, saturation=1, legend=False, ax=ax)
        theme.strip_overlay(ax, peak_agreement, "ptype", metric, order=order,
                            hue="ptype", palette=theme.PTYPE_COLORS,
                            dodge=True)
        ax.set_ylim(-0.03, 1.06)
        ax.set_xlabel("")
        ax.set_ylabel(ylab)
        theme.annot_n(ax, theme.group_sizes(peak_agreement, "ptype", order),
                      loc="lower left")
        theme.panel_label(ax, theme.panel_letter(k), cfg)
    theme.finalize(fig, "tight")
    save(fig, cfg, f"{prefix}_peak_agreement")


def plot_sample_agreement(sample_agreement: pd.DataFrame, cfg: dict[str, Any],
                          prefix: str = "C") -> None:
    if sample_agreement.empty:
        return
    style(cfg)
    comps = list(dict.fromkeys(sample_agreement["comparison"]))
    nrows, ncols = theme.grid_shape(len(comps), max_cols=2)
    fig, axes = plt.subplots(
        nrows, ncols, squeeze=False, sharey=True,
        figsize=(theme.width(cfg), 3.0 * nrows))
    flat = axes.ravel()
    for k, comp in enumerate(comps):
        ax = flat[k]
        sub = sample_agreement[sample_agreement["comparison"] == comp]
        sns.violinplot(data=sub, x="sample", y="exact_agreement",
                       hue="ptype", inner="box", cut=0, saturation=1,
                       palette=theme.PTYPE_COLORS, ax=ax, legend=(k == 0))
        ax.set_title(comp, fontsize=7.5)
        ax.set_xlabel("")
        ax.set_ylabel("exact agreement")
        ax.tick_params(axis="x", labelrotation=90, labelsize=4.5)
        ax.set_ylim(-0.03, 1.06)
        theme.panel_label(ax, theme.panel_letter(k), cfg)
        if k == 0:
            ax.legend(fontsize=6, ncol=2, loc="lower left")
    for ax in flat[len(comps):]:
        ax.set_visible(False)
    theme.finalize(fig, "tight")
    save(fig, cfg, f"{prefix}_sample_agreement")


def plot_frequency(freq: pd.DataFrame, cfg: dict[str, Any],
                   prefix: str = "C") -> None:
    if freq.empty:
        return
    style(cfg)
    sub = freq.copy()
    comps = list(dict.fromkeys(sub["comparison"]))
    ptypes = ["Amp", "Del"]
    fig, axes = plt.subplots(
        len(comps), len(ptypes), squeeze=False, sharex=True, sharey=True,
        figsize=theme.figsize(cfg, nrows=len(comps), panel_h=2.6))
    for i, comp in enumerate(comps):
        la, lb, _, _ = _sides(cfg, comp)
        for j, ptype in enumerate(ptypes):
            ax = axes[i][j]
            s = sub[(sub["comparison"] == comp) & (sub["ptype"] == ptype)]
            if s.empty:
                ax.set_visible(False)
                continue
            ax.scatter(s["a_freq"], s["b_freq"], s=16, alpha=0.8,
                       c=_pc(ptype), edgecolors="none")
            ax.plot([0, 1], [0, 1], "k--", lw=0.7, alpha=0.5)
            if len(s) >= 3 and np.std(s["a_freq"]) > 0:
                k, b = np.polyfit(s["a_freq"], s["b_freq"], 1)
                xs = np.linspace(0, 1, 50)
                ax.plot(xs, k * xs + b, color="#888888", lw=0.9)
            res_sp = stats.spearmanr(s["a_freq"], s["b_freq"])
            theme.annot_stat(ax, theme.fmt_ci(
                "\u03c1", float(res_sp.statistic), n=len(s),
                p=float(res_sp.pvalue)))
            ax.set_xlim(-0.02, 1.02)
            ax.set_ylim(-0.02, 1.02)
            ticks = [0.0, 0.25, 0.5, 0.75, 1.0]
            ax.set_xticks(ticks)
            ax.set_yticks(ticks)
            ax.tick_params(pad=2)
            if i == 0:
                ax.set_title(ptype)
                theme.panel_label(ax, theme.panel_letter(j), cfg)
            if j == 0:
                ax.set_ylabel(f"{comp}\n{lb}: frequency")
            if i == len(comps) - 1:
                ax.set_xlabel(f"{la}: frequency")
    theme.finalize(fig, "constrained", wspace=0.35, hspace=0.30)
    save(fig, cfg, f"{prefix}_frequency_scatter")


