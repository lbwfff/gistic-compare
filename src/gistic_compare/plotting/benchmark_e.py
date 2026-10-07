"""Figures for module E (FGA QC).

(Split from the original single-module ``plotting``.)
"""

from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ._core import *  # noqa: F401,F403

# ---------------------------------------------------------------------------
# E: FGA QC
# ---------------------------------------------------------------------------

def plot_fga(qc: pd.DataFrame, cfg: dict[str, Any],
             stats_df: pd.DataFrame | None = None) -> None:
    if qc.empty:
        return
    style(cfg)
    order = list(cfg["cohorts"].keys())
    pal = _palette(cfg)
    metrics = [("fga", "FGA (|log2| \u2265 0.2)"),
               ("n_segments", "segments / sample"),
               ("n_altered_peaks", "altered peaks / sample")]
    fig, axes = plt.subplots(1, 3, figsize=theme.figsize(cfg,
                                                         panel_h=3.8))
    for k, (ax, (metric, ylabel)) in enumerate(zip(axes, metrics)):
        if metric not in qc.columns:
            ax.set_visible(False)
            continue
        sub = qc.dropna(subset=[metric])
        present = [c for c in order if c in set(sub["cohort"])]
        sns.violinplot(data=sub, x="cohort", y=metric, order=present,
                       hue="cohort", palette=pal, inner="box", cut=0,
                       saturation=1, legend=False, ax=ax)
        theme.strip_overlay(ax, sub, "cohort", metric, order=present,
                            hue="cohort", palette=pal, dodge=True)
        ax.set_xlabel("")
        ax.set_ylabel(ylabel)
        # n travels with the tick labels: no collision with the distribution
        ax.set_xticks(range(len(present)))
        ax.set_xticklabels(
            [f"{theme.cohort_label(cfg, c)}\n"
             f"(n = {int((sub['cohort'] == c).sum())})" for c in present],
            rotation=0, ha="center", fontsize=6)
        # significance ladder: BH-adjusted pairwise Mann-Whitney
        if stats_df is not None and not stats_df.empty:
            mrows = stats_df[(stats_df["metric"] == metric)
                             & (stats_df["test"] == "mannwhitneyu")].copy()
            if not mrows.empty:
                mrows["_q"] = theme.bh_adjust(mrows["pvalue"].to_numpy(float))
                pairs = []
                for _, r in mrows.iterrows():
                    a, b = theme.comp_sides(cfg, r["comparison"])
                    if a in present and b in present:
                        pairs.append((min(present.index(a), present.index(b)),
                                      max(present.index(a), present.index(b)),
                                      float(r["_q"])))
                pairs.sort(key=lambda t: (t[1] - t[0], t[0]))
                ymax = float(np.nanmax(sub[metric].to_numpy(float)))
                ymin = float(np.nanmin(sub[metric].to_numpy(float)))
                span = max(ymax - ymin, 1e-9)
                y = ymax + 0.08 * span
                n_br = 0
                for lo_i, hi_i, q in pairs:
                    if hi_i - lo_i < 1:
                        continue
                    theme.sig_bracket(ax, lo_i, hi_i, y, q, h=0.025 * span)
                    y += 0.10 * span
                    n_br += 1
                if n_br:
                    ax.set_ylim(0, y)
            kw = stats_df[(stats_df["metric"] == metric)
                          & (stats_df["test"] == "kruskal_wallis")]
            if not kw.empty:
                # as an axes title: an in-axes corner would sit under the
                # significance ladder that spans the full panel width
                ax.set_title("Kruskal-Wallis " + theme.fmt_p(
                    float(kw["pvalue"].iloc[0])), fontsize=6.5)
        ax.set_ylim(bottom=0)
        ax.tick_params(axis="x", pad=2)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
        theme.panel_label(ax, theme.panel_letter(k), cfg)
    theme.finalize(fig, "tight")
    save(fig, cfg, "E_fga_qc")


def plot_fga_vs_peaks(qc: pd.DataFrame, cfg: dict[str, Any],
                      corr: pd.DataFrame | None = None) -> None:
    if qc.empty or "n_altered_peaks" not in qc.columns:
        return
    style(cfg)
    sub = qc.dropna(subset=["fga", "n_altered_peaks"])
    if sub.empty:
        return
    fig, ax = plt.subplots(figsize=theme.figsize(cfg, panel_h=3.4,
                                                 total_w=4.2))
    pal = _palette(cfg)
    rmap = {}
    if corr is not None and not corr.empty:
        for _, r in corr.iterrows():
            rmap[r["cohort"]] = (r["spearman_r"], r["pvalue"])
    for cid, g in sub.groupby("cohort"):
        rtxt = ""
        if cid in rmap:
            rr, _pp = rmap[cid]
            rtxt = f", \u03c1 = {rr:.2f}" if np.isfinite(rr) else ""
        ax.scatter(g["fga"], g["n_altered_peaks"], s=12, alpha=0.65,
                   label=(f"{theme.cohort_label(cfg, cid)} "
                          f"(n = {len(g)}{rtxt})"),
                   color=pal.get(cid, "#333333"), edgecolors="none")
    if len(sub) >= 3 and np.std(sub["fga"]) > 0:
        k, b = np.polyfit(sub["fga"], sub["n_altered_peaks"], 1)
        xs = np.linspace(sub["fga"].min(), sub["fga"].max(), 50)
        ax.plot(xs, k * xs + b, color="#555555", lw=1.0, ls="--")
    ax.set_xlabel("FGA")
    ax.set_ylabel("altered peaks / sample")
    ax.legend(fontsize=6, loc="upper left", ncol=1)
    theme.finalize(fig, "tight")
    save(fig, cfg, "E_fga_vs_peaks")


