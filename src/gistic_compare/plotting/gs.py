"""Figures for module GS (G-score differential).

(Split from the original single-module ``plotting``.)
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from ._core import *  # noqa: F401,F403

# ---------------------------------------------------------------------------
# GS: G-score differential (window / band)
# ---------------------------------------------------------------------------

def plot_gs_genome(df, cfg, direction, unit, prefix="GS"):
    if df is None or df.empty:
        return
    dcol, ycol = f"fdr_{direction.lower()}", f"dG_{direction.lower()}"
    if dcol not in df.columns:
        return
    d = df[np.isfinite(df[dcol])].copy()
    if d.empty:
        return
    style(cfg)
    d["x"] = _genome_x(d)
    d["sig"] = d[dcol] < 0.05
    col = theme.PTYPE_COLORS.get(direction, theme.OKABE["vermillion"])
    fig, ax = plt.subplots(figsize=theme.figsize(cfg, panel_h=2.2, total_w=7.0))
    _draw_chrom_guides(ax, d)
    ns, sg = d[~d["sig"]], d[d["sig"]]
    ax.scatter(ns["x"], ns[ycol], s=5, color="#CFCFCF", edgecolor="none",
               zorder=2)
    if not sg.empty:
        ax.scatter(sg["x"], sg[ycol], s=9, color=col, edgecolor="none",
                   zorder=3)
    ax.axhline(0, color="#777777", lw=0.6, ls=":", zorder=1)
    for _, r in d[d["paper_locus"] != ""].iterrows():
        ax.annotate(str(r["paper_locus"]), (float(r["x"]), float(r[ycol])),
                    fontsize=5, xytext=(0, 3), textcoords="offset points",
                    ha="center", color="#333333", zorder=4)
    bnds = _chrom_bounds(d)
    ax.set_xticks([b for b, _ in bnds])
    ax.set_xticklabels([c.replace("chr", "") for _, c in bnds], fontsize=5)
    ax.set_ylabel(f"\u0394G ({direction})  [A \u2212 B]")
    ax.set_title(
        f"{df['comparison'].iloc[0]}: G-score diff by {unit} \u2014 {direction}"
        f"  (coloured = BH FDR<0.05, n={int(d['sig'].sum())})", fontsize=7)
    ax.grid(axis="y", alpha=0.2)
    theme.finalize(fig, "tight")
    save(fig, cfg, f"{prefix}_{unit}_genome_{direction.lower()}")


def plot_gs_volcano(df, cfg, unit, prefix="GS"):
    if df is None or df.empty:
        return
    style(cfg)
    fig, axes = plt.subplots(1, 2, figsize=theme.figsize(
        cfg, panel_h=2.7, total_w=7.0))
    for ax, direction in zip(axes, ("Amp", "Del")):
        d = direction.lower()
        fcol = f"fdr_{d}"
        if fcol not in df.columns:
            ax.set_visible(False)
            continue
        s = df[np.isfinite(df[fcol])].copy()
        if s.empty:
            ax.set_visible(False)
            continue
        s["y"] = -np.log10(np.clip(s[fcol].astype(float), 1e-300, 1.0))
        sig = s[s[fcol] < 0.05]
        col = theme.PTYPE_COLORS.get(direction, theme.OKABE["vermillion"])
        ax.scatter(s["dG_" + d], s["y"], s=6, color="#CFCFCF", edgecolor="none")
        ax.scatter(sig["dG_" + d], sig["y"], s=9, color=col, edgecolor="none")
        ax.axvline(0, color="#777777", lw=0.6, ls=":")
        # one label per preset locus (its strongest point), collision-safe
        ann = (s[s["paper_locus"] != ""]
               .sort_values("y", ascending=False)
               .groupby("paper_locus", as_index=False).head(1))
        if not ann.empty:
            from ._core import _annotate_labeled_points
            _annotate_labeled_points(ax, ann["dG_" + d].to_numpy(dtype=float),
                                     ann["y"].to_numpy(dtype=float),
                                     ann["paper_locus"].tolist(),
                                     fontsize=5)
        ax.set_xlabel(f"\u0394G ({direction})  [A \u2212 B]")
        ax.set_ylabel("\u2212log10 FDR")
        ax.set_title(f"{direction}  (sig n={len(sig)})", fontsize=7)
        ax.grid(alpha=0.2)
    fig.suptitle(f"{df['comparison'].iloc[0]}: G-score differential by {unit}",
                 fontsize=8)
    theme.finalize(fig, "tight")
    save(fig, cfg, f"{prefix}_{unit}_volcano")


def figures_gs(res: dict, cfg) -> None:
    win, band = res.get("window"), res.get("band")
    for direction in ("Amp", "Del"):
        plot_gs_genome(win, cfg, direction, "window")
        plot_gs_genome(band, cfg, direction, "band")
    plot_gs_volcano(win, cfg, "window")
    plot_gs_volcano(band, cfg, "band")
