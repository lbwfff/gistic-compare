"""Figures for module P3 (per-peak logistic).

(Split from the original single-module ``plotting``.)
"""

from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ._core import *  # noqa: F401,F403


def plot_p3_or_forest(models: pd.DataFrame, cfg: dict[str, Any]) -> None:
    """Crude and covariate-adjusted OR per peak (top 20 by crude p).

    Adjusted markers sit on the same row as the crude marker, so they are
    drawn *underneath*: crude keeps the highest zorder, a larger size and a
    white halo, otherwise it disappears completely behind the adjusted dots.
    """
    if models is None or models.empty or "crude_or" not in models.columns:
        return
    m = models[np.isfinite(models["crude_or"].to_numpy(float))].copy()
    if m.empty:
        return
    style(cfg)
    # top peaks by crude significance
    if "crude_p" in m.columns:
        m = m.sort_values("crude_p")
    m = m.head(20).iloc[::-1]
    adj_cols = [c for c in m.columns if c.startswith("adj_or_")
                and not c.startswith("adj_or_lo") and not c.startswith("adj_or_hi")]
    fig, ax = plt.subplots(figsize=theme.figsize(
        cfg, panel_h=0.28 * len(m) + 1.2))
    y = np.arange(len(m))
    cmap = plt.get_cmap("Blues")

    # 1) adjusted CIs and points, underneath
    for j, ac in enumerate(adj_cols):
        col = cmap(0.45 + 0.3 * j)
        lo_c, hi_c = f"adj_or_lo_{ac[7:]}", f"adj_or_hi_{ac[7:]}"
        if lo_c in m.columns and hi_c in m.columns:
            for i, (_, r) in enumerate(m.iterrows()):
                if np.isfinite(r.get(lo_c, np.nan)) and np.isfinite(r.get(hi_c, np.nan)):
                    ax.plot([r[lo_c], r[hi_c]], [i, i], color=col,
                            lw=0.8, alpha=0.75, zorder=2)
        valid = np.isfinite(m[ac].to_numpy(float))
        ax.scatter(m.loc[valid, ac], y[valid], s=30, color=col, zorder=3,
                   edgecolors="none",
                   label=ac.replace("adj_or_", "adj: "))

    # 2) crude CI and points on top so they stay visible
    if {"crude_or_lo", "crude_or_hi"} <= set(m.columns):
        for i, (_, r) in enumerate(m.iterrows()):
            if np.isfinite(r["crude_or_lo"]) and np.isfinite(r["crude_or_hi"]):
                ax.plot([r["crude_or_lo"], r["crude_or_hi"]], [i, i],
                        color="#222222", lw=1.1, zorder=4)
    ax.scatter(m["crude_or"], y, s=15, color="#222222", zorder=5,
               edgecolors="white", linewidths=0.4, label="crude", marker="o")

    ax.axvline(1.0, color="#888888", ls="--", lw=0.8, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels(
        [_clean_label(r.get("cytoband")) or _clean_label(r.get("peak_key"))
         for _, r in m.iterrows()], fontsize=6)
    ax.set_xscale("log")
    ax.set_xlabel("OR (cohort effect on alteration)")
    # plain, readable OR ticks: matplotlib's log locator only labels decades,
    # which collapses to a lone "1" on a narrow range like 0.4-6
    lo = [v for v in m.get("crude_or_lo", pd.Series(dtype=float)).tolist()
          if isinstance(v, float) and np.isfinite(v)]
    hi = [v for v in m.get("crude_or_hi", pd.Series(dtype=float)).tolist()
          if isinstance(v, float) and np.isfinite(v)]
    if lo and hi:
        ax.set_xlim(min(lo) * 0.85, max(hi) * 1.15)
    xlo, xhi = ax.get_xlim()
    cand = [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 7.0,
            10.0, 20.0, 30.0, 50.0]
    ticks = [c for c in cand if xlo <= c <= xhi]
    if ticks:
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"{c:g}" for c in ticks])
    ax.xaxis.set_minor_formatter(NullFormatter())
    h, lab = ax.get_legend_handles_labels()
    if h:
        ax.legend(h, lab, frameon=False, fontsize=6, loc="lower right")
    theme.annot_n(ax, f"n = {len(m)} peaks (top by crude p)",
                  loc="upper left")
    theme.panel_label(ax, theme.panel_letter(0), cfg)
    theme.finalize(fig, "constrained")
    save(fig, cfg, "P3_or_forest")


def plot_p3_qq(models: pd.DataFrame, cfg: dict[str, Any]) -> None:
    if models is None or models.empty or "crude_p" not in models.columns:
        return
    p = models["crude_p"].to_numpy(float)
    p = p[np.isfinite(p) & (p > 0)]
    if len(p) < 5:
        return
    style(cfg)
    n = len(p)
    exp = -np.log10(np.sort((np.arange(1, n + 1) - 0.5) / n))
    obs = -np.log10(np.sort(p))
    fig, ax = plt.subplots(figsize=theme.figsize(
        cfg, panel_h=3.4, total_w=4.6))
    ax.scatter(exp, obs, s=14, alpha=0.7, color="#0072B2",
               edgecolors="none")
    lim = max(exp.max(), obs.max()) * 1.05
    ax.plot([0, lim], [0, lim], color="#888888", ls="--", lw=0.8)
    ax.set_xlabel("expected −log10(p)")
    ax.set_ylabel("observed −log10(p)")
    # crude lambda
    med = float(np.median(p))
    lam = float(np.median(stats.chi2.isf(p, 1)) / stats.chi2.isf(0.5, 1))
    theme.annot_stat(ax, f"λ = {lam:.3f}, n = {n}", loc="lower right")
    theme.panel_label(ax, theme.panel_letter(0), cfg)
    theme.finalize(fig, "tight")
    save(fig, cfg, "P3_qq")


def plot_p3_volcano(models: pd.DataFrame, cfg: dict[str, Any]) -> None:
    """Crude OR vs −log10 FDR volcano; optional adjusted overlay."""
    if models is None or models.empty or "crude_or" not in models.columns:
        return
    m = models.copy()
    m = m[np.isfinite(m["crude_or"].to_numpy(float))
          & (m["crude_or"].to_numpy(float) > 0)]
    if m.empty or "crude_fdr" not in m.columns:
        return
    style(cfg)
    x = np.log(m["crude_or"].to_numpy(float))
    y = -np.log10(np.clip(m["crude_fdr"].to_numpy(float), 1e-300, None))
    adj_cols = [c for c in m.columns
                if c.startswith("adj_or_")
                and not c.startswith(("adj_or_lo", "adj_or_hi"))]
    n_panels = 1 + len(adj_cols)
    # 1 panel → single column; ≥2 panels → full width with taller rows
    if n_panels == 1:
        fw, ph = 4.6, 3.4
    else:
        fw, ph = theme.FULL_WIDTH, 3.4
    fig, axes = plt.subplots(1, n_panels,
                             figsize=theme.figsize(cfg, panel_h=ph,
                                                   total_w=fw),
                             squeeze=False)
    axes = axes[0]
    panels = [("crude", x, y)] + [
        (c.replace("adj_or_", "adj: "),
         np.log(m[c].to_numpy(float)),
         -np.log10(np.clip(
             m[c.replace("adj_or_", "adj_fdr_")].to_numpy(float)
             if c.replace("adj_or_", "adj_fdr_") in m.columns
             else np.full(len(m), np.nan), 1e-300, None)))
        for c in adj_cols if c in m.columns
    ]
    for ax, (name, xv, yv) in zip(axes, panels):
        for ptype in ("Amp", "Del"):
            mask = (m["ptype"] == ptype).to_numpy() if "ptype" in m.columns \
                else np.ones(len(m), bool)
            ok = mask & np.isfinite(xv) & np.isfinite(yv)
            if not ok.any():
                continue
            ax.scatter(xv[ok], yv[ok], s=14, alpha=0.75,
                       color=_pc(ptype), label=ptype, edgecolors="none")
        ax.axhline(-np.log10(0.05), color="#888888", ls=":", lw=0.7)
        ax.axvline(0, color="#888888", ls="--", lw=0.7)
        ax.set_title(name, fontsize=7.5)
        ax.set_xlabel("log OR (cohort effect)")
        if ax is axes[0]:
            ax.set_ylabel("−log10 FDR")
            ax.legend(frameon=False, fontsize=6)
        n_sig = int(np.nansum(yv >= -np.log10(0.05)))
        theme.annot_stat(ax, f"n = {len(m)}, sig = {n_sig}", loc="upper right")
        ax.tick_params(labelsize=6)
    for k, ax in enumerate(axes):
        theme.panel_label(ax, theme.panel_letter(k), cfg)
    theme.finalize(fig, "constrained", wspace=0.3)
    save(fig, cfg, "P3_volcano")


def plot_p3_or_shift(models: pd.DataFrame,
                     cov_summary: pd.DataFrame,
                     cfg: dict[str, Any]) -> None:
    """Per-covariate ΔlogOR (adj − crude) distribution + summary."""
    if models is None or models.empty:
        return
    m = models.copy()
    if "crude_or" not in m.columns:
        return
    adj_cols = [c for c in m.columns
                if c.startswith("adj_or_")
                and not c.startswith(("adj_or_lo", "adj_or_hi"))]
    if not adj_cols:
        return
    style(cfg)
    crude = np.log(m["crude_or"].to_numpy(float))
    n = len(adj_cols)
    fig, axes = plt.subplots(1, n, figsize=theme.figsize(
        cfg, panel_h=3.0, total_w=3.4 * n), squeeze=False)
    axes = axes[0]
    for ax, ac in zip(axes, adj_cols):
        adj = np.log(m[ac].to_numpy(float))
        d = adj - crude
        ok = np.isfinite(d)
        ax.hist(d[ok], bins=30, color="#0072B2", alpha=0.8,
                edgecolor="white", linewidth=0.4)
        ax.axvline(0, color="#333333", ls="--", lw=0.8)
        med = float(np.median(d[ok])) if ok.any() else np.nan
        flips = int(np.sum(np.sign(crude[ok]) != np.sign(adj[ok])))
        ax.set_title(ac.replace("adj_or_", ""), fontsize=7.5)
        ax.set_xlabel("Δ log OR (adj − crude)")
        if ax is axes[0]:
            ax.set_ylabel("n peaks")
        theme.annot_stat(
            ax, f"med = {med:+.3f}, flips = {flips}, n = {int(ok.sum())}",
            loc="upper right", fontsize=5.5)
        ax.tick_params(labelsize=6)
    for k, ax in enumerate(axes):
        theme.panel_label(ax, theme.panel_letter(k), cfg)
    theme.finalize(fig, "constrained", wspace=0.3)
    save(fig, cfg, "P3_or_shift")


def figures_p3(res: dict, cfg) -> None:
    if res.get("skipped"):
        return
    plot_p3_or_forest(res.get("models"), cfg)
    plot_p3_qq(res.get("models"), cfg)
    plot_p3_volcano(res.get("models"), cfg)
    plot_p3_or_shift(res.get("models"), res.get("covariate_summary"), cfg)


