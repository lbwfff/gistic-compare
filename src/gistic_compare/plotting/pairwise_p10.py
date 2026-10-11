"""Figures for module P10 (burden-controlled permutation on peaks)."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ._core import *  # noqa: F401,F403


def plot_p10_volcano(peaks: pd.DataFrame, cfg, prefix: str = "P10") -> None:
    """Delta-frequency vs -log10(permutation p); burden-inflated peaks get
    an orange ring (naive-significant but **not** permutation-robust — the
    signature a global burden difference would manufacture)."""
    if peaks is None or peaks.empty:
        return
    style(cfg)
    fig, axes = plt.subplots(1, 2, figsize=theme.figsize(
        cfg, panel_h=2.8, total_w=7.0))
    for ax, direction in zip(axes, ("Amp", "Del")):
        d = peaks[(peaks["ptype"] == direction)
                  & np.isfinite(peaks["p_perm"])]
        if d.empty:
            ax.set_visible(False)
            continue
        y = -np.log10(np.clip(d["p_perm"].to_numpy(float), 1e-300, 1.0))
        x = d["delta_freq"].to_numpy(float)
        col = theme.ptype_color(direction)
        ns = d["q_perm"] >= 0.05
        ax.scatter(x[ns.to_numpy()], y[ns.to_numpy()], s=7, color="#CFCFCF",
                   edgecolor="none")
        sg = ~ns.to_numpy()
        ax.scatter(x[sg], y[sg], s=9, color=col, edgecolor="none")
        bi = d["burden_flag"].to_numpy(dtype=bool)
        if bi.any():
            ax.scatter(x[bi], y[bi], s=26, facecolors="none",
                       edgecolors=theme.OKABE["orange"], linewidths=0.9,
                       zorder=5)
        la, lb, _, _ = _sides(cfg, d["comparison"].iloc[0])
        ax.set_xlabel(f"\u0394freq [{lb} \u2212 {la}]")
        ax.set_ylabel("\u2212log10 p (permutation)")
        n_sig = int((~ns).sum())
        n_bi = int(bi.sum())
        ax.set_title(f"{direction}  (perm q<0.05: n={n_sig}; "
                     f"burden-inflated: {n_bi})", fontsize=7)
        ax.axvline(0, color="#777777", lw=0.6, ls=":")
        ax.grid(alpha=0.2)
    fig.suptitle(f"{peaks['comparison'].iloc[0]}: peaks, burden-controlled",
                 fontsize=8)
    if "burden_flag" in peaks.columns and peaks["burden_flag"].any():
        theme.annot_n(ax, "\u25cb = naive p<0.05 but not permutation-robust "
                      "(burden-inflated)",
                      loc="upper left", fontsize=5.0)
    theme.finalize(fig, "tight")
    save(fig, cfg, f"{prefix}_volcano")


def plot_p10_naive_vs_perm(peaks: pd.DataFrame, cfg,
                           prefix: str = "P10") -> None:
    """Naive vs permutation p per peak: points below the diagonal (x > y,
    i.e. smaller naive p) were inflated by the global burden difference;
    points above it are conservative under margin control."""
    if peaks is None or peaks.empty:
        return
    d = peaks[np.isfinite(peaks["p_naive"]) & np.isfinite(peaks["p_perm"])]
    if d.empty:
        return
    style(cfg)
    comp = d["comparison"].iloc[0]
    fig, ax = plt.subplots(figsize=theme.figsize(cfg, panel_h=2.9,
                                                 total_w=3.4))
    xn = -np.log10(np.clip(d["p_naive"].to_numpy(float), 1e-300, 1.0))
    yp = -np.log10(np.clip(d["p_perm"].to_numpy(float), 1e-300, 1.0))
    bi = d["burden_flag"].to_numpy(dtype=bool)
    ax.scatter(xn[~bi], yp[~bi], s=8, color="#8FA8C8", alpha=0.8,
               edgecolor="none")
    ax.scatter(xn[bi], yp[bi], s=12, color=theme.OKABE["orange"],
               edgecolor="none", label="burden-inflated")
    lim = max(0.5, float(np.nanmax([xn.max(), yp.max()])) * 1.05)
    ax.plot([0, lim], [0, lim], color="#999999", lw=0.7, ls="--")
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel("\u2212log10 p (naive, no margin control)", fontsize=7)
    ax.set_ylabel("\u2212log10 p (permutation)", fontsize=7)
    ax.set_title(f"{comp}: burden inflation per peak", fontsize=7)
    ax.legend(fontsize=5.5, loc="upper left")
    ax.grid(alpha=0.2)
    theme.finalize(fig, "tight")
    save(fig, cfg, f"{prefix}_naive_vs_perm")


def plot_p10_zhist(res: dict, cfg, prefix: str = "P10") -> None:
    """Permutation Z distribution vs observed Z (Yuan Fig 4C style)."""
    store = res.get("Z_perm") or {}
    peaks = res.get("peaks")
    if not store or peaks is None or peaks.empty:
        return
    style(cfg)
    comps = [c for c in dict.fromkeys(peaks["comparison"]) if c in store]
    nrows, ncols = theme.grid_shape(len(comps), max_cols=2)
    fig, axes = plt.subplots(nrows, ncols, figsize=theme.figsize(
        cfg, panel_h=2.4), squeeze=False)
    axes = axes.ravel()
    for k, comp in enumerate(comps):
        ax = axes[k]
        Z = np.asarray(store[comp], dtype=float).ravel()
        Z = Z[np.isfinite(Z)]
        if len(Z):
            ax.hist(Z, bins=41, density=True, color="#D9D9D9",
                    edgecolor="white", linewidth=0.3)
            from scipy import stats as _s
            xs = np.linspace(-4.5, 4.5, 200)
            ax.plot(xs, _s.norm.pdf(xs), color="#555555", lw=0.9,
                    label="N(0,1)")
        zo = peaks.loc[peaks["comparison"] == comp, "Z_obs"]
        zo = zo[np.isfinite(zo)].to_numpy(float)
        ax.scatter(zo, np.full(len(zo), ax.get_ylim()[1] * 0.02), s=5,
                   color=theme.OKABE["vermillion"], marker="|", zorder=5)
        ax.set_xlim(-4.5, 4.5)
        ax.set_xlabel("Z (per peak)", fontsize=7)
        ax.set_ylabel("density", fontsize=7)
        ax.set_title(f"{comp}", fontsize=7)
        ax.legend(fontsize=5.5, loc="upper left")
        ax.tick_params(labelsize=6.5)
    for k in range(len(comps), len(axes)):
        axes[k].set_visible(False)
    theme.finalize(fig, "tight")
    save(fig, cfg, f"{prefix}_zhist")


def figures_p10(res: dict, cfg) -> None:
    plot_p10_volcano(res.get("peaks"), cfg)
    plot_p10_naive_vs_perm(res.get("peaks"), cfg)
    plot_p10_zhist(res, cfg)
