"""Figures for module P8 (multivariate cohort separation)."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ._core import *  # noqa: F401,F403


def plot_p8_pca(pca: pd.DataFrame, cfg, sep: pd.DataFrame | None = None
                ) -> None:
    """PC1 vs PC2 ordination coloured by cohort (paired samples joined).

    Panels are square (equal aspect) and annotated with the PERMANOVA and
    CV-AUC statistics: the ordination is *unsupervised* — the first two PCs
    maximise total variance, not cohort contrast — so a significant
    classifier AUC with little visual separation is expected whenever the
    cohort effect is small relative to overall CN heterogeneity.
    """
    if pca is None or pca.empty:
        return
    style(cfg)
    comps = list(dict.fromkeys(pca["comparison"]))
    nrows, ncols = theme.grid_shape(len(comps), max_cols=2)
    panel = 3.4
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(panel * ncols, panel * nrows),
                             squeeze=False)
    axes = axes.ravel()
    pal = _palette(cfg)
    for k, comp in enumerate(comps):
        ax = axes[k]
        d = pca[pca["comparison"] == comp]
        la, lb, _, _ = _sides(cfg, comp)
        for cid in (d["cohort"].unique()):
            sub = d[d["cohort"] == cid]
            ax.scatter(sub["pc1"], sub["pc2"], s=10, alpha=0.8,
                       color=pal.get(cid, "#888888"),
                       label=theme.cohort_label(cfg, cid), edgecolor="none")
        if "pair_id" in d.columns and d["kind"].iloc[0] == "paired":
            for _, g in d.groupby("pair_id"):
                if len(g) == 2:
                    ax.plot(g["pc1"], g["pc2"], color="#BBBBBB", lw=0.5,
                            zorder=1)
        v1 = d["var_pc1"].iloc[0] * 100 if np.isfinite(d["var_pc1"].iloc[0]) \
            else 0
        v2 = d["var_pc2"].iloc[0] * 100 if np.isfinite(d["var_pc2"].iloc[0]) \
            else 0
        ax.set_xlabel(f"PC1 ({v1:.1f}%)", fontsize=7)
        ax.set_ylabel(f"PC2 ({v2:.1f}%)", fontsize=7)
        ax.set_title(f"{comp} ({_n_pair(cfg, comp)})", fontsize=7)
        ax.legend(fontsize=5.5, handletextpad=0.2, loc="upper left")
        if sep is not None and len(sep):
            s = sep[sep["comparison"] == comp]
            if len(s):
                r = s.iloc[0]
                theme.annot_stat(
                    ax,
                    f"PERMANOVA R\u00b2 = {r['permanova_R2']:.3f}, "
                    f"p = {r['permanova_p']:.3f}\n"
                    f"CV-AUC = {r['auc_oof']:.2f} "
                    f"(perm p = {r['auc_perm_p']:.3f})",
                    loc="lower right", fontsize=5.5)
        ax.grid(alpha=0.2)
        ax.set_aspect("equal")
    for k in range(len(comps), len(axes)):
        axes[k].set_visible(False)
    theme.finalize(fig, "tight")
    save(fig, cfg, "P8_pca")


def plot_p8_auc_perm(res: dict, cfg) -> None:
    """Permutation null of the cross-validated AUC, per comparison."""
    store = res.get("auc_perm") or {}
    sep = res.get("separation")
    if not store or sep is None or sep.empty:
        return
    style(cfg)
    comps = [c for c in dict.fromkeys(sep["comparison"]) if c in store]
    if not comps:
        return
    nrows, ncols = theme.grid_shape(len(comps), max_cols=2)
    fig, axes = plt.subplots(nrows, ncols, figsize=theme.figsize(
        cfg, panel_h=2.4), squeeze=False)
    axes = axes.ravel()
    for k, comp in enumerate(comps):
        ax = axes[k]
        aucs = np.asarray(store[comp], dtype=float)
        aucs = aucs[np.isfinite(aucs)]
        obs = sep.loc[sep["comparison"] == comp, "auc_oof"]
        obs = float(obs.iloc[0]) if len(obs) and np.isfinite(obs.iloc[0]) \
            else np.nan
        if len(aucs):
            ax.hist(aucs, bins=16, color="#D9D9D9", edgecolor="white",
                    linewidth=0.4)
        if np.isfinite(obs):
            ax.axvline(obs, color=theme.OKABE["vermillion"], lw=1.2)
            theme.annot_stat(ax, f"observed AUC = {obs:.2f}",
                             loc="upper left")
        ax.set_xlabel("cross-validated AUC (permuted labels)", fontsize=7)
        ax.set_ylabel("permutations", fontsize=7)
        ax.set_title(comp, fontsize=7)
        ax.tick_params(labelsize=6.5)
    for k in range(len(comps), len(axes)):
        axes[k].set_visible(False)
    theme.finalize(fig, "tight")
    save(fig, cfg, "P8_auc_permutation")


def _short_feat(name: str) -> str:
    """'chr14:103000000-104000000' -> 'chr14:103.0-104.0Mb'."""
    try:
        chrom, span = name.split(":")
        s, e = span.split("-")
        return f"{chrom}:{int(s)/1e6:.1f}-{int(e)/1e6:.1f}Mb"
    except ValueError:
        return name


def plot_p8_stability(stab: pd.DataFrame, cfg, n_show: int = 12) -> None:
    """Selection frequency of the most stable features, per comparison.

    Bar length = in how many of the bootstrap refits the feature landed in
    the top-k by |coefficient| (1.0 = always).  Bar colour = sign of the
    mean coefficient, i.e. which cohort the feature pushes the prediction
    towards.
    """
    if stab is None or stab.empty:
        return
    style(cfg)
    comps = list(dict.fromkeys(stab["comparison"]))
    nrows, ncols = theme.grid_shape(len(comps), max_cols=1)
    panel_h = min(0.32 * n_show + 1.1, 5.6)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.9, panel_h * nrows),
                             squeeze=False)
    axes = axes.ravel()
    for k, comp in enumerate(comps):
        ax = axes[k]
        d = stab[stab["comparison"] == comp].sort_values("rank")
        d = d.head(n_show).iloc[::-1]
        y = np.arange(len(d))
        colors = [theme.OKABE["vermillion"] if c > 0 else theme.OKABE["blue"]
                  for c in d["mean_coef"]]
        ax.barh(y, d["selection_freq"], color=colors, edgecolor="white",
                linewidth=0.3)
        ax.set_yticks(y)
        labels = [f"{a}  ({_short_feat(f)})" if a else _short_feat(f)
                  for f, a in zip(d["feature"], d["annotation"])]
        ax.set_yticklabels(labels, fontsize=5.8)
        ax.set_xlabel("selection frequency (bootstrap refits)", fontsize=7)
        la, lb, _, _ = _sides(cfg, comp)
        ax.set_title(f"{comp} \u2014 top {len(d)} stable features", fontsize=7)
        # direction legend: coefficient sign = predicted cohort
        import matplotlib.patches as mpatches
        ax.legend(handles=[
            mpatches.Patch(color=theme.OKABE["vermillion"],
                           label=f"higher in {lb}"),
            mpatches.Patch(color=theme.OKABE["blue"],
                           label=f"higher in {la}")],
            fontsize=5.5, loc="lower right", handletextpad=0.4)
        ax.tick_params(labelsize=6)
        ax.grid(axis="x", alpha=0.2)
    for k in range(len(comps), len(axes)):
        axes[k].set_visible(False)
    theme.finalize(fig, "tight")
    save(fig, cfg, "P8_stability")


def figures_p8(res: dict, cfg) -> None:
    plot_p8_pca(res.get("pca"), cfg, res.get("separation"))
    plot_p8_auc_perm(res, cfg)
    plot_p8_stability(res.get("stability"), cfg)
