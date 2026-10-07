"""Figures for module G (reference panel bench).

(Split from the original single-module ``plotting``.)
"""

from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ._core import *  # noqa: F401,F403

# ---------------------------------------------------------------------------
# G: paper WE benchmark
# ---------------------------------------------------------------------------

def plot_paper_genome(genome: pd.DataFrame, cfg: dict[str, Any]) -> None:
    """Bar chart: amp/del genome % per experiment vs paper reference lines."""
    if genome.empty:
        return
    style(cfg)
    control = cfg.get("control")
    genome = genome[genome["cohort"] != control] if control is not None else genome
    labels = (list(genome["label"]) if "label" in genome.columns
              else list(genome["cohort"]))
    tick_labels = []
    for lab, (_, r) in zip(labels, genome.iterrows()):
        na = r.get("n_amp_peaks", np.nan)
        nd = r.get("n_del_peaks", np.nan)
        if np.isfinite(na) and np.isfinite(nd):
            tick_labels.append(f"{lab}\n({int(na)}A/{int(nd)}D peaks)")
        else:
            tick_labels.append(str(lab))
    x = np.arange(len(genome))
    w = 0.35
    amp_c = theme.PTYPE_COLORS["Amp"]
    del_c = theme.PTYPE_COLORS["Del"]
    amp_ref = float(genome["paper_amp_genome_pct"].iloc[0])
    del_ref = float(genome["paper_del_genome_pct"].iloc[0])
    fig, ax = plt.subplots(figsize=theme.figsize(cfg, panel_h=3.4,
                                                 total_w=4.4))
    ax.bar(x - w / 2, genome["amp_genome_pct"], w, color=amp_c,
           label="Amp observed")
    ax.bar(x + w / 2, genome["del_genome_pct"], w, color=del_c,
           label="Del observed")
    ax.axhline(amp_ref, color=amp_c, ls="--", lw=1.2,
               label=f"paper Amp {amp_ref:g}%")
    ax.axhline(del_ref, color=del_c, ls="--", lw=1.2,
               label=f"paper Del {del_ref:g}%")
    ax.set_xticks(x)
    ax.set_xticklabels(tick_labels, rotation=20, ha="right")
    ax.set_ylabel("% of genome (sig peak union)")
    ax.legend(fontsize=6, loc="upper left", ncol=2)
    ymax = max(genome["del_genome_pct"].max(),
               genome["amp_genome_pct"].max(), del_ref, amp_ref)
    ax.set_ylim(0, ymax * 1.12)
    theme.finalize(fig, "tight")
    save(fig, cfg, "G_paper_genome_pct")


def plot_paper_cytoband_heatmap(cyto: pd.DataFrame, cfg: dict[str, Any]) -> None:
    """Presence of paper-claimed regions across experiment cohorts (no control)."""
    if cyto.empty:
        return
    style(cfg)
    control = cfg.get("control")
    cohorts = [c for c in cfg["cohorts"] if c != control]
    labels = []
    for c in cohorts:
        lab = theme.cohort_label(cfg, c)
        nn = _n(cfg, c)
        labels.append(f"{lab}\n(n = {nn})" if nn else lab)
    piv = (cyto.pivot_table(index=["side", "prefix"], columns="cohort",
                            values="present", aggfunc="first")
           .reindex(columns=cohorts)
           .fillna(False).astype(int))
    order = cyto.drop_duplicates(["side", "prefix"])[["side", "prefix"]]
    keys = list(order.itertuples(index=False, name=None))
    idx = [k for k in keys if k in piv.index]
    piv = piv.reindex(idx)
    piv.columns = labels
    ylabels = [f"{s} {p}" for s, p in piv.index]
    # presence coded by side: 0 absent, 1 Amp (orange), 2 Del (blue)
    mat = piv.copy()
    for row_idx in mat.index:
        side = row_idx[0]
        mat.loc[row_idx] = [
            0 if not v else (2 if side == "Del" else 1)
            for v in piv.loc[row_idx]
        ]
    from matplotlib.colors import BoundaryNorm, ListedColormap
    cmap = ListedColormap(["#F2F2F2", _pc("Amp"), _pc("Del")])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)
    fig, ax = plt.subplots(
        figsize=(5.0, 0.22 * len(piv) + 1.0))
    sns.heatmap(mat, cmap=cmap, norm=norm, cbar=False,
                linewidths=0.5, linecolor="#FFFFFF", ax=ax,
                xticklabels=True, yticklabels=ylabels)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha="right",
                       fontsize=6.5)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="y", labelsize=6)
    if "match_method" in cyto.columns:
        method = (cyto.pivot_table(index=["side", "prefix"], columns="cohort",
                                   values="match_method", aggfunc="first")
                  .reindex(index=piv.index, columns=cohorts))
        method.columns = labels
        for i in range(piv.shape[0]):
            for j in range(piv.shape[1]):
                if piv.iloc[i, j] >= 0.5:
                    mth = str(method.iloc[i, j])
                    if mth == "window":
                        ax.text(j + 0.5, i + 0.5, "W", ha="center",
                                va="center", fontsize=6, color="white")
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    ax.legend(
        handles=[
            Patch(facecolor="#F2F2F2", edgecolor="#999999", label="absent"),
            Patch(facecolor=_pc("Amp"), edgecolor="none", label="Amp peak"),
            Patch(facecolor=_pc("Del"), edgecolor="none", label="Del peak"),
            Line2D([], [], marker="$W$", color="none",
                   markeredgecolor="#333333", markersize=7,
                   label="W = window match"),
        ],
        loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False,
        fontsize=6,
    )
    fig.subplots_adjust(right=0.74, top=0.98)
    theme.finalize(fig, "none")
    save(fig, cfg, "G_paper_cytoband")


def plot_paper_driver_heatmap(driver: pd.DataFrame, cfg: dict[str, Any]) -> None:
    """Paper-named driver genes: presence per experiment (Amp/Del/both)."""
    if driver.empty:
        return
    style(cfg)
    control = cfg.get("control")
    cohorts = [c for c in cfg["cohorts"] if c != control]
    labels = []
    for c in cohorts:
        lab = theme.cohort_label(cfg, c)
        nn = _n(cfg, c)
        labels.append(f"{lab}\n(n = {nn})" if nn else lab)
    genes = list(dict.fromkeys(driver["gene"].tolist()))

    code = {"": 0, "Amp": 1, "Del": 2, "both": 3}
    mat = pd.DataFrame(0, index=genes, columns=cohorts, dtype=int)
    for _, r in driver.iterrows():
        if r["gene"] in mat.index and r["cohort"] in mat.columns:
            mat.loc[r["gene"], r["cohort"]] = code.get(str(r["amp_or_del"]), 0)
    mat.columns = labels

    from matplotlib.colors import BoundaryNorm, ListedColormap
    absent_c, amp_c, del_c = "#FFFFFF", _pc("Amp"), _pc("Del")
    both_c = theme.OKABE["reddish_purple"]
    cmap = ListedColormap([absent_c, amp_c, del_c, both_c])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cmap.N)

    fig, ax = plt.subplots(
        figsize=(theme.COL_WIDTH + 0.8, 0.26 * len(genes) + 1.0))
    sns.heatmap(mat, cmap=cmap, norm=norm, cbar=False,
                linewidths=0.5, linecolor="#CCCCCC", ax=ax,
                xticklabels=True, yticklabels=genes)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha="right",
                       fontsize=6.5)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="y", labelsize=6)
    from matplotlib.patches import Patch
    ax.legend(
        handles=[
            Patch(facecolor=absent_c, edgecolor="#999999", label="absent"),
            Patch(facecolor=amp_c, edgecolor="#999999", label="Amp"),
            Patch(facecolor=del_c, edgecolor="#999999", label="Del"),
            Patch(facecolor=both_c, edgecolor="#999999", label="both"),
        ],
        loc="upper left", fontsize=6, frameon=False,
        bbox_to_anchor=(1.02, 1.0),
    )
    fig.subplots_adjust(right=0.72, top=0.98)
    theme.finalize(fig, "none")
    save(fig, cfg, "G_paper_driver")


def plot_paper_top_freq(freq: pd.DataFrame, cfg: dict[str, Any]) -> None:
    """Grouped bars: observed max frequency of top paper regions vs threshold."""
    if freq.empty:
        return
    style(cfg)
    control = cfg.get("control")
    if control is not None:
        freq = freq[freq["cohort"] != control]
    regions = list(dict.fromkeys(zip(freq["side"], freq["region"])))
    cohorts = list(dict.fromkeys(freq["cohort"]))
    n_r, n_c = len(regions), len(cohorts)
    fig, ax = plt.subplots(
        figsize=(max(theme.COL_WIDTH, 0.85 * n_r + 1.4), 3.4))
    x = np.arange(n_r)
    w = 0.8 / max(n_c, 1)
    pal = _palette(cfg)
    for ci, cid in enumerate(cohorts):
        vals, thr = [], []
        for side, reg in regions:
            row = freq[(freq.cohort == cid) & (freq.side == side)
                       & (freq.region == reg)]
            if len(row) and np.isfinite(row["observed_freq_pct"].iloc[0]):
                vals.append(float(row["observed_freq_pct"].iloc[0]))
            else:
                vals.append(0.0)
            thr.append(float(row["paper_min_freq_pct"].iloc[0]) if len(row) else np.nan)
        offset = (ci - (n_c - 1) / 2) * w
        nn = _n(cfg, cid)
        lab = theme.cohort_label(cfg, cid)
        ax.bar(x + offset, vals, w * 0.9, color=pal.get(cid, "#333333"),
               label=f"{lab} (n = {nn})" if nn else lab)
        if ci == 0:
            for xi, t in enumerate(thr):
                if np.isfinite(t):
                    ax.hlines(t, xi - 0.45, xi + 0.45, colors="black",
                              linestyles="--", lw=1.0, zorder=5)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{s}\n{r}" for s, r in regions],
                       rotation=40, ha="right", fontsize=6.5)
    ax.set_ylabel("max frequency among samples (%)")
    top = max(ax.get_ylim()[1], 60.0) * 1.06
    ax.set_ylim(0, top)
    ax.legend(fontsize=6.5, loc="upper right", ncol=1, frameon=False)
    theme.finalize(fig, "tight")
    save(fig, cfg, "G_paper_top_freq")


