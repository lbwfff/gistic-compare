"""Figures for module D (gene sets).

(Split from the original single-module ``plotting``.)
"""

from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt
import pandas as pd

from ._core import *  # noqa: F401,F403

# ---------------------------------------------------------------------------
# D: gene sets
# ---------------------------------------------------------------------------

def plot_driver_heatmap(drivers: pd.DataFrame, cfg: dict[str, Any]) -> None:
    """Driver presence: Amp block (red) | Del block (blue); absent = white."""
    if drivers.empty:
        return
    style(cfg)
    det = drivers[drivers["detected_any"]].copy()
    if det.empty:
        det = drivers.copy()
    cohorts = list(cfg["cohorts"].keys())
    labels = [cfg["cohorts"][c].get("label", c) for c in cohorts]
    ns = [_n(cfg, c) for c in cohorts]
    xt = [f"{lab}\n(n = {n})" if n else lab for lab, n in zip(labels, ns)]
    mat_a = det.set_index("gene")[[f"{c}_amp" for c in cohorts]].astype(int)
    mat_d = det.set_index("gene")[[f"{c}_del" for c in cohorts]].astype(int)
    mat_a.columns = labels
    mat_d.columns = labels
    keep = (mat_a.sum(axis=1) + mat_d.sum(axis=1)) > 0
    mat_a = mat_a.loc[keep]
    mat_d = mat_d.loc[keep]
    order = sorted(
        mat_a.index,
        key=lambda g: (-(mat_a.loc[g].sum() + mat_d.loc[g].sum()), g),
    )
    mat_a = mat_a.loc[order]
    mat_d = mat_d.loc[order]
    if mat_a.empty:
        return

    mask_a = mat_a == 0
    mask_d = mat_d == 0

    fig = plt.figure(figsize=(theme.width(cfg),
                              0.24 * len(mat_a) + 1.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.5],
                          wspace=0.06, left=0.14, right=0.94,
                          top=0.92, bottom=0.08)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_d = fig.add_subplot(gs[0, 1])
    ax_leg = fig.add_subplot(gs[0, 2])
    ax_leg.axis("off")

    from matplotlib.colors import ListedColormap
    amp_c = _pc("Amp")
    del_c = _pc("Del")

    sns.heatmap(mat_a, mask=mask_a.to_numpy(),
                cmap=ListedColormap(["#FFFFFF", amp_c]), cbar=False,
                linewidths=0.5, linecolor="#DDDDDD", ax=ax_a, vmin=0, vmax=1,
                yticklabels=True, xticklabels=True)
    ax_a.text(0.0, 1.02, "Amplification", transform=ax_a.transAxes,
              color=amp_c, fontsize=7.5, ha="left", va="bottom")
    ax_a.set_ylabel("gene", fontsize=7.5)
    ax_a.tick_params(axis="y", labelsize=6)
    ax_a.set_xticklabels(xt, fontsize=6.5)
    ax_a.set_xlabel("")

    sns.heatmap(mat_d, mask=mask_d.to_numpy(),
                cmap=ListedColormap(["#FFFFFF", del_c]), cbar=False,
                linewidths=0.5, linecolor="#DDDDDD", ax=ax_d, vmin=0, vmax=1,
                yticklabels=False, xticklabels=True)
    ax_d.text(0.0, 1.02, "Deletion", transform=ax_d.transAxes,
              color=del_c, fontsize=7.5, ha="left", va="bottom")
    ax_d.set_ylabel("")
    ax_d.set_yticklabels([])
    ax_d.set_xticklabels(xt, fontsize=6.5)
    ax_d.set_xlabel("")

    from matplotlib.patches import Patch
    ax_leg.legend(
        handles=[
            Patch(facecolor=amp_c, edgecolor="gray", label="amp present"),
            Patch(facecolor=del_c, edgecolor="gray", label="del present"),
            Patch(facecolor="white", edgecolor="gray", label="absent"),
        ],
        loc="center left", fontsize=6.5, frameon=False,
    )
    theme.finalize(fig, "none")
    save(fig, cfg, "D_driver_heatmap")


def plot_gene_jaccard(jaccard_df: pd.DataFrame, cfg: dict[str, Any]) -> None:
    if jaccard_df.empty:
        return
    style(cfg)
    g = sns.catplot(data=jaccard_df, x="comparison", y="jaccard",
                    hue="ptype", col="level", kind="bar",
                    height=3.4, aspect=0.92,
                    palette=theme.PTYPE_COLORS, edgecolor="white",
                    linewidth=0.4, errorbar=None)
    comps = list(dict.fromkeys(jaccard_df["comparison"]))
    tick_lab = []
    for c in comps:
        s = jaccard_df[jaccard_df["comparison"] == c]
        if len(s) and "n_a" in s.columns:
            tick_lab.append(f"{c}\n({int(s['n_a'].iloc[0])}/"
                            f"{int(s['n_b'].iloc[0])})")
        else:
            tick_lab.append(c)
    g.set_xticklabels(tick_lab, rotation=20, ha="right", fontsize=6.5)
    g.set_axis_labels("comparison  (n = set A / set B)", "Jaccard")
    for ax in g.axes.flat:
        ax.set_ylim(0, max(0.4, ax.get_ylim()[1] * 1.1))
        for cont in ax.containers:
            ax.bar_label(cont, fmt="%.2f", fontsize=5.5, padding=1,
                         color="#444444")
    for k, ax in enumerate(g.axes.flat):
        theme.panel_label(ax, theme.panel_letter(k), cfg)
    if g.legend is not None:
        g.legend.set_title("")
    theme.finalize(g.fig, "tight")
    save(g.fig, cfg, "D_set_jaccard")


