"""Figures for module P6 (arm-level broad landscape)."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ._core import *  # noqa: F401,F403

_ARM_ORDER = [f"{n}{s}" for n in range(1, 23) for s in ("p", "q")]


def _arm_sort_key(arm: str) -> tuple[int, int]:
    a = str(arm)
    try:
        return int(a[:-1]), 0 if a[-1] == "p" else 1
    except ValueError:
        return 99, 2


def plot_p6_arm_landscape(diff: pd.DataFrame, cfg) -> None:
    """Diverging arm-level delta-frequency bars, one panel per direction."""
    if diff is None or diff.empty:
        return
    style(cfg)
    comps = list(dict.fromkeys(diff["comparison"]))
    nrows, _ = theme.grid_shape(len(comps), max_cols=1)
    fig, axes = plt.subplots(nrows, 2, figsize=theme.figsize(
        cfg, panel_h=3.1, total_w=7.0), squeeze=False)
    for row, comp in enumerate(comps):
        d = diff[diff["comparison"] == comp]
        la, lb, ca, cb = _sides(cfg, comp)
        for col, direction in enumerate(("Amp", "Del")):
            ax = axes[row][col]
            dd = d.set_index("arm").reindex(
                [a for a in _ARM_ORDER if a in set(d["arm"])])
            if dd.empty:
                ax.set_visible(False)
                continue
            delta = dd[f"delta_freq_{direction.lower()}"].to_numpy(float)
            fdr = dd[f"fdr_{direction.lower()}"].to_numpy(float)
            x = np.arange(len(dd))
            colors = [ca if v > 0 else cb for v in delta]
            ax.bar(x, delta, width=0.75, color=colors, edgecolor="white",
                   linewidth=0.3)
            ax.axhline(0, color="#555555", lw=0.6)
            # BH-significant arms get a star
            for xi, (v, f) in enumerate(zip(delta, fdr)):
                if np.isfinite(f) and f < 0.05:
                    ax.text(xi, v + (0.012 if v >= 0 else -0.012), "*",
                            ha="center",
                            va="bottom" if v >= 0 else "top",
                            fontsize=6, color="#111111")
            ax.set_xticks(x)
            ax.set_xticklabels(list(dd.index), fontsize=4.5, rotation=90)
            ax.set_ylabel(f"\u0394freq ({direction})\n[{la} \u2212 {lb}]",
                          fontsize=6.5)
            n_sig = int(np.sum(np.isfinite(fdr) & (fdr < 0.05)))
            ax.set_title(f"{direction} \u2014 {comp}  (FDR<0.05: n={n_sig})",
                         fontsize=7)
            ax.tick_params(axis="y", labelsize=6)
            ax.grid(axis="y", alpha=0.2)
            lim = max(0.02, np.nanmax(np.abs(delta)) * 1.3)
            ax.set_ylim(-lim, lim)
    theme.finalize(fig, "tight")
    save(fig, cfg, "P6_arm_landscape")


def plot_p6_aneuploidy(aneu: pd.DataFrame, tests: pd.DataFrame, cfg) -> None:
    """Broad-event burden per sample: box + jittered points per cohort.

    One figure row per comparison; the comparison id is carried in each
    panel title so nothing depends on a suptitle surviving the export.
    """
    if aneu is None or aneu.empty:
        return
    style(cfg)
    metrics = [("n_broad_amp", "broad amplifications per sample"),
               ("n_broad_del", "broad deletions per sample"),
               ("broad_frac", "fraction of arms broadly altered")]
    comps = list(dict.fromkeys(aneu["comparison"]))
    nrows, _ = theme.grid_shape(len(comps), max_cols=1)
    fig, axes = plt.subplots(nrows, len(metrics), figsize=theme.figsize(
        cfg, panel_h=2.5, total_w=7.0), squeeze=False)
    pal = _palette(cfg)
    for row, comp in enumerate(comps):
        d = aneu[aneu["comparison"] == comp]
        t = (tests[tests["comparison"] == comp] if tests is not None
             and len(tests) else pd.DataFrame())
        order = list(dict.fromkeys(d["cohort"]))
        for col, (metric, label) in enumerate(metrics):
            ax = axes[row][col]
            data = [d.loc[d["cohort"] == c, metric].dropna().to_numpy(float)
                    for c in order]
            if all(len(x) == 0 for x in data):
                ax.set_visible(False)
                continue
            # matplotlib box positions are 1-based; draw the jittered points
            # at the same positions (seaborn's categorical axis would put
            # them at 0-based positions and misalign the panels)
            bp = ax.boxplot(data, tick_labels=order, patch_artist=True,
                            widths=0.5, medianprops=dict(color="#333333"))
            for patch, c in zip(bp["boxes"], order):
                patch.set_facecolor(pal.get(c, "#888888"))
                patch.set_alpha(0.45)
            rng = np.random.default_rng(row * 10 + col)
            for xi, (c, vals) in enumerate(zip(order, data), start=1):
                if len(vals) == 0:
                    continue
                ax.scatter(rng.normal(xi, 0.07, len(vals)), vals, s=7,
                           color=pal.get(c, "#888888"), alpha=0.65,
                           edgecolor="none", zorder=4)
            if not t.empty and f"{metric}_p_mw" in t.columns:
                p = t.iloc[0][f"{metric}_p_mw"]
                theme.annot_stat(ax, theme.fmt_p(
                    float(p) if np.isfinite(p) else None), loc="upper left")
            theme.annot_n(ax, theme.group_sizes(
                d, "cohort", order=order), loc="upper right", fontsize=5.5)
            ax.set_title(f"{label} \u2014 {comp}", fontsize=7)
            ax.tick_params(axis="both", labelsize=6.5)
            ax.grid(axis="y", alpha=0.2)
    theme.finalize(fig, "tight")
    save(fig, cfg, "P6_aneuploidy")


def figures_p6(res: dict, cfg) -> None:
    plot_p6_arm_landscape(res.get("arm_diff"), cfg)
    plot_p6_aneuploidy(res.get("aneuploidy"), res.get("aneuploidy_tests"), cfg)
