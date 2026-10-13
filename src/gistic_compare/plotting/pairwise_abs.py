"""Figures for module ABS (absolute copy-number comparisons)."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from ._core import *  # noqa: F401,F403


def plot_abs_wgd(wgd: pd.DataFrame, cfg) -> None:
    """WGD+ fraction per cohort with Wilson CIs and the Fisher p."""
    if wgd is None or wgd.empty:
        return
    style(cfg)
    main = wgd[~wgd["comparison"].str.contains("high-purity", na=False)]
    comps = list(dict.fromkeys(main["comparison"]))
    nrows, _ = theme.grid_shape(len(comps), max_cols=2)
    fig, axes = plt.subplots(nrows, 2, figsize=theme.figsize(
        cfg, panel_h=2.6, total_w=6.2), squeeze=False)
    axes = axes.ravel()
    for k, comp in enumerate(comps):
        ax = axes[k]
        r = main[main["comparison"] == comp].iloc[0]
        la, lb, _, _ = _sides(cfg, comp)
        heights = [r["wgd_frac_a"], r["wgd_frac_b"]]
        lo, hi = theme.wilson_ci(
            [r["wgd_n_a"], r["wgd_n_b"]], [r["n_a"], r["n_b"]])
        theme.bar_with_ci(ax, [0, 1], heights, lo, hi,
                          color=[theme.palette(cfg).get(c, "#888888")
                                 for c in (comp.split("_vs_")[0],
                                           comp.split("_vs_")[1])])
        ax.set_xticks([0, 1])
        ax.set_xticklabels([theme.cohort_label(cfg, comp.split("_vs_")[0]),
                            theme.cohort_label(cfg, comp.split("_vs_")[1])],
                           fontsize=7)
        theme.annot_stat(ax, f"Fisher p = {r['p_fisher']:.3g}, "
                             f"OR = {r['odds_ratio']}", loc="upper left")
        ax.set_ylabel("WGD+ fraction", fontsize=7)
        ax.set_ylim(0, max(0.1, float(np.nanmax(hi)) * 1.25))
        ax.set_title(comp, fontsize=7)
        ax.tick_params(labelsize=6.5)
        ax.grid(axis="y", alpha=0.2)
    for k in range(len(comps), len(axes)):
        axes[k].set_visible(False)
    theme.finalize(fig, "tight")
    save(fig, cfg, "ABS_wgd_freq")


def plot_abs_metric_dists(metrics: pd.DataFrame, scores: pd.DataFrame, cfg
                          ) -> None:
    """purity / ploidy / LOH / CN-LOH distributions per cohort."""
    if metrics is None or metrics.empty or scores is None or scores.empty:
        return
    style(cfg)
    pal = _palette(cfg)
    comps = list(dict.fromkeys(metrics["comparison"]))
    nrows, _ = theme.grid_shape(len(comps), max_cols=1)
    fig, axes = plt.subplots(nrows, 4, figsize=theme.figsize(
        cfg, panel_h=2.4, total_w=7.0), squeeze=False)
    panels = [("purity", "purity"), ("ploidy_abs", "absolute ploidy"),
              ("loh_frac", "LOH fraction"),
              ("cnloh_frac", "copy-neutral LOH")]
    for row, comp in enumerate(comps):
        d = scores[scores["cohort"].isin(comp.split("_vs_"))]
        t = metrics[metrics["comparison"] == comp]
        order = list(dict.fromkeys(d["cohort"]))
        for col, (metric, label) in enumerate(panels):
            ax = axes[row][col]
            data = [d.loc[d["cohort"] == c, metric].dropna().to_numpy(float)
                    for c in order]
            if all(len(x) == 0 for x in data):
                ax.set_visible(False)
                continue
            bp = ax.boxplot(data, tick_labels=order, patch_artist=True,
                            widths=0.5, medianprops=dict(color="#333333"))
            for patch, c in zip(bp["boxes"], order):
                patch.set_facecolor(pal.get(c, "#888888"))
                patch.set_alpha(0.45)
            rng = np.random.default_rng(row * 4 + col)
            for xi, vals in enumerate(data, start=1):
                if len(vals) == 0:
                    continue
                ax.scatter(rng.normal(xi, 0.07, len(vals)), vals, s=5,
                           color=pal.get(order[xi - 1], "#888888"),
                           alpha=0.55, edgecolor="none", zorder=4)
            if len(t) and "p_mw" in t.columns:
                r = t[t["metric"] == metric]
                if len(r):
                    pv = r["p_mw"].iloc[0]
                    theme.annot_stat(
                        ax, theme.fmt_p(float(pv)
                                        if np.isfinite(pv) else None),
                        loc="upper left", fontsize=6.0)
            ax.set_title(label, fontsize=6.8)
            ax.tick_params(labelsize=6)
            ax.grid(axis="y", alpha=0.2)
    fig.suptitle("Module ABS \u2014 absolute-CN metrics per sample",
                 fontsize=8)
    theme.finalize(fig, "tight")
    save(fig, cfg, "ABS_metric_dists")


def plot_abs_purity_adjusted(peaks: pd.DataFrame, cfg) -> None:
    """Crude vs purity-adjusted OR per peak: does the cohort effect survive
    purity adjustment?"""
    if peaks is None or peaks.empty:
        return
    d = peaks[np.isfinite(peaks["crude_or"])
              & np.isfinite(peaks["adj_or_purity"])]
    if d.empty:
        return
    style(cfg)
    comp = d["comparison"].iloc[0]
    fig, ax = plt.subplots(figsize=theme.figsize(cfg, panel_h=2.9,
                                                 total_w=3.5))
    sig = d["adj_q_purity"] < 0.05 if "adj_q_purity" in d.columns \
        else pd.Series(False, index=d.index)
    ax.scatter(d.loc[~sig, "crude_or"], d.loc[~sig, "adj_or_purity"],
               s=10, color="#8FA8C8", alpha=0.8, edgecolor="none")
    ax.scatter(d.loc[sig, "crude_or"], d.loc[sig, "adj_or_purity"], s=14,
               color=theme.OKABE["vermillion"], edgecolor="none",
               label="adj q < 0.05")
    lim = max(0.5, float(np.nanmax([d["crude_or"].max(),
                                    d["adj_or_purity"].max()])) * 1.05)
    ax.plot([0, lim], [0, lim], color="#999999", lw=0.7, ls="--")
    ax.axvline(1, color="#777777", lw=0.6, ls=":")
    ax.axhline(1, color="#777777", lw=0.6, ls=":")
    ax.set_xlabel("crude OR (cohort)", fontsize=7)
    ax.set_ylabel("purity-adjusted OR", fontsize=7)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_title(f"{comp}: does the effect survive purity adjustment?",
                 fontsize=7)
    ax.legend(fontsize=5.5, loc="upper left")
    ax.grid(alpha=0.2)
    theme.finalize(fig, "tight")
    save(fig, cfg, "ABS_purity_adjusted")


_WGD_CLASS_COLORS = {"0x": "#BBBBBB", "1x": theme.OKABE["blue"],
                     ">=2x": theme.OKABE["vermillion"]}
_WGD_CLASS_LABEL = {"0x": "no WGD", "1x": "1x WGD", ">=2x": "$\\geq$2x WGD"}


def _wgd_informative(cn_tables: dict) -> pd.DataFrame:
    parts = [t[t.get("qc_tier", "").astype(str) == "informative"]
             for t in cn_tables.values()]
    parts = [p.assign(cohort=cid) for cid, p in zip(cn_tables, parts)
             if len(p)]
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def plot_wgd_class_stacked(cn_tables: dict, cfg) -> None:
    """WGD-status composition per cohort, strict 3-class vs loose binary.

    Strict: nMinor >= 2 signature (parsimonious whole-genome-doubling
    evidence; (2,1)-type segments are ambiguous and excluded).  Loose:
    the legacy binary wgd_status (includes ambiguous (2,1) segments —
    sensitive but not specific).  Non-informative samples (purity >= 0.99)
    are excluded and their counts annotated.
    """
    d = _wgd_informative(cn_tables)
    if d.empty:
        return
    d = d.copy()
    d["wgd_status"] = pd.to_numeric(d["wgd_status"], errors="coerce").map(
        {0: "WGD-", 1: "WGD+"})
    style(cfg)
    cohs = list(cn_tables.keys())
    fig, axes = plt.subplots(1, 2, figsize=theme.figsize(
        cfg, total_w=6.4, panel_h=3.3))

    ax = axes[0]
    cats = ["0x", "1x", ">=2x"]
    tab = (pd.crosstab(d["cohort"], d["wgd_class"])
           .reindex(index=cohs, columns=cats).fillna(0).astype(int))
    frac = tab.div(tab.sum(axis=1), axis=0)
    bottom = np.zeros(len(cohs))
    for c in cats:
        v = frac[c].to_numpy()
        ax.bar(range(len(cohs)), v, bottom=bottom, color=_WGD_CLASS_COLORS[c],
               width=0.62, label=_WGD_CLASS_LABEL[c], edgecolor="white",
               linewidth=0.5)
        for i, (vv, bb) in enumerate(zip(v, bottom)):
            if vv > 0.04:
                ax.text(i, bb + vv / 2, f"{vv:.0%}\n({tab[c].iloc[i]})",
                        ha="center", va="center", fontsize=6,
                        color="white" if c != "0x" else "#333333")
        bottom += v
    chi2, p_chi, _, _ = stats.chi2_contingency(tab.to_numpy())
    excl = {c: len(t) - int(tab.loc[c].sum())
            for c, t in zip(cohs, cn_tables.values())}
    ax.set_xticks(range(len(cohs)))
    ax.set_xticklabels([f"{c}\n(n={int(tab.loc[c].sum())}, "
                        f"excl. {excl[c]})" for c in cohs], fontsize=7)
    ax.set_ylabel("Fraction of samples")
    ax.set_ylim(0, 1.0)
    ax.set_title(f"Strict (nMinor$\\geq$2 signature)  "
                 f"$\\chi^2$={chi2:.2f}, p={p_chi:.2f}", fontsize=7.5)
    ax.legend(frameon=False, fontsize=6, loc="upper right")
    ax.tick_params(labelsize=6.5)

    ax = axes[1]
    cats2 = ["WGD-", "WGD+"]
    tab2 = (pd.crosstab(d["cohort"], d["wgd_status"])
            .reindex(index=cohs, columns=cats2).fillna(0).astype(int))
    frac2 = tab2.div(tab2.sum(axis=1), axis=0)
    bottom = np.zeros(len(cohs))
    cols2 = {"WGD-": "#BBBBBB", "WGD+": theme.OKABE["vermillion"]}
    for c in cats2:
        v = frac2[c].to_numpy()
        ax.bar(range(len(cohs)), v, bottom=bottom, color=cols2[c],
               width=0.62, label=c, edgecolor="white", linewidth=0.5)
        for i, (vv, bb) in enumerate(zip(v, bottom)):
            if vv > 0.04:
                ax.text(i, bb + vv / 2, f"{vv:.0%}\n({tab2[c].iloc[i]})",
                        ha="center", va="center", fontsize=6,
                        color="white" if c != "WGD-" else "#333333")
        bottom += v
    odds, p_f = stats.fisher_exact(tab2.to_numpy())
    ax.set_xticks(range(len(cohs)))
    ax.set_xticklabels([f"{c}\n(n={int(tab2.loc[c].sum())})" for c in cohs],
                       fontsize=7)
    ax.set_ylabel("Fraction of samples")
    ax.set_ylim(0, 1.0)
    ax.set_title(f"Loose (incl. ambiguous (2,1) segs)  Fisher p={p_f:.2f}",
                 fontsize=7.5)
    ax.legend(frameon=False, fontsize=6, loc="upper right")
    ax.tick_params(labelsize=6.5)

    total = sum(len(t) for t in cn_tables.values())
    fig.suptitle(f"WGD status by cohort (assessable samples only, n={len(d)}; "
                 f"excluding {total - len(d)} purity=1.0 non-informative)",
                 fontsize=8)
    theme.finalize(fig, "constrained", wspace=0.25)
    save(fig, cfg, "WGD_class_stacked")


def plot_ploidy_rank(cn_tables: dict, cfg) -> None:
    """Ploidy rank plot per cohort (ASCAT-release style): samples sorted by
    ploidy, cohort median in red, points coloured by WGD class."""
    d = _wgd_informative(cn_tables)
    if d.empty:
        return
    style(cfg)
    cohs = list(cn_tables.keys())
    fig, ax = plt.subplots(figsize=theme.figsize(cfg, total_w=4.0,
                                                 panel_h=3.4))
    for i, c in enumerate(cohs):
        ax.axvspan(i - 0.40, i + 0.40, color="#D9D9D9" if i % 2 else "white",
                   zorder=0)
    for i, c in enumerate(cohs):
        g = d[d["cohort"] == c].sort_values("ploidy_abs").reset_index(
            drop=True)
        g["xr"] = i - 0.36 + 0.72 * (np.arange(1, len(g) + 1) / len(g))
        for cls in ("0x", "1x", ">=2x"):
            s = g[g["wgd_class"] == cls]
            if s.empty:
                continue
            ax.scatter(s["xr"], s["ploidy_abs"], s=9,
                       color=_WGD_CLASS_COLORS[cls], edgecolors="none",
                       zorder=3, label=_WGD_CLASS_LABEL[cls] if i == 0 else None)
        med = float(g["ploidy_abs"].median())
        ax.plot([i - 0.40, i + 0.40], [med, med], color="red", lw=1.3,
                zorder=4)
        ax.text(i + 0.28, med + 0.10, f"{med:.2f}", ha="center", fontsize=6.5,
                color="red", zorder=5)
    ax.axhline(2.0, color="#555555", lw=0.5, ls=":", zorder=1)
    ax.set_xticks(range(len(cohs)))
    ax.set_xticklabels([c.upper() for c in cohs], fontsize=8)
    ax.set_xlim(-0.55, len(cohs) - 0.45)
    ax.set_ylabel("Ploidy", fontsize=8)
    ax.tick_params(labelsize=6.5)
    ax.legend(frameon=False, fontsize=6.5, loc="upper left",
              title="WGD class", title_fontsize=6.5)
    theme.finalize(fig, "constrained")
    save(fig, cfg, "ABS_ploidy_rank")


def figures_abs(res: dict, cfg) -> None:
    plot_abs_wgd(res.get("wgd"), cfg)
    plot_abs_metric_dists(res.get("metrics"), res.get("scores"), cfg)
    plot_abs_purity_adjusted(res.get("purity_peaks"), cfg)
    plot_wgd_class_stacked(res.get("cn_tables"), cfg)
    plot_ploidy_rank(res.get("cn_tables"), cfg)
