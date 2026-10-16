"""Figures for module ABS (absolute copy-number comparisons)."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from ..config import out_table
from ._core import *  # noqa: F401,F403


def plot_abs_wgd(wgd: pd.DataFrame, cfg) -> None:
    """WGD+ fraction per cohort with Wilson CIs and the Fisher p."""
    if wgd is None or wgd.empty:
        return
    style(cfg)
    main = wgd[~wgd["comparison"].str.contains("high-purity", na=False)]
    comps = list(dict.fromkeys(main["comparison"]))
    nrows, ncols = theme.grid_shape(len(comps), max_cols=2)
    fig, axes = plt.subplots(nrows, ncols, figsize=theme.figsize(
        cfg, panel_h=2.6, total_w=6.2 if ncols == 2 else 3.5), squeeze=False)
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
    # stratified rows ("<comp> [WGD+ only]" ...) are drawn by
    # plot_loh_wgd_strata; the boxplots here show the pooled comparison
    metrics = metrics[~metrics["comparison"].str.contains("\\[", na=False)]
    if metrics.empty:
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
                        loc="top center", fontsize=6.0)
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
    bar_w = 0.72                       # tighter gap between the two cohorts

    ax = axes[0]
    cats = ["0x", "1x", ">=2x"]
    tab = (pd.crosstab(d["cohort"], d["wgd_class"])
           .reindex(index=cohs, columns=cats).fillna(0).astype(int))
    frac = tab.div(tab.sum(axis=1), axis=0)
    bottom = np.zeros(len(cohs))
    for c in cats:
        v = frac[c].to_numpy()
        ax.bar(range(len(cohs)), v, bottom=bottom, color=_WGD_CLASS_COLORS[c],
               width=bar_w, label=_WGD_CLASS_LABEL[c], edgecolor="white",
               linewidth=0.5)
        for i, (vv, bb) in enumerate(zip(v, bottom)):
            if vv > 0.04:
                ax.text(i, bb + vv / 2, f"{vv:.0%}\n({tab[c].iloc[i]})",
                        ha="center", va="center", fontsize=6,
                        color="white" if c != "0x" else "#333333")
        bottom += v
    try:
        chi2, p_chi, _, _ = stats.chi2_contingency(tab.to_numpy())
        chi_txt = f"$\\chi^2$={chi2:.2f}, p={p_chi:.2f}"
    except ValueError:      # a class is empty in one cohort: no test
        chi_txt = "no test (empty class)"
    excl = {c: len(t) - int(tab.loc[c].sum())
            for c, t in zip(cohs, cn_tables.values())}
    ax.set_xticks(range(len(cohs)))
    ax.set_xticklabels([f"{c}\n(n={int(tab.loc[c].sum())}, "
                        f"excl. {excl[c]})" for c in cohs], fontsize=7)
    ax.set_ylabel("Fraction of samples")
    ax.set_ylim(0, 1.0)
    ax.set_title(f"Strict (nMinor$\\geq$2 signature)  {chi_txt}",
                 fontsize=7.5)
    ax.tick_params(labelsize=6.5)
    handles, labels = ax.get_legend_handles_labels()

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
               width=bar_w, label=c, edgecolor="white", linewidth=0.5)
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
    ax.tick_params(labelsize=6.5)

    total = sum(len(t) for t in cn_tables.values())
    fig.suptitle(f"WGD status by cohort (assessable samples only, n={len(d)}; "
                 f"excluding {total - len(d)} purity=1.0 non-informative)",
                 fontsize=8)
    # single shared legend to the right of both panels: a panel-internal
    # legend sits on group B's stacked bar at this bar count
    fig.legend(handles, labels, loc="outside center right", frameon=False,
               fontsize=6.5, title="WGD class", title_fontsize=6.5)
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


_HRD_METRICS = (("HRD_sum", "HRD-sum"), ("HRD_LOH", "HRD-LOH"),
                ("TAI", "Telomeric AI"), ("LST", "LST"))


def plot_hrd_dists(hrd: pd.DataFrame, cfg) -> None:
    """scarHRD score distributions, CUE vs CHB: density (left) + ECDF
    (right) per metric.  Informative samples with purity >= 0.30 only —
    low-purity samples inflate the scores via segmentation noise
    (HRD_sum ~ n_segments rho = 0.9, ~ purity rho = -0.6 in this cohort).
    """
    d = hrd.copy()
    d["cohort"] = d["cohort"].str.lower()
    s = d[(d.get("qc_tier", "").astype(str) == "informative")
          & (d["purity"] >= 0.30)].copy()
    if s.empty:
        return
    style(cfg)
    col = {"chb": theme.OKABE["blue"], "cue": theme.OKABE["vermillion"]}
    cohort_lab = {c: c.upper() for c in col}
    for metric, mlab in _HRD_METRICS:
        if metric not in s.columns:
            continue
        a = s[s.cohort == "chb"][metric].to_numpy(float)
        b = s[s.cohort == "cue"][metric].to_numpy(float)
        fig, axes = plt.subplots(1, 2, figsize=theme.figsize(
            cfg, total_w=6.8, panel_h=3.2))
        ax = axes[0]
        grid = np.linspace(0, max(a.max(), b.max()) + 5, 300)
        for arr, c in ((a, "chb"), (b, "cue")):
            if len(arr) > 2 and arr.std() > 0:
                kde = stats.gaussian_kde(arr, bw_method=0.35)
                ax.fill_between(grid, kde(grid), alpha=0.25, color=col[c])
                ax.plot(grid, kde(grid), color=col[c], lw=1.2,
                        label=f"{cohort_lab[c]} (n={len(arr)}, "
                              f"med={np.median(arr):.0f})")
                ax.plot(arr, np.full_like(arr, -0.001), "|", color=col[c],
                        ms=5, alpha=0.5)
        if metric == "HRD_sum":
            ax.axvline(42, color="#555555", lw=0.7, ls=":")
            ax.text(42.5, ax.get_ylim()[1] * 0.02, "HRD+ 42", fontsize=5.5,
                    rotation=90, va="bottom", color="#555555")
        p_mw = stats.mannwhitneyu(a, b)[1]
        ax.set_xlabel(mlab)
        ax.set_ylabel("Density")
        ax.set_title(f"density  |  MW p={p_mw:.2f}", fontsize=7.5)
        ax.legend(frameon=False, fontsize=6)
        ax.tick_params(labelsize=6.5)
        ax = axes[1]
        for arr, c in ((a, "chb"), (b, "cue")):
            xs = np.sort(arr)
            ys = np.arange(1, len(xs) + 1) / len(xs)
            ax.step(xs, ys, where="post", color=col[c], lw=1.2,
                    label=cohort_lab[c])
        ks = stats.ks_2samp(a, b)
        ax.annotate(f"KS D={ks.statistic:.2f}\np={ks.pvalue:.2f}",
                    xy=(0.97, 0.06), xycoords="axes fraction", ha="right",
                    fontsize=6.5)
        ax.set_xlabel(mlab)
        ax.set_ylabel("Empirical CDF")
        ax.set_title("ECDF", fontsize=7.5)
        ax.legend(frameon=False, fontsize=6, loc="lower right")
        ax.tick_params(labelsize=6.5)
        fig.suptitle(f"{mlab} distribution, CUE vs CHB (informative, "
                     f"purity$\\geq$0.30, n={len(s)})", fontsize=8)
        theme.finalize(fig, "constrained", wspace=0.25)
        save(fig, cfg, f"ABS_hrd_dist_{metric}")


_LOH_METRICS = (("loh_frac", "LOH fraction"),
                ("cnloh_frac", "copy-neutral LOH"))


def plot_loh_wgd_strata(metrics: pd.DataFrame, scores: pd.DataFrame, cfg
                        ) -> None:
    """LOH by cohort *within* WGD strata (one comparison per figure).

    The pooled test mixes a WGD+ stratum (minor CN >= 1 is common, so LOH
    segments are rarer) with a WGD- stratum; its result therefore follows
    the cohorts' WGD composition and can even reverse sign (true in the
    CUE/CHB data).  Boxes are per cohort within each stratum, with the
    stratified Mann-Whitney p; the pooled p is in the panel title for
    contrast.  The figure pins itself to the first stratified comparison so
    the annotations can never mix comparisons.
    """
    if metrics is None or metrics.empty or scores is None or scores.empty:
        return
    if "wgd_status" not in scores.columns:
        return
    strat = metrics[metrics["comparison"].str.contains("\\[", na=False)]
    if strat.empty:
        return
    comp = str(strat["comparison"].iloc[0]).split(" [")[0]
    ca, cb = comp.split("_vs_")[:2]
    cohs = [ca, cb]
    strat = strat[strat["comparison"].str.startswith(comp)]
    if "cohort" in scores.columns:
        scores = scores[scores["cohort"].isin(cohs)]
    style(cfg)
    pal = _palette(cfg)
    fig, axes = plt.subplots(1, 2, figsize=theme.figsize(
        cfg, total_w=6.8, panel_h=3.0))
    for ax, (metric, label) in zip(axes, _LOH_METRICS):
        pooled = metrics[(metrics["metric"] == metric)
                         & (metrics["comparison"] == comp)]
        xs_all: list[float] = []
        labs_all: list[str] = []
        for gi, (sval, slab) in enumerate(((0, "WGD-"), (1, "WGD+"))):
            data = []
            for ci, c in enumerate(cohs):
                v = scores.loc[(scores["cohort"] == c)
                               & (scores["wgd_status"] == sval),
                               metric].dropna().to_numpy(float)
                data.append(v)
            pos = [gi * 3.0 + ci * 0.9 for ci in range(len(cohs))]
            if all(len(v) == 0 for v in data):
                continue
            bp = ax.boxplot(data, positions=pos, widths=0.62,
                            patch_artist=True,
                            medianprops=dict(color="#333333"))
            for patch, c in zip(bp["boxes"], cohs):
                patch.set_facecolor(pal.get(c, "#888888"))
                patch.set_alpha(0.45)
            rng = np.random.default_rng(gi * 7 + int(sval))
            for xi, v, c in zip(pos, data, cohs):
                if len(v) == 0:
                    continue
                ax.scatter(rng.normal(xi, 0.07, len(v)), v, s=5,
                           color=pal.get(c, "#888888"), alpha=0.55,
                           edgecolor="none", zorder=4)
            xs_all += pos
            labs_all += [c.upper() for c in cohs]
            r = strat[(strat["metric"] == metric)
                      & (strat["comparison"].str.contains(f"{slab} only",
                                                          regex=False))]
            txt = slab
            if len(r):
                rr = r.iloc[0]
                txt = (f"{slab}:  {cohs[0]} {rr['median_a']:.3f} vs "
                       f"{cohs[1]} {rr['median_b']:.3f}\n"
                       f"MW p={rr['p_mw']:.2f} (n={rr['n_a']}/{rr['n_b']})")
            ax.annotate(txt, xy=(float(np.mean(pos)), 1.005),
                        xycoords=("data", "axes fraction"), ha="center",
                        va="bottom", fontsize=5.6, color="#333333",
                        linespacing=1.4)
        ax.axvline(1.35, color="#BBBBBB", lw=0.6, ls=":")
        ptxt = ""
        if len(pooled):
            ptxt = (f"  \u2014  pooled: p={pooled['p_mw'].iloc[0]:.2f} "
                    f"({pooled['median_a'].iloc[0]:.3f} vs "
                    f"{pooled['median_b'].iloc[0]:.3f})")
        ax.set_xticks(xs_all)
        ax.set_xticklabels(labs_all, fontsize=6.5)
        ax.set_ylabel(label, fontsize=7)
        ax.set_title(f"{label}, by WGD stratum{ptxt}", fontsize=7, pad=22)
        ax.tick_params(labelsize=6)
        ax.grid(axis="y", alpha=0.2)
    fig.suptitle(f"Module ABS \u2014 {comp}: LOH by cohort within WGD strata "
                 "(the pooled comparison is confounded by WGD composition)",
                 fontsize=8)
    theme.finalize(fig, "tight")
    save(fig, cfg, "ABS_loh_wgd_strata")


_LOH_CLASS_COLORS = {
    "del+LOH": theme.OKABE["vermillion"],
    "del-only": theme.OKABE["sky_blue"],
    "LOH-only": theme.OKABE["orange"],
    "neither": "#DDDDDD",
}


def plot_loh_del_cross(cross: pd.DataFrame, cfg) -> None:
    """LOH vs GISTIC's del call at the same peaks.

    Left: per peak, GISTIC's del frequency (x) against the fraction of the
    cohort with LOH at that peak (y) — points above the diagonal are LOH
    the del track under-reports, sized by peak width.  Right: composition
    of the four classes (del+LOH / del-only / LOH-only / neither) for each
    cohort's most LOH-only peaks plus the cohort total.  One comparison per
    figure (the first in the table).
    """
    if cross is None or cross.empty:
        return
    comp = str(cross["comparison"].iloc[0])
    cross = cross[cross["comparison"] == comp]
    if cross.empty:
        return
    style(cfg)
    pal = _palette(cfg)
    present = list(dict.fromkeys(cross["cohort"]))
    cohs = [c for c in comp.split("_vs_") if c in present] or present
    fig, axes = plt.subplots(1, 2, figsize=theme.figsize(
        cfg, total_w=7.0, panel_h=3.4),
        gridspec_kw=dict(width_ratios=[1.0, 1.15]))

    ax = axes[0]
    lim = float(np.nanmax([cross["gistic_freq"].max(),
                           cross["loh_frac_any"].max()])) * 1.12
    lim = max(lim, 0.1)
    ax.plot([0, lim], [0, lim], color="#999999", lw=0.7, ls="--")
    for c in cohs:
        d = cross[cross["cohort"] == c]
        ax.scatter(d["gistic_freq"], d["loh_frac_any"],
                   s=np.clip(np.log10(d["kb"].to_numpy(float) + 10) * 9, 6, 90),
                   color=pal.get(c, "#888888"), alpha=0.75,
                   edgecolor="none", label=c.upper())
    top = cross.nlargest(4, "loh_only_frac")
    for _, r in top.iterrows():
        ax.annotate(f"{r['cohort'].upper()} {r['cytoband']}",
                    xy=(r["gistic_freq"], r["loh_frac_any"]),
                    xytext=(3, 3), textcoords="offset points", fontsize=5.2,
                    color="#444444")
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel("GISTIC del frequency (amplitude call)", fontsize=7)
    ax.set_ylabel("LOH frequency at the same peak", fontsize=7)
    ax.set_title("above the line = LOH the del call misses\n"
                 "(marker size $\\propto$ peak width)", fontsize=6.8)
    ax.legend(frameon=False, fontsize=6, loc="upper left")
    ax.tick_params(labelsize=6)
    ax.grid(alpha=0.2)

    ax = axes[1]
    rows = []
    for c in cohs:
        d = cross[cross["cohort"] == c]
        n = float(d["n_samples"].sum())
        rows.append(dict(label=f"{c.upper()}  all Del peaks",
                         frac=[d["del_loh"].sum() / n,
                               d["del_only"].sum() / n,
                               d["loh_only"].sum() / n,
                               d["neither"].sum() / n]))
        for _, r in d.nlargest(6, "loh_only_frac").iterrows():
            m = float(r["n_samples"])
            rows.append(dict(
                label=f"  {r['cytoband']} ({r['kb']:.0f} kb)",
                frac=[r["del_loh"] / m, r["del_only"] / m,
                      r["loh_only"] / m, r["neither"] / m]))
    order = ["del+LOH", "del-only", "LOH-only", "neither"]
    ypos = np.arange(len(rows))[::-1]
    left = np.zeros(len(rows))
    for k, cls in enumerate(order):
        v = np.array([r["frac"][k] for r in rows])
        ax.barh(ypos, v, left=left, height=0.72,
                color=_LOH_CLASS_COLORS[cls], edgecolor="white",
                linewidth=0.4, label=cls)
        for y, vv, ll in zip(ypos, v, left):
            if vv >= 0.09:
                ax.text(ll + vv / 2, y, f"{vv:.0%}", ha="center",
                        va="center", fontsize=5,
                        color="white" if cls in ("del+LOH", "LOH-only")
                        else "#333333")
        left += v
    ax.set_yticks(ypos)
    ax.set_yticklabels([r["label"] for r in rows], fontsize=5.6)
    ax.set_xlim(0, 1.0)
    ax.set_xlabel("fraction of the cohort's samples", fontsize=7)
    ax.set_title("composition of the four classes per peak\n"
                 "(samples $\\times$ Del peaks, q$\\leq$0.25)",
                 fontsize=6.8)
    handles, labels = ax.get_legend_handles_labels()
    ax.tick_params(labelsize=6)
    ax.grid(axis="x", alpha=0.2)
    fig.legend(handles, labels, loc="outside center right", frameon=False,
               fontsize=6.5, title="class", title_fontsize=6.5)
    fig.suptitle(f"Module ABS \u2014 {comp}: LOH vs GISTIC deletion at the "
                 "same loci", fontsize=8)
    theme.finalize(fig, "constrained", wspace=0.35)
    save(fig, cfg, "ABS_loh_del_cross")


def figures_abs(res: dict, cfg) -> None:
    plot_abs_wgd(res.get("wgd"), cfg)
    plot_abs_metric_dists(res.get("metrics"), res.get("scores"), cfg)
    plot_loh_wgd_strata(res.get("metrics"), res.get("scores"), cfg)
    plot_abs_purity_adjusted(res.get("purity_peaks"), cfg)
    plot_loh_del_cross(res.get("loh_del"), cfg)
    plot_wgd_class_stacked(res.get("cn_tables"), cfg)
    plot_ploidy_rank(res.get("cn_tables"), cfg)
    hrd = None
    try:
        hrd = pd.read_csv(out_table(cfg, "ABS_hrd_scores.tsv"), sep="\t")
    except (OSError, ValueError):
        hrd = None
    if hrd is not None and len(hrd):
        plot_hrd_dists(hrd, cfg)
