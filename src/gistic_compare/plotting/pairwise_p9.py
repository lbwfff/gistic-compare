"""Figures for module P9 (three-level SCNA burden + wGII)."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from ..pairwise_burden import METRICS
from ._core import *  # noqa: F401,F403


def plot_p9_dists(scores: pd.DataFrame, tests: pd.DataFrame, cfg) -> None:
    """Per-metric burden distributions, one figure per comparison:
    box + jittered points, MW p annotated, cohort colours from the palette."""
    if scores is None or scores.empty:
        return
    style(cfg)
    pal = _palette(cfg)
    if "comparison" not in tests.columns or tests.empty:
        comps = [(" vs ".join(sorted(set(scores["cohort"]))), scores)]
    elif "comparison" in scores.columns:
        comps = [(c, scores[scores["comparison"] == c])
                 for c in dict.fromkeys(tests["comparison"])]
    else:
        comps = [(c, scores[scores["cohort"].isin(c.split("_vs_"))])
                 for c in dict.fromkeys(tests["comparison"])]
    for comp, d in comps:
        t = tests[tests["comparison"] == comp] if len(tests) \
            else pd.DataFrame()
        order = [c for c in dict.fromkeys(d["cohort"])]
        fig, axes = plt.subplots(1, 5, figsize=theme.figsize(
            cfg, panel_h=2.5, total_w=7.0))
        for mi, (ax, (metric, label)) in enumerate(zip(axes, METRICS)):
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
            rng = np.random.default_rng(mi)
            for xi, (c, vals) in enumerate(zip(order, data), start=1):
                if len(vals) == 0:
                    continue
                ax.scatter(rng.normal(xi, 0.07, len(vals)), vals, s=6,
                           color=pal.get(c, "#888888"), alpha=0.6,
                           edgecolor="none", zorder=4)
            if len(t):
                r = t[t["metric"] == metric]
                if len(r) and "p_mw" in r.columns:
                    pv = r["p_mw"].iloc[0]
                    theme.annot_stat(
                        ax, theme.fmt_p(float(pv)
                                        if np.isfinite(pv) else None),
                        loc="top center", fontsize=6.0)
            ax.set_title(label, fontsize=6.8)
            ax.tick_params(labelsize=6)
            ax.grid(axis="y", alpha=0.2)
        fig.suptitle(f"{comp}: per-sample burden", fontsize=8)
        theme.finalize(fig, "tight")
        save(fig, cfg, "P9_burden_dists"
             if len(comps) == 1 else f"P9_burden_dists_{comp}")


def plot_p9_forest(tests: pd.DataFrame, cfg) -> None:
    """Cohort coefficient (beta) with 95% CI per metric, per comparison."""
    if tests is None or tests.empty:
        return
    style(cfg)
    comps = list(dict.fromkeys(tests["comparison"]))
    nrows, _ = theme.grid_shape(len(comps), max_cols=1)
    fig, axes = plt.subplots(nrows, 1, figsize=theme.figsize(
        cfg, panel_h=1.9, total_w=5.2), squeeze=False)
    pal = _palette(cfg)
    for row, comp in enumerate(comps):
        ax = axes[row][0]
        d = tests[tests["comparison"] == comp].copy()
        d = d[d["metric"].isin([m for m, _ in METRICS])]
        d = d.iloc[::-1]
        y = np.arange(len(d))
        la, lb, _, _ = _sides(cfg, comp)
        for yi, (_, r) in zip(y, d.iterrows()):
            color = pal.get(r["cohort_higher"], "#555555")
            ax.plot([r["ci_lo"], r["ci_hi"]], [yi, yi], color=color, lw=1.4,
                    solid_capstyle="round", zorder=3)
            ax.scatter([r["beta"]], [yi], s=22, color=color, zorder=4)
        ax.axvline(0, color="#555555", lw=0.7, ls=":")
        ax.set_yticks(y)
        ax.set_yticklabels([str(x) for x in d["label"]], fontsize=6.5)
        ax.set_xlabel(f"cohort coefficient \u03b2  "
                      f"(per unit of being {lb})", fontsize=7)
        ax.set_title(f"{comp}", fontsize=7)
        ax.tick_params(labelsize=6)
        ax.grid(axis="x", alpha=0.2)
    theme.finalize(fig, "tight")
    save(fig, cfg, "P9_effect_forest")


def plot_p9_decomposition(scores: pd.DataFrame, tests: pd.DataFrame, cfg
                          ) -> None:
    """Broad-scale vs focal-scale burden scatter (pooled percentiles)."""
    if scores is None or scores.empty:
        return
    style(cfg)
    if tests is not None and len(tests) and "comparison" in tests.columns:
        comp = list(dict.fromkeys(tests["comparison"]))[0]
    else:
        comp = " vs ".join(sorted(set(scores["cohort"])))
    ca, cb = theme.comp_sides(cfg, comp)
    pal = _palette(cfg)
    fig, ax = plt.subplots(figsize=theme.figsize(cfg, panel_h=2.9,
                                                 total_w=3.6))
    for cid in (ca, cb):
        d = scores[scores["cohort"] == cid]
        if "comparison" in scores.columns:
            d = d[d["comparison"] == comp]
        ax.scatter(d["S_arm_star"] + d["S_chrom_star"], d["S_focal_star"],
                   s=12, alpha=0.7, color=pal.get(cid, "#888888"),
                   edgecolor="none",
                   label=theme.cohort_label(cfg, cid))
    ax.set_xlabel("broad-scale burden (S arm* + S chrom*)", fontsize=7)
    ax.set_ylabel("focal burden (S focal*)", fontsize=7)
    ax.set_title(f"{comp}: which scale separates the cohorts?", fontsize=7)
    ax.legend(fontsize=6, handletextpad=0.2)
    ax.grid(alpha=0.2)
    theme.finalize(fig, "tight")
    save(fig, cfg, "P9_level_decomposition")


def plot_p9_sensitivity(sens: pd.DataFrame, verdict: pd.DataFrame, cfg
                        ) -> None:
    """Beta of S_total across the parameter grid, one panel per parameter.

    A verdict printed per panel: whether the sign of beta is consistent in
    every combination (the conclusion is parameterisation-robust)."""
    if sens is None or sens.empty:
        return
    style(cfg)
    st = sens[sens["metric"] == "S_total"]
    if st.empty:
        return
    params = ["arm_min_frac", "arm_match_diff", "focal_thr",
              "focal_broad_overlap_max"]
    fig, axes = plt.subplots(1, 4, figsize=theme.figsize(
        cfg, panel_h=2.5, total_w=7.0))
    for pi, (ax, param) in enumerate(zip(axes, params)):
        levels = sorted(st[param].unique())
        rng = np.random.default_rng(pi + 7)
        for li, lev in enumerate(levels, start=1):
            vals = st.loc[st[param] == lev, "beta"].to_numpy(float)
            ax.scatter(rng.normal(li, 0.06, len(vals)), vals, s=5,
                       color="#4C72B0", alpha=0.4, edgecolor="none")
            ax.scatter([li], [np.median(vals)], s=26, color="#B22222",
                       zorder=4)
        ax.axhline(0, color="#555555", lw=0.6, ls=":")
        ax.set_xticks(range(1, len(levels) + 1))
        ax.set_xticklabels([str(v) for v in levels], fontsize=6)
        ax.set_xlabel(param, fontsize=6.5)
        if ax is axes[0]:
            ax.set_ylabel("\u03b2 (S total)", fontsize=7)
        ax.tick_params(labelsize=6)
        ax.grid(axis="y", alpha=0.2)
    ok = (verdict is not None and len(verdict)
          and verdict.loc[verdict["metric"] == "S_total",
                          "sign_consistent"].iloc[0])
    fig.suptitle("P9 sensitivity \u2014 \u03b2 sign consistent in all "
                 f"combos: {bool(ok)}", fontsize=8)
    theme.finalize(fig, "tight")
    save(fig, cfg, "P9_sensitivity")


def _nice_step(span: float, target: int = 4) -> float:
    """A "nice" axis step (1/2/2.5/5/10 x 10**k) covering ``span`` in ~target."""
    import math

    span = abs(float(span))
    if not np.isfinite(span) or span <= 0:
        return 0.1
    mag = 10.0 ** math.floor(math.log10(span / target))
    for m in (1.0, 2.0, 5.0, 10.0):
        if m * mag >= span / target:
            return m * mag
    return 10.0 * mag


def plot_p9_events(freq: pd.DataFrame, cfg, level: str) -> None:
    """Cohort-level arm- or chromosome-event frequency, mirrored around zero.

    For each unit (arm / chromosome) an up (Amp) and a down (Del) pair of
    bars is drawn; the two cohorts are the two members of each pair and are
    told apart by colour lightness (dark = cohort A, light = cohort B) while
    the hue keeps the global Amp = vermillion / Del = blue convention.  A
    green marker at the baseline flags a cohort difference (two-sided Fisher
    across the unit, BH-FDR below the configured threshold; upper = Amp,
    lower = Del).  One figure per comparison.

    ``level`` is ``"arm"`` or ``"chrom"``.
    """
    if freq is None or freq.empty or "level" not in freq.columns:
        return
    d = freq[freq["level"] == level]
    if d.empty:
        return
    style(cfg)
    amp_dark = theme.OKABE["vermillion"]
    del_dark = theme.OKABE["blue"]
    amp_light = _blend(amp_dark, 0.55)
    del_light = _blend(del_dark, 0.55)
    sig_color = theme.OKABE["bluish_green"]
    rcfg = ((cfg.get("analysis") or {}).get("regional") or {})
    fdr_thr = float(rcfg.get("fdr", 0.05))

    comps = list(dict.fromkeys(d["comparison"]))
    for comp in comps:
        sub = d[d["comparison"] == comp]
        units = [str(u) for u in dict.fromkeys(sub["unit"])]
        n = len(units)
        if n == 0:
            continue
        uchrom = (sub.drop_duplicates("unit")
                     .set_index("unit")["chrom"].to_dict())
        val: dict[tuple[str, str], tuple[float, float, bool]] = {}
        for r in sub.itertuples():
            fa = float(r.freq_a) if np.isfinite(r.freq_a) else 0.0
            fb = float(r.freq_b) if np.isfinite(r.freq_b) else 0.0
            val[(str(r.unit), str(r.direction))] = (fa, fb, bool(r.sig))
        la, lb, _, _ = _sides(cfg, comp)

        fig, ax = plt.subplots(figsize=theme.figsize(cfg, panel_h=3.0,
                                                     total_w=7.0))
        x = np.arange(n)
        w, off = 0.34, 0.18
        maxf = 0.0
        for i, u in enumerate(units):
            fa, fb, _ = val.get((u, "Amp"), (0.0, 0.0, False))
            da, db, _ = val.get((u, "Del"), (0.0, 0.0, False))
            maxf = max(maxf, fa, fb, da, db)
            ax.bar(i - off, fa, width=w, color=amp_dark, edgecolor="white",
                   linewidth=0.3, zorder=3)
            ax.bar(i + off, fb, width=w, color=amp_light, edgecolor="white",
                   linewidth=0.3, zorder=3)
            ax.bar(i - off, -da, width=w, color=del_dark, edgecolor="white",
                   linewidth=0.3, zorder=3)
            ax.bar(i + off, -db, width=w, color=del_light, edgecolor="white",
                   linewidth=0.3, zorder=3)

        # alternating chromosome shading
        groups: list[tuple[int, int]] = []
        start, prev = 0, uchrom.get(units[0])
        for i, u in enumerate(units):
            c = uchrom.get(u)
            if c != prev:
                groups.append((start, i - 1))
                start, prev = i, c
        groups.append((start, n - 1))
        for gi, (s, e) in enumerate(groups):
            if gi % 2 == 0:
                ax.axvspan(s - 0.5, e + 0.5, color="#F0F0F0", zorder=0)

        ax.axhline(0, color="#555555", lw=0.7, zorder=4)
        L = max(0.1, maxf * 1.12)
        ax.set_ylim(-L, L)
        step = _nice_step(L)
        tick_vals = np.arange(0.0, L + step / 2, step)
        ax.set_yticks(list(tick_vals) + list(-tick_vals))

        # cohort-difference markers (baseline: upper = Amp, lower = Del)
        eps = 0.03 * L
        for i, u in enumerate(units):
            if val.get((u, "Amp"), (0, 0, False))[2]:
                ax.plot(i, eps, marker="o", ms=3.2, color=sig_color,
                        markeredgecolor="white", markeredgewidth=0.4,
                        zorder=6)
            if val.get((u, "Del"), (0, 0, False))[2]:
                ax.plot(i, -eps, marker="o", ms=3.2, color=sig_color,
                        markeredgecolor="white", markeredgewidth=0.4,
                        zorder=6)

        ax.set_xticks(x)
        labels = [u.replace("chr", "") if level == "chrom" else u
                  for u in units]
        ax.set_xticklabels(labels, fontsize=4.6)
        ax.set_xlim(-0.6, n - 0.4)
        ax.set_ylabel("frequency", fontsize=7)
        ax.set_xlabel("chromosome" if level == "chrom" else "chromosome arm",
                      fontsize=7)
        ax.yaxis.set_major_formatter(FuncFormatter(
            lambda v, _: f"{abs(v):.2f}"))
        ax.tick_params(axis="y", labelsize=6)
        ax.grid(axis="y", alpha=0.15, zorder=0)
        ax.set_axisbelow(True)

        handles = [
            Patch(facecolor=amp_dark, edgecolor="white",
                  label=f"Amp \u2014 {la}"),
            Patch(facecolor=amp_light, edgecolor="white",
                  label=f"Amp \u2014 {lb}"),
            Patch(facecolor=del_dark, edgecolor="white",
                  label=f"Del \u2014 {la}"),
            Patch(facecolor=del_light, edgecolor="white",
                  label=f"Del \u2014 {lb}"),
            Line2D([], [], marker="o", ls="", color=sig_color,
                   markeredgecolor="white", markeredgewidth=0.4,
                   label=f"cohort difference FDR < {fdr_thr:g}"),
        ]
        ax.legend(handles=handles, loc="upper center",
                  bbox_to_anchor=(0.5, -0.16), ncol=5, fontsize=6,
                  frameon=False, handlelength=1.0, columnspacing=1.0,
                  handletextpad=0.4)
        n_amp = int((sub[sub["direction"] == "Amp"]["sig"]).sum())
        n_del = int((sub[sub["direction"] == "Del"]["sig"]).sum())
        fig.suptitle(f"{comp}: {level}-level event frequency  "
                     f"(cohort-difference FDR<{fdr_thr:g}: "
                     f"amp n={n_amp}, del n={n_del})", fontsize=8)
        theme.finalize(fig, "tight")
        save(fig, cfg, f"P9_{level}_events"
             if len(comps) == 1 else f"P9_{level}_events_{comp}")


def figures_p9(res: dict, cfg) -> None:
    plot_p9_dists(res.get("scores"), res.get("tests"), cfg)
    plot_p9_forest(res.get("tests"), cfg)
    plot_p9_decomposition(res.get("scores"), res.get("tests"), cfg)
    plot_p9_events(res.get("events"), cfg, "arm")
    plot_p9_events(res.get("events"), cfg, "chrom")
    if res.get("sensitivity") is not None:
        plot_p9_sensitivity(res.get("sensitivity"),
                            res.get("sensitivity_verdict"), cfg)
