"""Figures for module P1 (peak set differences).

(Split from the original single-module ``plotting``.)
"""

from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ._core import *  # noqa: F401,F403
from .benchmark_a import plot_three_way
from .benchmark_b import plot_track_correlations

# ---------------------------------------------------------------------------
# P1 / P2 / P3: pairwise workflow
# ---------------------------------------------------------------------------

def plot_p1_jaccard(jaccard: pd.DataFrame, cfg: dict[str, Any]) -> None:
    if jaccard.empty:
        return
    style(cfg)
    order = ["moderate", "strict", "loose"]
    sub = jaccard[jaccard["tier"].isin(order)]
    if sub.empty:
        return
    g = sns.catplot(
        data=sub, x="comparison", y="jaccard", hue="tier",
        col="ptype", kind="bar", order=sorted(sub["comparison"].unique()),
        hue_order=order, height=2.7, aspect=0.95,
        palette="Blues_r", errorbar=None, edgecolor="white",
    )
    comps = sorted(sub["comparison"].unique())
    tick = []
    for c in comps:
        s = sub[sub["comparison"] == c]
        if len(s) and "n_a" in s.columns:
            tick.append(f"{c}\n({int(s['n_a'].iloc[0])}/"
                        f"{int(s['n_b'].iloc[0])}, "
                        f"{int(s['n_pairs'].iloc[0])})")
        else:
            tick.append(c)
    g.set_xticklabels(tick, rotation=20, ha="right", fontsize=6)
    g.set_axis_labels("comparison  (n = peaks A / peaks B, matched)",
                      "Jaccard")
    for k, ax in enumerate(g.axes.flat):
        theme.panel_label(ax, theme.panel_letter(k), cfg)
    g.legend.set_title("tier")
    theme.finalize(g.figure, "constrained")
    save(g.figure, cfg, "P1_jaccard")


def plot_p1_delta_freq(delta: pd.DataFrame, delta_stats: pd.DataFrame,
                       cfg: dict[str, Any]) -> None:
    if delta is None or delta.empty:
        return
    style(cfg)
    sub = delta.copy()
    n_comp = sub["comparison"].nunique()
    fig, axes = plt.subplots(1, 2, sharey=True,
                             figsize=theme.figsize(
                                 cfg, panel_h=3.0 + 0.35 * max(0, n_comp - 1)))
    for k, (ax, ptype) in enumerate(zip(axes, ("Amp", "Del"))):
        g = sub[sub["ptype"] == ptype]
        if g.empty:
            ax.set_visible(False)
            continue
        comps = sorted(g["comparison"].unique())
        data_box = [g.loc[g["comparison"] == c, "delta_freq"].to_numpy(float)
                    for c in comps]
        tick = []
        for c in comps:
            la, lb, _, _ = _sides(cfg, c)
            tick.append(f"{la}\nvs {lb}")
        bp = ax.boxplot(data_box, tick_labels=tick, widths=0.55,
                        patch_artist=True, showfliers=False)
        col = _pc(ptype)
        for patch in bp["boxes"]:
            patch.set_facecolor(col)
            patch.set_alpha(0.45)
        rng = np.random.default_rng(42)
        for i, arr in enumerate(data_box):
            x = rng.normal(i + 1, 0.045, len(arr))
            ax.scatter(x, arr, s=10, alpha=0.7, color=col,
                       edgecolors="none", zorder=3)
        ax.axhline(0, color="#666666", ls="--", lw=0.8)
        ax.set_title(ptype, fontsize=8)
        ax.set_ylabel("Δ frequency (A − B)")
        if not delta_stats.empty:
            sg = delta_stats[delta_stats["ptype"] == ptype]
            stats_txt = []
            for _, r in sg.iterrows():
                n_pr = int(r["n_pairs"]) if "n_pairs" in sg.columns else None
                kind = _kind(cfg, r["comparison"])
                stats_txt.append(
                    f"{r['comparison']}"
                    + (f" [{kind}]" if kind else "")
                    + f": med {r['median_delta_freq']:+.2f}, "
                    f"{theme.fmt_p(r.get('pvalue', np.nan))}"
                    + (f", n = {n_pr}" if n_pr is not None else ""))
            if stats_txt:
                # reserve headroom first: the block is wide and would
                # otherwise be drawn across the box and whiskers
                y0, y1 = ax.get_ylim()
                span = max(y1 - y0, 1e-9)
                ax.set_ylim(y0 - 0.06 * span,
                            y1 + 0.06 * span
                            + 0.26 * span * max(1, len(stats_txt)))
                theme.annot_stat(ax, "\n".join(stats_txt), loc="upper left",
                                 fontsize=5.5)
        theme.panel_label(ax, theme.panel_letter(k), cfg)
    theme.finalize(fig, "tight")
    save(fig, cfg, "P1_delta_freq")


def plot_p1_upset(membership: dict[str, pd.DataFrame],
                  cfg: dict[str, Any]) -> None:
    """Multi-cohort peak-key membership UpSet when >= 3 cohorts are listed."""
    if len(cfg.get("cohorts", {})) < 3:
        return
    rows = []
    for ptype, df in (membership or {}).items():
        if df is None or df.empty or "peak_key" not in df.columns:
            continue
        for key, g in df.groupby("peak_key"):
            rows.append(dict(type=ptype,
                             cohorts=",".join(sorted(set(g["cohort"])))))
    if not rows:
        return
    upset_input = pd.DataFrame(rows)
    if upset_input.empty:
        return
    plot_three_way(upset_input, cfg, prefix="P1",
                   name_pattern="{prefix}_upset_{ptype_l}")


def plot_p1_track_corr(corr: pd.DataFrame, cfg: dict[str, Any]) -> None:
    """Pearson + Spearman heatmaps (P1-prefixed, no file rename)."""
    if corr is None or corr.empty:
        return
    for method in ("pearson", "spearman"):
        plot_track_correlations(corr, cfg, method=method, prefix="P1")














def plot_p1_unique_genome(unique: pd.DataFrame, cfg: dict[str, Any]) -> None:
    """Unique (one-sided) peaks along the genome: x = position, y = −log10 q."""
    if unique is None or unique.empty:
        return
    style(cfg)
    u = unique.copy()
    u["x"] = _genome_x(u)
    comps = sorted(u["comparison"].unique())
    # side a = first cohort of comparison, side b = second
    fig, axes = plt.subplots(
        len(comps), 2, squeeze=False,
        figsize=theme.figsize(cfg, total_w=theme.FULL_WIDTH,
                              panel_h=3.0 * max(1, len(comps))))
    # genome-wide chromosome tick positions (all chrs, not just observed)
    from ..config import CHROMS, HG19_CHROM_SIZES
    ticks, labels = [], []
    acc = 0.0
    for c in CHROMS:
        ticks.append(acc + HG19_CHROM_SIZES[c] / 2)
        labels.append(c.replace("chr", ""))
        acc += HG19_CHROM_SIZES[c]
    x_max = acc
    for i, comp in enumerate(comps):
        ca, cb = comp.split("_vs_") if "_vs_" in comp else (comp, "")
        for j, (side, who) in enumerate((("a", ca), ("b", cb))):
            ax = axes[i, j]
            g = u[(u["comparison"] == comp) & (u["side"] == side)]
            if g.empty:
                ax.set_visible(False)
                continue
            _draw_chrom_guides(ax, g)
            xs_all, ys_all, labs_all = [], [], []
            for ptype in ("Amp", "Del"):
                s = g[g["ptype"] == ptype]
                if s.empty:
                    continue
                ax.scatter(s["x"], s["neglogq"], s=14, alpha=0.8,
                           color=_pc(ptype), label=ptype,
                           edgecolors="none", zorder=3)
                xs_all.extend(s["x"].tolist())
                ys_all.extend(s["neglogq"].tolist())
                labs_all.extend(
                    _label_series(s, ["cytoband"], "peak_id"))
            # label peaks (cytoband), strongest first, budgeted
            _annotate_labeled_points(ax, xs_all, ys_all, labs_all,
                                     fontsize=5.5, min_sep_pt=20,
                                     max_labels=_top_k(cfg))
            ax.set_xlim(0, x_max * 1.01)
            ax.set_ylim(bottom=0)
            ax.set_title(f"only in {theme.cohort_label(cfg, who)}",
                         fontsize=7.5)
            theme.annot_n(ax, f"n = {len(g)} peaks", loc="upper right")
            if i == len(comps) - 1:
                ax.set_xticks(ticks)
                ax.set_xticklabels(labels, fontsize=5.5, rotation=90)
                ax.set_xlabel("chromosome")
            else:
                ax.set_xticklabels([])
            if j == 0:
                ax.set_ylabel("−log10 q")
            if i == 0 and j == 0:
                ax.legend(frameon=False, fontsize=6, loc="upper left")
            ax.tick_params(labelsize=6)
            theme.panel_label(ax, theme.panel_letter(i * 2 + j), cfg)
    theme.finalize(fig, "constrained", hspace=0.4, wspace=0.16)
    save(fig, cfg, "P1_unique_genome")


def plot_p1_matched_scatter(matched: pd.DataFrame, cfg: dict[str, Any]) -> None:
    """Matched-pair freq–freq and −log10 q–q scatters (A vs B)."""
    if matched is None or matched.empty:
        return
    style(cfg)
    m = matched.copy()
    m["a_nlq"] = -np.log10(np.clip(m["a_q"].to_numpy(float), 1e-300, None))
    m["b_nlq"] = -np.log10(np.clip(m["b_q"].to_numpy(float), 1e-300, None))
    comps = list(dict.fromkeys(m["comparison"]))
    ptypes = ["Amp", "Del"]
    metrics = [("a_freq", "b_freq", "frequency altered", "\u03c1"),
               ("a_nlq", "b_nlq", "\u2212log10 q", "\u03c1")]
    fig, axes = plt.subplots(
        len(comps), 2, squeeze=False, sharex=False, sharey=False,
        figsize=theme.figsize(cfg, total_w=theme.FULL_WIDTH,
                              panel_h=3.3 * len(comps)))
    for i, comp in enumerate(comps):
        la_lab, lb_lab, _, _ = _sides(cfg, comp)
        for k, (xc, yc, lab, sym) in enumerate(metrics):
            ax = axes[i, k]
            g = m[m["comparison"] == comp]
            if g.empty:
                ax.set_visible(False)
                continue
            for ptype in ptypes:
                s = g[g["ptype"] == ptype]
                if s.empty:
                    continue
                ax.scatter(s[xc], s[yc], s=18, alpha=0.8, color=_pc(ptype),
                           label=ptype, edgecolors="none", zorder=3)
            lo = float(np.nanmin(g[[xc, yc]].to_numpy(float)))
            hi = float(np.nanmax(g[[xc, yc]].to_numpy(float)))
            pad = 0.05 * (hi - lo + 1e-9)
            ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad],
                    color="#888888", ls="--", lw=0.7, zorder=1)
            if k == 0:
                ax.set_xlim(-0.02, 1.02)
                ax.set_ylim(-0.02, 1.02)
            xv = g[xc].to_numpy(float)
            yv = g[yc].to_numpy(float)
            ok = np.isfinite(xv) & np.isfinite(yv)
            if ok.sum() >= 3:
                r, p = stats.spearmanr(xv[ok], yv[ok])
                theme.annot_stat(
                    ax, theme.fmt_ci(sym, float(r), n=int(ok.sum()),
                                     p=float(p)),
                    loc="lower right")
            if i == 0:
                ax.set_title(f"{lab} | {comp}", fontsize=7.5)
                theme.panel_label(ax, theme.panel_letter(k), cfg)
            if i == len(comps) - 1:
                ax.set_xlabel(f"{la_lab}: {lab}")
            if k == 0 and i == 0:
                ax.legend(frameon=False, fontsize=6, loc="upper left")
            if k == 0:
                ax.set_ylabel(f"{lb_lab}: {lab}")
            ax.tick_params(labelsize=6)
    theme.finalize(fig, "constrained", hspace=0.3, wspace=0.25)
    save(fig, cfg, "P1_matched_scatter")


def plot_p1_set_summary(jaccard: pd.DataFrame, cfg: dict[str, Any]) -> None:
    """Detected-peak budget per comparison, primary (moderate) tier.

    The two bars of a comparison are cohort A and cohort B, drawn in their
    cohort colours and split into peaks shared by both sides (solid) and
    peaks unique to that side (hatched). Each bar is annotated with
    ``matched/unique``.
    """
    if jaccard is None or jaccard.empty:
        return
    style(cfg)
    primary = jaccard[jaccard["tier"] == "moderate"].copy()
    if primary.empty:
        primary = jaccard.copy()
    primary["only_a"] = primary["n_a"] - primary["n_pairs"]
    primary["only_b"] = primary["n_b"] - primary["n_pairs"]
    comps = sorted(primary["comparison"].unique())
    fig, axes = plt.subplots(1, 2, sharey=True,
                             figsize=theme.figsize(cfg, panel_h=3.4))
    for k, (ax, ptype) in enumerate(zip(axes, ("Amp", "Del"))):
        g = (primary[primary["ptype"] == ptype]
             .set_index("comparison").reindex(comps))
        if g.empty or g["n_a"].isna().all():
            ax.set_visible(False)
            continue
        x = np.arange(len(comps))
        w = 0.30
        gap = 0.06
        ymax = float(np.nanmax(g[["n_a", "n_b"]].to_numpy(dtype=float)))
        tick_pos, tick_lab = [], []
        for xi, comp in enumerate(comps):
            la, lb, ca, cb = _sides(cfg, comp)
            row = g.loc[comp]
            matched = float(row["n_pairs"])
            for di, (uniq_col, side_col, side_lab) in enumerate(
                    (("only_a", ca, la), ("only_b", cb, lb))):
                uniq = float(row[uniq_col])
                xpos = x[xi] + (di - 0.5) * (w + gap)
                ax.bar(xpos, matched, w, color=side_col, edgecolor="white",
                       linewidth=0.4)
                ax.bar(xpos, uniq, w, bottom=matched, color=side_col,
                       alpha=0.30, edgecolor="#555555", linewidth=0.4,
                       hatch="//")
                ax.text(xpos, matched + uniq + ymax * 0.025,
                        f"{int(matched)}/{int(uniq)}", ha="center",
                        va="bottom", fontsize=5.5)
                tick_pos.append(xpos)
                tick_lab.append(f"{side_lab}\n({comp})")
        ax.set_xticks(tick_pos)
        ax.set_xticklabels(tick_lab, rotation=30, ha="right", fontsize=6)
        ax.set_title(ptype, fontsize=8)
        if k == 0:
            ax.set_ylabel("n peaks")
        ax.set_xlim(-0.6, len(comps) - 0.4)
        ax.set_ylim(0, ymax * 1.22)
        # legend from proxies so it is identical in both panels
        from matplotlib.patches import Patch
        ax.legend(handles=[
            Patch(facecolor="#BBBBBB", edgecolor="white", label="matched"),
            Patch(facecolor="#BBBBBB", alpha=0.30, edgecolor="#555555",
                  hatch="//", label="unique to that cohort"),
        ], frameon=False, fontsize=6, loc="upper left")
        theme.panel_label(ax, theme.panel_letter(k), cfg)
    theme.finalize(fig, "tight")
    save(fig, cfg, "P1_set_summary")


def plot_p1_delta_track(delta_grid: pd.DataFrame, cfg: dict[str, Any]) -> None:
    """Genome-wide Δ(−log10 q) = A − B filled tracks (red >0, blue <0)."""
    if delta_grid is None or delta_grid.empty:
        return
    style(cfg)
    df = delta_grid.reset_index(drop=True)
    if {"chrom", "start", "end"} > set(df.columns):
        pass
    cols = [c for c in df.columns
            if c not in ("chrom", "start", "end") and "_" in c]
    # columns named {comp}_{ptype} where comp contains _vs_
    amp_cols = [c for c in cols if c.endswith("_Amp")]
    del_cols = [c for c in cols if c.endswith("_Del")]
    if not amp_cols and not del_cols:
        return
    x = np.arange(len(df))
    chrom_bounds = []
    prev = None
    for i, c in enumerate(df["chrom"]):
        if c != prev:
            chrom_bounds.append((i, c))
            prev = c

    for ptype, use_cols in (("Amp", amp_cols), ("Del", del_cols)):
        if not use_cols:
            continue
        n = len(use_cols)
        fig, axes = plt.subplots(
            n, 1, figsize=theme.figsize(cfg, nrows=n, panel_h=2.4),
            sharex=True, squeeze=False)
        vmax = np.nanmax([np.nanmax(np.abs(df[c].to_numpy(float)))
                          for c in use_cols]) or 1.0
        for i, col in enumerate(use_cols):
            ax = axes[i, 0]
            comp = col.rsplit("_", 1)[0]
            la, lb, ca, cb = _sides(cfg, comp)
            vals = df[col].to_numpy(float)
            pos = np.where(vals > 0, vals, 0.0)
            neg = np.where(vals < 0, vals, 0.0)
            ax.fill_between(x, 0, pos, color=ca, alpha=0.9,
                            linewidth=0, zorder=2, label=f"{la} higher")
            ax.fill_between(x, 0, neg, color=cb, alpha=0.9,
                            linewidth=0, zorder=2, label=f"{lb} higher")
            ax.set_ylim(-vmax * 1.05, vmax * 1.05)
            ax.axhline(0, color="#666666", lw=0.6, zorder=3)
            ax.set_ylabel(comp, fontsize=6.5)
            if i == 0:
                ax.legend(frameon=False, fontsize=6, loc="upper right",
                          ncol=2)
            theme.panel_label(ax, theme.panel_letter(i), cfg)
            for k, (pos_i, _) in enumerate(chrom_bounds):
                if k % 2 == 1:
                    end = (chrom_bounds[k + 1][0] if k + 1 < len(chrom_bounds)
                           else len(df))
                    ax.axvspan(pos_i, end, color="#F2F2F2", zorder=0)
                if k > 0:
                    ax.axvline(pos_i, color="#CCCCCC", lw=0.5, zorder=1)
            ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
            ax.tick_params(labelsize=6)
        axes[-1, 0].set_xticks([p for p, _ in chrom_bounds])
        axes[-1, 0].set_xticklabels(
            [c.replace("chr", "") for _, c in chrom_bounds], fontsize=6)
        axes[-1, 0].set_xlabel("chromosome")
        fig.supylabel(f"Δ −log10 q ({ptype}: A − B)", fontsize=8)
        theme.finalize(fig, "constrained", hspace=0.35)
        save(fig, cfg, f"P1_delta_track_{ptype.lower()}")


def figures_p1(res: dict, cfg) -> None:
    plot_p1_jaccard(res.get("jaccard"), cfg)
    plot_p1_delta_freq(res.get("delta"), res.get("delta_stats"), cfg)
    plot_p1_track_corr(res.get("correlations"), cfg)
    plot_p1_upset(res.get("membership"), cfg)
    plot_p1_unique_genome(res.get("unique"), cfg)
    plot_p1_matched_scatter(res.get("matched"), cfg)
    plot_p1_set_summary(res.get("jaccard"), cfg)
    plot_p1_delta_track(res.get("delta_grid"), cfg)
    from .oncoprint import figures_oncoprint
    figures_oncoprint(res, cfg)
    from .oncoplot import figures_oncoplot as _figures_oncoplot
    _figures_oncoplot(res, cfg)


