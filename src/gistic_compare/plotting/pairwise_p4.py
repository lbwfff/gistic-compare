"""Figures for P4/P5 (regional differential).

(Split from the original single-module ``plotting``.)
"""

from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ._core import *  # noqa: F401,F403

# ---------------------------------------------------------------------------
# P4: regional between-cohort differential
# ---------------------------------------------------------------------------

def _plot_p4_genome(bands: pd.DataFrame, cfg: dict[str, Any],
                    criterion: str, prefix: str = "P4") -> None:
    """Cytoband-level A−B Δfreq, significance by one criterion.

    criterion = "fisher" (per-event Fisher p<0.05) or "mw" (MW FDR<0.05).
    Two panels per comparison (Amp, Del); y = Δ frequency (A − B).
    """
    if bands is None or bands.empty:
        return
    from ..config import CHROM_RANK
    from ..presets import from_cfg as preset_from_cfg
    style(cfg)
    preset = preset_from_cfg(cfg)
    mark = "*" if criterion == "fisher" else "^"
    legtxt = "Fisher p<0.05" if criterion == "fisher" else "MW FDR<0.05"
    paper = [(n, c, lo, hi, d) for (n, c, lo, hi, d, w) in preset["loci"]]
    comps = list(dict.fromkeys(bands["comparison"]))
    nrows = 2 * len(comps)
    fig, axes = plt.subplots(
        nrows, 1, squeeze=False, sharex=True,
        figsize=theme.figsize(cfg, nrows=nrows, panel_h=2.1))
    ri = 0
    for comp in comps:
        la, lb, cola, colb = _sides(cfg, comp)
        g = bands[bands["comparison"] == comp].copy()
        g["_r"] = g["chrom"].map(CHROM_RANK)
        g = g.sort_values(["_r", "start"]).reset_index(drop=True)
        x = np.arange(len(g))
        bounds, prev = [], None
        for k, c in enumerate(g["chrom"]):
            if c != prev:
                bounds.append((k, c)); prev = c
        for ptype in ("Amp", "Del"):
            ax = axes[ri][0]
            ri += 1
            df = (g[f"freq_{ptype.lower()}_a"].to_numpy(float)
                  - g[f"freq_{ptype.lower()}_b"].to_numpy(float))
            if criterion == "fisher":
                pv = g[f"p_{ptype.lower()}"].to_numpy(float)
            else:
                pv = g["p_mw_fdr"].to_numpy(float)
            sig = np.isfinite(pv) & (pv < 0.05)
            pos = sig & (df > 0)
            neg = sig & (df < 0)
            for k, (pos_i, _) in enumerate(bounds):
                if k > 0:
                    ax.axvline(pos_i, color="#E5E5E5", lw=0.7, zorder=0)
            ax.axhline(0, color="#666666", lw=0.6, zorder=2)
            ax.scatter(x[~sig], df[~sig], s=8, color="#CCCCCC", zorder=3,
                       label="n.s.")
            ax.scatter(x[pos], df[pos], s=18, color=cola, zorder=4,
                       label=f"{la} more ({legtxt})")
            ax.scatter(x[neg], df[neg], s=18, color=colb, zorder=4,
                       label=f"{lb} more ({legtxt})")
            band_marks = [(nm, c, l, h) for (nm, c, l, h, dd) in paper
                          if dd == ptype]
            for k, (name, chrom, lo, hi) in enumerate(band_marks):
                m = ((g["chrom"] == chrom) & (g["end"] > lo)
                     & (g["start"] < hi)).to_numpy()
                if not m.any():
                    continue
                lo_i = int(np.argmax(m))
                hi_i = int(len(m) - np.argmax(m[::-1]) - 1)
                ax.axvspan(lo_i - 0.5, hi_i + 0.5, color="#BDBDBD",
                           alpha=0.35, zorder=1)
                mid = (lo_i + hi_i) // 2
                if not bool(sig[mid]):
                    continue
                dv = float(df[mid])
                col = cola if dv > 0 else colb
                up = dv > 0
                dx = 10 if k % 2 == 0 else -10
                dy = 13 if up else -13
                ax.annotate(f"{name}{mark}", xy=(x[mid], dv), xytext=(dx, dy),
                            textcoords="offset points",
                            ha="left" if dx > 0 else "right",
                            va="bottom" if up else "top",
                            fontsize=5.5, color=col,
                            arrowprops=dict(arrowstyle="-", lw=0.5, color=col),
                            annotation_clip=True)
            ax.set_title(f"{comp} | {ptype}", fontsize=8)
            ax.set_ylabel(f"Δ freq ({la} − {lb})")
            ax.tick_params(labelsize=6)
            ax.set_xlim(-1, len(g))
            ymax = max(0.06, float(np.nanmax(np.abs(df))) * 1.35)
            ax.set_ylim(-ymax, ymax)
            npair = _n_pair(cfg, comp)
            theme.annot_n(ax, f"{npair} samples" if npair else "",
                          loc="lower right")
            theme.panel_label(ax, theme.panel_letter(ri - 1), cfg)
        axes[ri - 1][0].set_xticks([p_ for p_, _ in bounds])
        axes[ri - 1][0].set_xticklabels([c.replace("chr", "") for _, c in bounds],
                                        fontsize=6)
        axes[ri - 1][0].set_xlabel("chromosome")
    h, l = axes[0][0].get_legend_handles_labels()
    if h:
        axes[0][0].legend(h, l, frameon=False, fontsize=5.5, loc="upper left")
    theme.finalize(fig, "constrained", hspace=0.16)
    save(fig, cfg, f"{prefix}_genome_diff_{criterion}")


def plot_p4_genome_diff_fisher(bands, cfg, prefix="P4"):
    _plot_p4_genome(bands, cfg, "fisher", prefix=prefix)


def plot_p4_genome_diff_mw(bands, cfg, prefix="P4"):
    _plot_p4_genome(bands, cfg, "mw", prefix=prefix)


def plot_p4_null(nullcal, summary, cfg, prefix="P4"):
    """Observed vs random-split null count of FDR-significant 1 Mb windows."""
    if nullcal is None or nullcal.empty:
        return
    style(cfg)
    coh = list(nullcal["cohort"])
    x = np.arange(len(coh))
    pal = _palette(cfg)
    # observed significant windows for the comparison each cohort takes part in
    obs = []
    sum_n = (summary["n_sig_windows"].to_numpy(dtype=float)
             if summary is not None and not summary.empty
             and "n_sig_windows" in summary.columns else np.array([]))
    comps_sum = (list(summary["comparison"]) if summary is not None
                 and not summary.empty else [])
    for c in coh:
        val = np.nan
        idx = [i for i, comp in enumerate(comps_sum)
               if c in theme.comp_sides(cfg, comp)]
        if idx and len(sum_n):
            val = float(np.nanmax(sum_n[idx]))
        obs.append(val)
    obs = np.asarray(obs, dtype=float)
    null_mean = nullcal["sig_windows_mean"].to_numpy(dtype=float)
    null_hi = nullcal["sig_windows_p95"].to_numpy(dtype=float)
    w = 0.34
    fig, ax = plt.subplots(figsize=theme.figsize(cfg, panel_h=3.0,
                                                 total_w=4.8))
    for xi, c in enumerate(coh):
        if np.isfinite(obs[xi]):
            ax.bar(xi - w / 2 - 0.01, obs[xi], w, color=pal.get(c, "#333333"),
                   edgecolor="white", linewidth=0.4)
            if obs[xi] > 0:
                ax.text(xi - w / 2 - 0.01, obs[xi], f"{int(obs[xi])}",
                        ha="center", va="bottom", fontsize=5.5,
                        color="#444444")
    ax.bar(x + w / 2 + 0.01, null_mean, w, color="#CCCCCC",
           edgecolor="#999999", linewidth=0.4)
    ax.errorbar(x + w / 2 + 0.01, null_mean,
                yerr=[np.zeros(len(x)),
                      np.clip(null_hi - null_mean, 0, None)],
                fmt="none", ecolor="#666666", elinewidth=0.8, capsize=2.5,
                capthick=0.8, zorder=5)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{theme.cohort_label(cfg, c)}\n"
                        f"(half of {int(r.n_samples)})"
                        for c, r in zip(coh, nullcal.itertuples())],
                       fontsize=6.5)
    ax.set_ylabel("FDR-significant 1 Mb windows")
    ymax = 0.0
    for arr in (obs, null_mean, null_hi):
        if len(arr) and np.isfinite(arr).any():
            ymax = max(ymax, float(np.nanmax(arr)))
    if ymax <= 0:
        # an all-zero panel otherwise autoscales to a nonsensical ±range
        ax.set_ylim(0, 1)
        theme.annot_stat(ax, "no FDR-significant windows\nin any split",
                         loc="upper left", fontsize=6.5)
    else:
        ax.set_ylim(0, ymax * 1.28)
    from matplotlib.patches import Patch
    ax.legend(handles=[
        Patch(facecolor="#BBBBBB", edgecolor="white", label="observed"),
        Patch(facecolor="#CCCCCC", edgecolor="#999999",
              label="null (random split, mean, p95)"),
    ], frameon=False, fontsize=6.5, loc="upper right")
    theme.panel_label(ax, theme.panel_letter(0), cfg)
    theme.finalize(fig, "tight")
    save(fig, cfg, f"{prefix}_null_calibration")


def plot_p4_paper_loci(paper: pd.DataFrame, cfg: dict[str, Any],
                       prefix: str = "P4") -> None:
    """Alteration frequency at the paper's prior loci, cohort A vs cohort B."""
    if paper is None or paper.empty:
        return
    style(cfg)
    loci = list(dict.fromkeys(paper["locus"]))
    comp = str(paper["comparison"].iloc[0])
    la, lb, ca, cb = _sides(cfg, comp)
    fig, axes = plt.subplots(1, len(loci), figsize=theme.figsize(
        cfg, panel_h=3.2, total_w=theme.FULL_WIDTH), sharey=True)
    if len(loci) == 1:
        axes = [axes]
    w = 0.34
    # with sharey the limits must be common to every panel, otherwise the last
    # panel sets them and taller panels' labels/CIs escape the axes box
    gtop = 0.05
    for _, r in paper.iterrows():
        for v in (r["freq_amp_a"], r["freq_amp_b"],
                  r["freq_del_a"], r["freq_del_b"]):
            if np.isfinite(v):
                gtop = max(gtop, float(v))
        for v, n_ in ((r["freq_amp_a"], r["n_a"]), (r["freq_amp_b"], r["n_b"]),
                      (r["freq_del_a"], r["n_a"]), (r["freq_del_b"], r["n_b"])):
            if np.isfinite(v) and np.isfinite(n_) and n_ > 0:
                _lo, _hi = theme.wilson_ci([v * n_], [n_])
                if np.isfinite(_hi[0]):
                    gtop = max(gtop, float(_hi[0]))
    ylim_top = gtop * 1.50
    for k, (ax, loc) in enumerate(zip(axes, loci)):
        r = paper[paper["locus"] == loc].iloc[0]
        na = float(r.get("n_a", np.nan))
        nb = float(r.get("n_b", np.nan))
        groups = [("Amp", float(r["freq_amp_a"]), float(r["freq_amp_b"])),
                  ("Del", float(r["freq_del_a"]), float(r["freq_del_b"]))]
        for gi, (pt, fa, fb) in enumerate(groups):
            la_lo, la_hi = theme.wilson_ci([fa * na], [na])
            lb_lo, lb_hi = theme.wilson_ci([fb * nb], [nb])
            xa, xb = gi - w / 2 - 0.01, gi + w / 2 + 0.01
            theme.bar_with_ci(
                ax, [xa], [fa], [float(la_lo[0])], [float(la_hi[0])],
                width=w, color=ca, edgecolor="white", linewidth=0.4,
                label=f"{la} (n = {int(na)})" if gi == 0 else None)
            theme.bar_with_ci(
                ax, [xb], [fb], [float(lb_lo[0])], [float(lb_hi[0])],
                width=w, color=cb, edgecolor="white", linewidth=0.4,
                label=f"{lb} (n = {int(nb)})" if gi == 0 else None)
            for xpos, v, hi in ((xa, fa, float(la_hi[0])),
                                (xb, fb, float(lb_hi[0]))):
                if v > 0:
                    ax.text(xpos, max(v, hi), f"{v:.3f}", ha="center",
                            va="bottom", fontsize=5, color="#444444")
        ax.set_xticks(np.arange(len(groups)))
        ax.set_xticklabels([g[0] for g in groups])
        ax.set_title(loc, fontsize=8)
        ax.set_ylim(0, ylim_top)
        ax.grid(False)
        ax.tick_params(labelsize=6.5)
        theme.annot_stat(ax, "\n".join([
            "MW " + theme.fmt_p(float(r["p_mw"])),
            f"Amp {theme.fmt_p(float(r['p_amp']))}, "
            f"Del {theme.fmt_p(float(r['p_del']))}",
        ]), loc="upper left", fontsize=5.0)
        theme.panel_label(ax, theme.panel_letter(k), cfg)
    axes[0].set_ylabel("fraction altered")
    # one figure-level legend: per-axes legends collide with the p-value text
    h, lab = axes[0].get_legend_handles_labels()
    if h:
        fig.legend(h, lab, frameon=False, fontsize=6.5, ncol=2,
                   loc="upper center", bbox_to_anchor=(0.5, 0.99))
    # paper expectation belongs in a footnote, not in the panel titles
    if {"paper_direction", "expected_cohort"} <= set(paper.columns):
        foot = "; ".join(
            f"{r['locus']} {r['paper_direction']} {r['expected_cohort']}"
            for _, r in paper.drop_duplicates("locus").iterrows())
        fig.text(0.5, 0.015, "paper expectation \u2014 " + foot,
                 fontsize=5.5, color="#555555", ha="center", va="bottom")
    # reserve a strip for the legend and the footnote; tight_layout alone
    # ignores figure-level artists and would let them overlap the axes
    fig.tight_layout(rect=(0.0, 0.06, 1.0, 0.90))
    save(fig, cfg, f"{prefix}_paper_loci")


def _paper_rows(bands: pd.DataFrame, comp: str,
                cfg: dict[str, Any] | None = None) -> list[dict]:
    from ..presets import from_cfg as preset_from_cfg
    preset = preset_from_cfg(cfg)
    b = bands[bands["comparison"] == comp].drop_duplicates("band").set_index("band")
    rows = []
    for name, star in preset["bands"]:
        if name not in b.index:
            continue
        r = b.loc[name]
        rows.append(dict(name=name, star=star,
                         da=float(r["diff_amp"]), dd=float(r["diff_del"]),
                         pa=r["p_amp"], pd=r["p_del"],
                         mw=r["p_mw"], mwf=r["p_mw_fdr"]))
    return rows


def _paper_set_ylabels(ax, rows: list[dict]) -> None:
    ax.set_yticklabels([r["name"] for r in rows], fontsize=6)
    for t, r in zip(ax.get_yticklabels(), rows):
        if r["star"]:
            t.set_fontweight("bold"); t.set_color("#B2182B")


def plot_p4_paper_list_heatmap_fisher(bands, cfg, prefix="P4"):
    """Paper Table-2 bands × {Amp,Del}: Δ freq (A−B) + Fisher p<0.05 marker."""
    if bands is None or bands.empty:
        return
    style(cfg)
    comp = str(bands["comparison"].iloc[0])
    rows = _paper_rows(bands, comp, cfg)
    if not rows:
        return
    mat = np.array([[r["da"], r["dd"]] for r in rows], dtype=float)
    txt = [["*" if (pd.notna(r["pa"]) and r["pa"] < 0.05) else "",
            "*" if (pd.notna(r["pd"]) and r["pd"] < 0.05) else ""]
           for r in rows]
    fig, ax = plt.subplots(figsize=theme.figsize(
        cfg, panel_h=max(5.0, 0.25 * len(rows) + 1.8), total_w=4.6))
    cmap = sns.diverging_palette(240, 25, as_cmap=True)
    im = ax.imshow(mat, cmap=cmap, aspect="auto", vmin=-0.2, vmax=0.2)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Amp", "Del"], fontsize=7)
    ax.set_yticks(range(len(rows))); _paper_set_ylabels(ax, rows)
    for i in range(len(rows)):
        for j in range(2):
            if txt[i][j]:
                ax.text(j, i, txt[i][j], ha="center", va="center",
                        fontsize=11, color="#222222", fontweight="bold")
    for i in range(len(rows) + 1):
        ax.axhline(i - 0.5, color="white", lw=0.8)
    for j in (-0.5, 0.5, 1.5):
        ax.axvline(j, color="white", lw=0.8)
    ax.set_title(f"{comp}: paper Table-2 bands", fontsize=8)
    ax.set_xlabel(""); ax.set_ylabel("")
    cbar = fig.colorbar(im, ax=ax, shrink=0.5)
    cbar.set_label("Δ freq (A − B)", fontsize=7)
    fig.supxlabel("* Fisher p<0.05; bold red = ethnic-specific",
                  fontsize=5.5, color="#555555")
    theme.finalize(fig, "tight")
    save(fig, cfg, f"{prefix}_paper_list_heatmap_fisher")


def plot_p4_paper_list_heatmap_mw(bands, cfg, prefix="P4"):
    """Paper Table-2 bands × Mann-Whitney: −log10(FDR) + FDR<0.05 marker."""
    if bands is None or bands.empty:
        return
    style(cfg)
    comp = str(bands["comparison"].iloc[0])
    rows = _paper_rows(bands, comp, cfg)
    if not rows:
        return
    val = np.array([[-np.log10(max(float(r["mwf"]), 1e-300))]
                    for r in rows], dtype=float)
    txt = [["*" if (pd.notna(r["mwf"]) and r["mwf"] < 0.05) else ""]
           for r in rows]
    vmax = max(1.4, float(np.nanmax(val)))
    fig, ax = plt.subplots(figsize=theme.figsize(
        cfg, panel_h=max(5.0, 0.25 * len(rows) + 1.8), total_w=3.2))
    im = ax.imshow(val, cmap="Blues", aspect="auto", vmin=0, vmax=vmax)
    ax.set_xticks([0]); ax.set_xticklabels(["MW"], fontsize=7)
    ax.set_yticks(range(len(rows))); _paper_set_ylabels(ax, rows)
    for i in range(len(rows)):
        if txt[i][0]:
            ax.text(0, i, txt[i][0], ha="center", va="center",
                    fontsize=11, color="#222222", fontweight="bold")
    for i in range(len(rows) + 1):
        ax.axhline(i - 0.5, color="white", lw=0.8)
    for j in (-0.5, 0.5):
        ax.axvline(j, color="white", lw=0.8)
    ax.set_title(f"{comp}: paper Table-2 bands", fontsize=8)
    ax.set_xlabel(""); ax.set_ylabel("")
    cbar = fig.colorbar(im, ax=ax, shrink=0.5)
    cbar.set_label("−log10 FDR (MW)", fontsize=7)
    fig.supxlabel("* FDR<0.05; bold red = ethnic-specific",
                  fontsize=5.5, color="#555555")
    theme.finalize(fig, "tight")
    save(fig, cfg, f"{prefix}_paper_list_heatmap_mw")


def figures_p4(res: dict, cfg) -> None:
    bands = res.get("bands")
    plot_p4_genome_diff_fisher(bands, cfg)
    plot_p4_genome_diff_mw(bands, cfg)
    plot_p4_paper_loci(res.get("paper"), cfg)
    plot_p4_paper_list_heatmap_fisher(bands, cfg)
    plot_p4_paper_list_heatmap_mw(bands, cfg)
    plot_p4_null(res.get("null"), res.get("summary"), cfg)


def figures_p5(res: dict, cfg) -> None:
    """Same figures as P4, but on GISTIC's marker matrix (P5_ prefix)."""
    bands = res.get("bands")
    plot_p4_genome_diff_fisher(bands, cfg, prefix="P5")
    plot_p4_genome_diff_mw(bands, cfg, prefix="P5")
    plot_p4_paper_loci(res.get("paper"), cfg, prefix="P5")
    plot_p4_paper_list_heatmap_fisher(bands, cfg, prefix="P5")
    plot_p4_paper_list_heatmap_mw(bands, cfg, prefix="P5")
    plot_p4_null(res.get("null"), res.get("summary"), cfg, prefix="P5")


