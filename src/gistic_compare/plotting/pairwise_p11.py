"""Figures for module P11 (genome-wide focal q-value landscape)."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from ..config import CHROMS, HG19_CHROM_SIZES
from ._core import *  # noqa: F401,F403

_THR_COLOR = theme.OKABE["bluish_green"]
_GENOME_LEN = float(sum(HG19_CHROM_SIZES.values()))
# The q axis is linear in -log10(q) up to q = 10**-LINTHRESH and compressed
# (symlog) beyond: a single 10**-26 deletion must not squash the peaks between
# 10**-2 and 10**-8 that carry the biology.
_SYMLOG_LINTHRESH = 5.0


def _cohort_shades(n: int) -> list[float]:
    """Blend fractions from the full hue (first cohort) towards white (last).

    Two cohorts -> 0.0 (dark) and 0.5 (light): the darker one is drawn as a
    line, the lighter one as a line plus a fill, so the reader can tell them
    apart where the tracks overlap.
    """
    if n <= 1:
        return [0.0]
    return [0.5 * i / (n - 1) for i in range(n)]


def _q_ticks(maxv: float, thr: float) -> tuple[list[float], list[str]]:
    """Explicit q-value ticks: the threshold plus 10**-k up to ``maxv``.

    The tick *labels* show the q magnitude (``10^-2`` ...), while the axis
    positions come from the (compressed) -log10(q) scale.
    """
    cand = ([thr, 2, 4, 6, 8, 10]
            + [20, 30, 50, 70, 100, 200, 300, 500, 700, 1000])
    ticks, labels = [], []
    for t in cand:
        if t > maxv + 1e-9:
            continue
        if abs(t - thr) < 1e-9:
            ticks.append(float(thr))
            labels.append(f"{10 ** (-thr):g}")
        else:
            k = int(round(t))
            ticks.append(float(t))
            labels.append(f"$10^{{-{k}}}$")
    return ticks, labels


def plot_p11_focal_landscape(scores: dict, peaks: pd.DataFrame | None, cfg,
                             q_threshold: float = 0.25) -> None:
    """Amplifications above / deletions below the q threshold, per cohort.

    Two equal-height panels share the genome x-axis: the amp panel rises from
    the baseline ``q = 1`` (``-log10(q) = 0``), the del panel mirrors it
    downward.  The y axis is **not** linear: it is linear in -log10(q) down to
    ``q = 10**-5`` and compressed beyond, so a single 10**-26 deletion cannot
    hide the peaks near the threshold.  Each cohort is a line profile; the
    lighter (last) cohort is additionally filled so the reader can tell the
    two apart at a glance, and the hue still encodes the direction (Amp =
    vermillion, Del = blue).  Peak genes (q below the threshold) are labelled
    in the panel of their direction.
    """
    if not scores:
        return
    style(cfg)
    cohorts = [c for c in cfg["cohorts"] if c in scores]
    if not cohorts:
        return
    thr = float(-np.log10(q_threshold))

    d: dict[str, pd.DataFrame] = {}
    for cid in cohorts:
        s = scores[cid].copy()
        s["x"] = _genome_x(s)
        d[cid] = s

    def _tmax(col: str) -> float:
        vals = [float(np.nanmax(d[c].loc[d[c]["type"] == col,
                                       "neg_log10_q"]))
                for c in cohorts if (d[c]["type"] == col).any()]
        vals = [v for v in vals if np.isfinite(v)]
        return max(vals + [thr * 1.3])

    max_amp, max_del = _tmax("Amp"), _tmax("Del")
    panel_h = 2.2
    fig, (ax_a, ax_d) = plt.subplots(
        2, 1, sharex=True, figsize=(theme.width(cfg), panel_h * 2.0),
        height_ratios=[1.0, 1.0])

    for ax in (ax_a, ax_d):
        ax.set_yscale("symlog", base=10, linthresh=_SYMLOG_LINTHRESH)
        ax.minorticks_off()
        ax.set_xlim(0, _GENOME_LEN)
        _draw_chrom_guides(ax, d[cohorts[0]])

    shades = _cohort_shades(len(cohorts))
    fill_idx = len(cohorts) - 1          # lighter cohort is filled, others lines
    for j, cid in enumerate(cohorts):
        s = d[cid]
        for col, ax in (("Amp", ax_a), ("Del", ax_d)):
            sub = s[s["type"] == col]
            if sub.empty:
                continue
            sub = sub.sort_values("x")
            x = sub["x"].to_numpy(dtype=float)
            y = np.clip(sub["neg_log10_q"].to_numpy(dtype=float), 0.0, None)
            c = _blend(theme.OKABE["vermillion"] if col == "Amp"
                       else theme.OKABE["blue"], shades[j])
            if j == fill_idx:
                ax.fill_between(x, 0.0, y, step="mid", color=c, alpha=0.85,
                                linewidth=0, zorder=2)
            ax.plot(x, y, color=c, lw=0.6, solid_joinstyle="round",
                    zorder=3 + j)

    ax_a.set_ylim(0.0, max_amp * 1.06)
    ax_d.set_ylim(max_del * 1.06, 0.0)          # inverted: q small = deeper
    for ax in (ax_a, ax_d):
        ax.axhline(thr, color=_THR_COLOR, ls="--", lw=0.9, zorder=6)

    # q-value ticks (labels show the q magnitude, not -log10)
    for ax, mx in ((ax_a, max_amp), (ax_d, max_del)):
        tk, lb = _q_ticks(mx, thr)
        ax.set_yticks(tk)
        ax.set_yticklabels(lb, fontsize=5.5)
        ax.tick_params(axis="y", length=2)
    ax_a.text(0.004, thr, f"q = {q_threshold:g}  ",
              transform=ax_a.get_yaxis_transform(), va="bottom", ha="left",
              fontsize=5.5, color=_THR_COLOR)
    ax_a.set_ylabel("Amp\nq value", fontsize=7)
    ax_d.set_ylabel("Del\nq value", fontsize=7)

    # chromosome strip between the panels: the ticks belong to the *del*
    # panel's top edge (pointing into it, like a GISTIC landscape strip) and
    # its bottom edge carries nothing at all
    bnds = _chrom_bounds(d[cohorts[0]])
    ax_a.set_xticks([b for b, _ in bnds])
    ax_a.set_xticklabels([c.replace("chr", "") for _, c in bnds], fontsize=5)
    ax_a.tick_params(axis="x", bottom=False, labelbottom=False, top=False,
                     labeltop=False)
    ax_d.xaxis.set_ticks_position("top")
    ax_d.tick_params(axis="x", labeltop=True, labelbottom=False, top=True,
                     bottom=False, direction="in")
    ax_d.spines["top"].set_visible(True)
    ax_d.spines["bottom"].set_visible(False)

    # peak-gene labels: strongest union peaks below the threshold
    if peaks is not None and not peaks.empty and "gene" in peaks.columns:
        p = peaks[(peaks["q"] < q_threshold) & peaks["chrom"].isin(CHROMS)
                  & peaks["q"].gt(0)].copy()
        p = p.sort_values("q").drop_duplicates("gene", keep="first")
        top_k = _top_k(cfg)
        for col, ax in (("Amp", ax_a), ("Del", ax_d)):
            sub = p[p["type"] == col].sort_values("q").head(top_k)
            if sub.empty:
                continue
            sub = sub.assign(x=_genome_x(sub))
            _annotate_labeled_points(
                ax, sub["x"].to_numpy(dtype=float),
                (-np.log10(sub["q"])).to_numpy(dtype=float),
                sub["gene"].tolist(), fontsize=5.2, max_labels=top_k)

    handles = []
    for j, cid in enumerate(cohorts):
        c = _blend(theme.OKABE["vermillion"], shades[j])
        if j == fill_idx:
            handles.append(Patch(facecolor=c, edgecolor=c,
                                 label=theme.cohort_label(cfg, cid)))
        else:
            handles.append(Line2D([], [], color=c, lw=1.2,
                                  label=theme.cohort_label(cfg, cid)))
    ax_d.legend(handles=handles, loc="upper center",
                bbox_to_anchor=(0.5, -0.10), ncol=min(len(cohorts), 4),
                fontsize=6, frameon=False, handlelength=1.0)
    fig.suptitle("focal q-value landscape \u2014 "
                 + " vs ".join(theme.cohort_label(cfg, c) for c in cohorts),
                 fontsize=8)
    theme.finalize(fig, "tight")
    save(fig, cfg, "P11_focal_landscape")


def figures_p11(res: dict, cfg) -> None:
    plot_p11_focal_landscape(res.get("scores"), res.get("peaks"), cfg,
                             float(res.get("q_threshold", 0.25)))
