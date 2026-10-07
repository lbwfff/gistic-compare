"""Union-peak OncoPrint: per-sample GISTIC calls across cohorts.

Rows are union peak keys (P1 moderate-match chains), columns are samples
grouped by cohort; cells show the 0/1/2 amplitude-threshold call.  Grey
means the cohort has **no member peak** for that key (GISTIC did not call
a peak there) — distinct from a white "called, no alteration" cell.

Layout is fully manual (no tight_layout) because cohort labels, panel type
tags and the footnote live outside the axes.
"""

from __future__ import annotations

import traceback
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap

from .. import theme
from ._core import save, style

# amplitude-threshold colours, index = call value 0/1/2
_AMP_COLORS = ["#ffffff", "#f2a68c", "#d55e00"]
_DEL_COLORS = ["#ffffff", "#8ab8d8", "#0072b2"]
_MISSING = "#e3e3e3"

_MAX_LABEL_ROWS_FULL = 45   # below this: label every row
_LABEL_MAXLEN = 14


def _row_labels(meta: pd.DataFrame, rows: pd.Index) -> list[str]:
    labels = []
    for k in rows:
        r = meta.loc[k]
        s = str(r["cytoband"]) or f"{r['chrom']}:{int(r['start']) // 10**6}Mb"
        labels.append(s[:_LABEL_MAXLEN])
    if len(rows) > _MAX_LABEL_ROWS_FULL:
        # thin out: keep the first row of every chromosome run
        keep: set[int] = set()
        prev = None
        for i, k in enumerate(rows):
            c = meta.loc[k, "chrom"]
            if c != prev:
                keep.add(i)
                prev = c
        labels = [labels[i] if i in keep else "" for i in range(len(rows))]
    return labels


def _draw_panel(ax, ax_delta, matrix: pd.DataFrame, meta: pd.DataFrame,
                rows: pd.Index, blocks: list[dict], colors: list[str],
                cola: str, colb: str) -> None:
    sub = matrix.loc[rows]
    cmap = ListedColormap(colors)
    cmap.set_bad(_MISSING)
    ax.imshow(sub.to_numpy(dtype=float), aspect="auto", cmap=cmap,
              vmin=-0.5, vmax=2.5, interpolation="nearest")

    # cohort block separators (main heatmap only; the delta strip has its
    # own x-axis in Δfreq units and no sample columns)
    x0 = 0
    for bi, b in enumerate(blocks):
        n = len(b["columns"])
        if bi:
            ax.axvline(x0 - 0.5, color="white", lw=2.0)
        ax.axvline(x0 + n - 0.5, color="white", lw=2.0)
        x0 += n

    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(_row_labels(meta, rows), fontsize=5)
    ax.set_xticks([])
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)

    # delta-frequency side bars (cohort colours, not Amp/Del colours).
    # A cohort without a member peak (NaN) counts as frequency 0 — that is
    # what makes a cohort-specific peak show its full frequency difference.
    fa = np.nan_to_num(meta.loc[rows, "freq_a"].to_numpy(float))
    fb = np.nan_to_num(meta.loc[rows, "freq_b"].to_numpy(float))
    delta = fa - fb
    yy = np.arange(len(rows))
    ax_delta.barh(yy, delta, height=0.85,
                  color=[cola if d > 0 else colb for d in delta],
                  edgecolor="none")
    lim = max(0.05, float(np.nanmax(np.abs(delta))) * 1.15) if len(rows) else 0.05
    ax_delta.set_xlim(-lim, lim)
    ax_delta.set_ylim(len(rows) - 0.5, -0.5)
    ax_delta.axvline(0, color="#666666", lw=0.6)
    ax_delta.set_xticks([])
    ax_delta.set_yticks([])
    ax_delta.tick_params(length=0)
    for s in ax_delta.spines.values():
        s.set_visible(False)


def plot_oncoprint(onco: dict[str, Any], cfg: dict[str, Any]) -> None:
    """Draw the OncoPrint figure for one comparison."""
    matrix: pd.DataFrame = onco["matrix"]
    meta: pd.DataFrame = onco["meta"]
    blocks: list[dict] = onco["blocks"]
    if matrix is None or matrix.empty or meta.empty or len(blocks) < 2:
        return

    style(cfg)
    pal = theme.palette(cfg)

    comp = onco["comparison"]
    name = "P1_oncoprint" if not comp else f"P1_oncoprint_{comp}"
    la = blocks[0]["label"]
    lb = blocks[1]["label"]
    cola = pal.get(onco["a"], "#555555")
    colb = pal.get(onco["b"], "#777777")

    ptypes = [(t, m) for t, m in
              (("Amp", _AMP_COLORS), ("Del", _DEL_COLORS))
              if (meta["ptype"] == t).any()]
    if not ptypes:
        return

    fig_h = sum(min(6.0, max(1.2, 0.11 * int((meta["ptype"] == t).sum()) + 0.7))
                for t, _ in ptypes) + 1.15
    fig = plt.figure(figsize=(theme.capped_width(cfg, 7.0), fig_h))
    gs = fig.add_gridspec(
        len(ptypes), 2, width_ratios=[1, 0.13], wspace=0.04,
        height_ratios=[min(6.0, max(1.2, 0.11 * int((meta["ptype"] == t).sum())
                                    + 0.7)) for t, _ in ptypes],
        top=1.0 - 0.55 / fig_h, bottom=0.32 / fig_h,
        left=0.17, right=0.985)

    for gi, (ptype, colors) in enumerate(ptypes):
        rows = meta.index[meta["ptype"] == ptype]
        ax = fig.add_subplot(gs[gi, 0])
        ax_d = fig.add_subplot(gs[gi, 1])
        _draw_panel(ax, ax_d, matrix, meta, rows, blocks, colors, cola, colb)
        ax_d.set_xlabel("Δfreq", fontsize=6, labelpad=1)
        ax.annotate(ptype, xy=(-0.01, 1.02),
                    xycoords="axes fraction", ha="right", va="bottom",
                    fontsize=6.5, color=theme.ptype_color(ptype),
                    annotation_clip=False)
        if gi == 0:
            x0 = 0
            for b in blocks:
                n = len(b["columns"])
                ax.annotate(f"{b['label']}  (n={n})",
                            xy=(x0 + n / 2, 1.06),
                            xycoords=("data", "axes fraction"),
                            ha="center", va="bottom", fontsize=6.5,
                            fontweight="bold", annotation_clip=False)
                x0 += n

    fig.text(0.17, 0.06 / fig_h,
             f"{la} vs {lb} — 0/1/2 amplitude-threshold calls on union peaks "
             "(P1 moderate matching);  grey = no significant peak called in "
             "that cohort;  right bars: Δfreq (A − B; a cohort lacking the "
             "peak counts as freq 0) in cohort colours",
             fontsize=5.5, color="#555555")
    save(fig, cfg, name)


def figures_oncoprint(res: dict, cfg) -> None:
    """Entry point used by figures_p1; draws one figure per comparison."""
    for onco in res.get("onco", []) or []:
        try:
            plot_oncoprint(onco, cfg)
        except Exception:
            print(f"  [warn] oncoprint failed for {onco.get('comparison')}:")
            traceback.print_exc()
