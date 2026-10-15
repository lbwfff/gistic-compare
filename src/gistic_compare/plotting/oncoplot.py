"""Union-peak OncoPlot: per-sample GISTIC calls across cohorts.

Rows are **union peak loci**: the P1 moderate-match chains are grouped so an
Amplification chain and a Deletion chain whose intervals overlap (reciprocal
>= 50 % of the shorter) share a single row — that way one peak can carry both
the Amp and the Del alteration state, which is what the right-hand stacked
frequencies need.  Cells show the alteration state (vermillion = amplified,
blue = deleted) or white (no alteration, *or* the cohort never called a peak
there — the two are deliberately not distinguished here).

Rows are ranked by combined alteration frequency (top N only); the two
cohorts' sample columns sit side by side in one heat, split by a white
separator; the top strip counts each column's amp/del calls and the right
panel carries the per-cohort amp/del stacked frequencies.

Layout is fully manual (no tight_layout) because cohort labels, the
"altered in ..." header and the footnote live outside the axes.
"""

from __future__ import annotations

import traceback
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

from .. import theme
from ._core import save, style

# cell states: white = no alteration / no peak called there
_AMP = theme.OKABE["vermillion"]
_DEL = theme.OKABE["blue"]
_STATE_COLORS = ["#ffffff", _AMP, _DEL]

_TOP_ROWS = 10        # rows shown, ranked by combined alteration frequency
_MERGE_RECIP = 0.5    # amp/del chains sharing >= this of the shorter -> one locus
_LABEL_MAXLEN = 18


def _merge_loci(meta: pd.DataFrame) -> list[list[str]]:
    """How the peak chains group into locus rows."""
    recs = meta.reset_index().to_dict("records")
    used: set[str] = set()
    groups: list[list[str]] = []
    for i, a in enumerate(recs):
        key = a["peak_key"]
        if key in used:
            continue
        used.add(key)
        grp = [key]
        for b in recs[i + 1:]:
            bk = b["peak_key"]
            if bk in used or b["chrom"] != a["chrom"]:
                continue
            lo = max(a["start"], b["start"])
            hi = min(a["end"], b["end"])
            ov = max(0, hi - lo)
            shorter = max(1, min(a["end"] - a["start"],
                                 b["end"] - b["start"]))
            if ov / shorter < _MERGE_RECIP:
                continue
            used.add(bk)
            grp.append(bk)
        groups.append(grp)
    return groups


def _locus_label(meta: pd.DataFrame, keys: list[str]) -> str:
    band = str(meta.loc[keys[0], "cytoband"]).strip()
    if not band:
        band = (f"{meta.loc[keys[0], 'chrom']}:"
                f"{int(meta.loc[keys[0], 'start']) // 10**6}Mb")
    bands = {str(meta.loc[k, "cytoband"]).strip() for k in keys}
    label = band if len(bands) == 1 else f"{band}+"
    return label[:_LABEL_MAXLEN]


def _locus_table(meta: pd.DataFrame, matrix: pd.DataFrame,
                 blocks: list[dict]) -> pd.DataFrame:
    """Per-locus cell values and per-cohort amp/del frequencies."""
    rows: list[dict[str, Any]] = []
    for keys in _merge_loci(meta):
        keys = [k for k in keys if k in matrix.index]
        if not keys:
            continue
        sub = matrix.loc[keys]
        amp_members = [k for k in keys if str(meta.loc[k, "ptype"]) == "Amp"]
        del_members = [k for k in keys if str(meta.loc[k, "ptype"]) == "Del"]
        if not amp_members and not del_members:
            continue
        amps = ((sub.loc[amp_members] > 0).any(axis=0) if amp_members
                else pd.Series(False, index=sub.columns))
        dels = ((sub.loc[del_members] > 0).any(axis=0) if del_members
                else pd.Series(False, index=sub.columns))
        finite = sub.notna().any(axis=0)
        vals = pd.Series(np.nan, index=sub.columns, dtype=float)
        vals[finite] = 0.0
        vals[amps] = 1.0
        vals[dels & ~amps] = 2.0

        counts: dict[str, float] = {}
        for bi, b in enumerate(blocks):
            cols = [c for c in b["columns"] if c in vals.index]
            n = len(cols) or 1
            v = vals.loc[cols]
            na_, nd_ = int((v == 1).sum()), int((v == 2).sum())
            counts[f"n_amp_{bi}"] = float(na_)
            counts[f"n_del_{bi}"] = float(nd_)
            counts[f"n_{bi}"] = float(n)
            counts[f"amp_{bi}"] = na_ / n
            counts[f"del_{bi}"] = nd_ / n
        denom = counts["n_0"] + counts["n_1"]
        first = meta.loc[keys[0]]
        row = dict(
            cell=vals,
            locus_key="|".join(keys),
            cytoband=_locus_label(meta, keys),
            chrom=str(first["chrom"]), start=int(first["start"]),
            end=int(max(meta.loc[k, "end"] for k in keys)),
            mixed=bool(amp_members and del_members))
        row.update(counts)
        row["freq_union"] = float(
            (counts["n_amp_0"] + counts["n_del_0"]
             + counts["n_amp_1"] + counts["n_del_1"]) / denom)
        rows.append(row)
    if not rows:
        return pd.DataFrame()
    frame = pd.DataFrame(rows)
    frame = frame.sort_values(
        ["freq_union", "chrom", "start"], ascending=[False, True, True])
    return frame.reset_index(drop=True)


def _draw_top_strip(ax, cell: pd.DataFrame, blocks: list[dict]) -> None:
    """Per-column stacked count of amp / del peaks (the TMB analogue)."""
    a = (cell == 1).sum(axis=0).to_numpy(float)
    d = (cell == 2).sum(axis=0).to_numpy(float)
    x = np.arange(len(a))
    ax.bar(x, a, width=1.0, color=_AMP, edgecolor="none")
    ax.bar(x, d, bottom=a, width=1.0, color=_DEL, edgecolor="none")
    tot = float((a + d).max())
    ax.set_xlim(-0.5, len(a) - 0.5)
    ax.set_ylim(0, max(3.0, tot * 1.30))
    hi = max(1, int(round(tot)))
    ax.set_yticks([0, hi])
    ax.set_yticklabels([0, hi], fontsize=5.5)
    ax.set_ylabel("calls/sample", fontsize=5.5, rotation=0, ha="right",
                  va="center", labelpad=32)
    ax.grid(axis="y", alpha=0.15)
    ax.set_xticks([])
    ax.minorticks_off()
    x0 = 0
    for bi, b in enumerate(blocks):
        n = int(b["_n"])
        if bi:
            ax.axvline(x0 - 0.5, color="white", lw=2.0)
        ax.axvline(x0 + n - 0.5, color="white", lw=2.0)
        x0 += n
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def _draw_heat(ax, cell: pd.DataFrame, blocks: list[dict]):
    ax.imshow(cell.to_numpy(dtype=float),
              aspect="auto", cmap=ListedColormap(_STATE_COLORS), vmin=-0.5,
              vmax=2.5, interpolation="nearest")
    x = 0
    for bi, b in enumerate(blocks):
        n = int(b["_n"])
        if bi:
            ax.axvline(x - 0.5, color="white", lw=2.0)
        ax.axvline(x + n - 0.5, color="white", lw=2.0)
        x += n
    ax.set_xticks([])
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)


def _draw_right_bars(ax, frame: pd.DataFrame, blocks: list[dict],
                     cola: str, colb: str) -> None:
    """Per-locus stacked amp/del frequencies, one thin bar per cohort."""
    yy = np.arange(len(frame))
    h = 0.30
    for bi, cohort_c in ((0, cola), (1, colb)):
        ycen = yy - (0.17 if bi == 0 else -0.17)
        amp = frame[f"amp_{bi}"].to_numpy(float)
        dele = frame[f"del_{bi}"].to_numpy(float)
        ax.barh(ycen, amp, height=h, color=_AMP, edgecolor=cohort_c,
                linewidth=0.7)
        ax.barh(ycen, dele, left=amp, height=h, color=_DEL,
                edgecolor=cohort_c, linewidth=0.7)
        for y, f_ in zip(ycen, amp + dele):
            if f_ <= 0.002:
                continue
            ax.text(f_ + 0.02, y, f"{f_:.2f}", fontsize=4.6, va="center",
                    ha="left", color="#454545")
    ax.set_xlim(0, 1.0)
    ax.set_yticks([])
    ax.set_ylim(len(frame) - 0.5, -0.5)
    ax.set_xticks([0, 0.5, 1.0])
    ax.set_xticklabels(["0", ".5", "1"], fontsize=5.5)
    ax.tick_params(length=2, pad=1)
    ax.set_xlabel("freq", fontsize=6, labelpad=1)
    ax.legend([Patch(facecolor="#ffffff", edgecolor=cola, linewidth=1.2),
               Patch(facecolor="#ffffff", edgecolor=colb, linewidth=1.2)],
              [blocks[0]["label"], blocks[1]["label"]], frameon=False,
              fontsize=6, loc="lower center", bbox_to_anchor=(0.5, 1.06),
              handlelength=1.4, handletextpad=0.3, columnspacing=0.6, ncol=1,
              title="cohort edge", title_fontsize=6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def plot_oncoplot(onco: dict[str, Any], cfg: dict[str, Any]) -> None:
    """Draw the OncoPlot figure for one comparison."""
    matrix: pd.DataFrame = onco["matrix"]
    meta: pd.DataFrame = onco["meta"]
    blocks: list[dict] = onco["blocks"]
    if matrix is None or matrix.empty or meta.empty or len(blocks) < 2:
        return

    style(cfg)
    pal = theme.palette(cfg)
    comp = onco["comparison"]
    name = "P1_oncoplot" if not comp else f"P1_oncoplot_{comp}"
    la = blocks[0]["label"]
    lb = blocks[1]["label"]
    cola = pal.get(onco["a"], "#555555")
    colb = pal.get(onco["b"], "#777777")
    blocks = [dict(b) for b in blocks]      # local copy: do not mutate input
    for b in blocks:
        b["_n"] = len(b["columns"])

    frame = _locus_table(meta, matrix, blocks)
    if frame.empty:
        return
    frame = frame.head(_TOP_ROWS).reset_index(drop=True)
    cell = pd.DataFrame(list(frame["cell"].values), columns=matrix.columns)
    n_samples = int(cell.shape[1])
    total = int((cell > 0).any(axis=0).sum())
    labels = []
    for _, r in frame.iterrows():
        base = str(r["cytoband"])
        if r["mixed"]:
            base += \
                f"  ({r['chrom']}:{int(r['start']) // 10**6}Mb)"
        labels.append(base)

    # waterfall: inside each cohort block the samples *without* the top
    # locus come first (left) and the altered ones last (right); ties are
    # ordered by total calls descending so the shapes fall off smoothly
    calls_per_sample = cell.gt(0).sum(axis=0)
    top_hit = cell.iloc[0].gt(0)
    for b in blocks:
        cols = [c for c in b["columns"] if c in cell.columns]
        cols = sorted(cols, key=lambda c: (
            not bool(top_hit.get(c, False)), -float(calls_per_sample.get(c, 0))))
        b["columns"] = cols
        b["_n"] = len(cols)
    cell = cell.loc[:, [c for b in blocks for c in b["columns"]]]

    strip_h, header_h, note_h = 0.85, 0.46, 0.0
    heat_h = 0.34 * len(frame) + 0.45
    fig_h = strip_h + header_h + note_h + heat_h
    fig = plt.figure(figsize=(theme.capped_width(cfg, 7.0), fig_h))
    gs = fig.add_gridspec(
        2, 2, width_ratios=[1, 0.36], wspace=0.05,
        height_ratios=[strip_h, heat_h],
        top=1.0 - header_h / fig_h, bottom=note_h / fig_h,
        left=0.19, right=0.97)
    ax_top = fig.add_subplot(gs[0, 0])
    ax = fig.add_subplot(gs[1, 0])
    ax_b = fig.add_subplot(gs[1, 1])

    _draw_top_strip(ax_top, cell, blocks)
    _draw_heat(ax, cell, blocks)
    ax.set_yticks(range(len(frame)))
    ax.set_yticklabels(labels, fontsize=5)
    _draw_right_bars(ax_b, frame, blocks, cola, colb)

    x = 0
    for b in blocks:
        n = int(b["_n"])
        ax_top.annotate(f"{b['label']}  (n={n})",
                        xy=(x + n / 2, 1.04),
                        xycoords=("data", "axes fraction"),
                        ha="center", va="bottom", fontsize=6.5,
                        fontweight="bold", annotation_clip=False)
        x += n
    pct = 100.0 * total / max(1, n_samples)
    fig.text(0.19, 1.0 - 0.08 / fig_h,
             f"Altered in {total} ({pct:.1f}%) of {n_samples} samples "
             f"(top-{len(frame)} peak loci by combined frequency)",
             fontsize=7.5, fontweight="bold", color="#111111",
             va="top", ha="left")
    save(fig, cfg, name)


def figures_oncoplot(res: dict, cfg) -> None:
    """Entry point used by figures_p1; draws one figure per comparison."""
    for onco in res.get("onco", []) or []:
        try:
            plot_oncoplot(onco, cfg)
        except Exception:
            print(f"  [warn] oncoplot failed for {onco.get('comparison')}:")
            traceback.print_exc()
