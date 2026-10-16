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
    """How the peak chains group into locus rows.

    The grouping runs on the **narrow peaks** (Peak Limits): two chains whose
    focal cores overlap count as one locus even though their wide-shaped
    flanks differ.  Rows without a narrow span fall back to their wide one.
    """
    recs = []
    for key, m in zip(meta.index, meta.to_dict("records")):
        rec = dict(m, peak_key=key)
        if "narrow_start" not in rec or pd.isna(rec.get("narrow_start")):
            rec["narrow_start"], rec["narrow_end"] = rec["start"], rec["end"]
        recs.append(rec)
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
            lo = max(a["narrow_start"], b["narrow_start"])
            hi = min(a["narrow_end"], b["narrow_end"])
            ov = max(0, hi - lo)
            shorter = max(1, min(a["narrow_end"] - a["narrow_start"],
                                 b["narrow_end"] - b["narrow_start"]))
            if ov / shorter < _MERGE_RECIP:
                continue
            used.add(bk)
            grp.append(bk)
        groups.append(grp)
    return groups


def _span(bp: float) -> str:
    """Genomic span in its most readable unit (kb below 1 Mb, else Mb)."""
    kb = abs(float(bp)) / 1000.0
    return f"{kb / 1000.0:.1f}Mb" if kb >= 1000.0 else f"{kb:.0f}kb"


def _locus_label(row: dict[str, Any]) -> str:
    """`cytoband (narrow-span / wide-span)` — both granularities on one row.

    A merged amp+del row reports the union of its members' narrow cores and
    the union of their wide peaks; when narrow == wide (a strictly focal
    peak that GISTIC could not extend) the single span is shown instead.
    """
    band = str(row.get("cytoband", "")).strip()
    if not band:
        band = f"{row['chrom']}:{int(row['start']) // 10**6}Mb"
    nar = abs(int(row.get("narrow_end", 0)) - int(row.get("narrow_start", 0)))
    wid = abs(int(row.get("end", 0)) - int(row.get("start", 0)))
    if abs(nar - wid) < 1000.0:
        text = _span(nar)          # narrow == wide: a cannot-extend focal peak
    else:
        text = f"{_span(nar)}/{_span(wid)}"
    return f"{band} ({text})"[:_LABEL_MAXLEN + 12]


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
        nar_s: list[int] = []
        nar_e: list[int] = []
        for k in keys:
            m = meta.loc[k]
            has_n = ("narrow_start" in m.index and "narrow_end" in m.index
                     and pd.notna(m["narrow_start"])
                     and pd.notna(m["narrow_end"]))
            nar_s.append(int(m["narrow_start"] if has_n else m["start"]))
            nar_e.append(int(m["narrow_end"] if has_n else m["end"]))
        first = meta.loc[keys[0]]
        row = dict(
            cell=vals,
            locus_key="|".join(keys),
            cytoband=str(first["cytoband"]),
            chrom=str(first["chrom"]),
            start=int(min(meta.loc[k, "start"] for k in keys)),
            end=int(max(meta.loc[k, "end"] for k in keys)),
            narrow_start=int(min(nar_s)), narrow_end=int(max(nar_e)),
            mixed=bool(amp_members and del_members))
        row["locus_label"] = _locus_label(row)
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


def _draw_top_strip(ax, cell: pd.DataFrame, blocks: list[dict],
                    vmax: float | None = None, rot90: bool = False) -> None:
    """Per-column stacked count of amp / del peaks (the TMB analogue).

    ``vmax`` pins the y scale so several strips can share it (the mirrored
    figure passes the max over both cohorts); ``rot90`` sets the label
    upright next to the ticks (v2) instead of the wide horizontal label the
    single-panel figure uses.
    """
    a = (cell == 1).sum(axis=0).to_numpy(float)
    d = (cell == 2).sum(axis=0).to_numpy(float)
    x = np.arange(len(a))
    ax.bar(x, a, width=1.0, color=_AMP, edgecolor="none")
    ax.bar(x, d, bottom=a, width=1.0, color=_DEL, edgecolor="none")
    tot = float((a + d).max()) if vmax is None else float(vmax)
    scale = max(3.0, tot * 1.30)
    ax.set_xlim(-0.5, len(a) - 0.5)
    ax.set_ylim(0, scale)
    hi = max(1, int(round(tot)))
    ax.set_yticks([0, hi])
    ax.set_yticklabels([0, hi], fontsize=5.5)
    if rot90:
        ax.set_ylabel("calls/sample", fontsize=6, rotation=90, labelpad=2)
    else:
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
              vmax=2.5, interpolation="nearest")    # light cell frames: a white cell then reads as an (empty) cell rather
    # than as blank canvas
    ax.set_axisbelow(False)
    ax.set_xticks(np.arange(-0.5, cell.shape[1] + 0.5, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, cell.shape[0] + 0.5, 1), minor=True)
    ax.grid(which="minor", color="#DEDEDE", lw=0.25, alpha=0.9)
    ax.tick_params(which="minor", length=0)
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


def _draw_right_bars(ax, frame: pd.DataFrame, blocks: list[dict]) -> None:
    """One panel, one stacked amp/del bar per cohort (upper = cohort A,
    lower = cohort B), with the combined frequency printed beside each bar.

    The cohort is encoded by **position** (upper/lower), never by colour:
    the vermillion/blue fills carry the Amp/Del direction, so no colour is
    doing double duty here.
    """
    yy = np.arange(len(frame))
    h = 0.30
    for bi, ycen in ((0, yy - 0.17), (1, yy + 0.17)):
        amp = frame[f"amp_{bi}"].to_numpy(float)
        dele = frame[f"del_{bi}"].to_numpy(float)
        ax.barh(ycen, amp, height=h, color=_AMP, edgecolor="white",
                linewidth=0.3)
        ax.barh(ycen, dele, left=amp, height=h, color=_DEL,
                edgecolor="white", linewidth=0.3)
        for y, f_ in zip(ycen, amp + dele):
            ax.text(f_ + 0.02, y, f"{f_:.2f}", fontsize=4.6, va="center",
                    ha="left", color="#454545")
    ax.set_xlim(0, 1.0)
    ax.set_yticks([])
    ax.set_ylim(len(frame) - 0.5, -0.5)
    ax.set_xticks([0, 0.5, 1.0])
    ax.set_xticklabels(["0", ".5", "1"], fontsize=5.5)
    ax.tick_params(length=2, pad=1)
    ax.set_xlabel("freq", fontsize=6, labelpad=1)
    ax.grid(axis="x", alpha=0.15)
    ax.legend([Patch(facecolor="#ffffff", edgecolor="#888888", linewidth=0.9)],
              [f"upper = {blocks[0]['label']}, "
               f"lower = {blocks[1]['label']}"],
              frameon=False, fontsize=5.6, loc="lower center",
              bbox_to_anchor=(0.5, 1.06), handlelength=1.2,
              handletextpad=0.4)
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
    labels = [str(v) for v in frame["locus_label"]]

    # waterfall: samples ordered by their alteration *pattern* over the
    # ranked loci, most-altered patterns first (left); ties fall back to the
    # number of calls
    for b in blocks:
        b["columns"] = _waterfall_order(cell, b["columns"], descending=True)
        b["_n"] = len(b["columns"])
    cell = cell.loc[:, [c for b in blocks for c in b["columns"]]]

    strip_h, header_h = 0.55, 0.46
    heat_h = 0.34 * len(frame) + 0.45
    fig_h = strip_h + header_h + heat_h
    fig = plt.figure(figsize=(theme.capped_width(cfg, 7.0), fig_h))
    gs = fig.add_gridspec(
        2, 2, width_ratios=[1, 0.30], wspace=0.06,
        height_ratios=[strip_h, heat_h],
        top=1.0 - header_h / fig_h, bottom=0.22 / fig_h,
        left=0.19, right=0.97)
    ax_top = fig.add_subplot(gs[0, 0])
    ax = fig.add_subplot(gs[1, 0])
    ax_b = fig.add_subplot(gs[1, 1])

    _draw_top_strip(ax_top, cell, blocks)
    _draw_heat(ax, cell, blocks)
    ax.set_yticks(range(len(frame)))
    ax.set_yticklabels(labels, fontsize=5)
    _draw_right_bars(ax_b, frame, blocks)

    x = 0
    for b in blocks:
        n = int(b["_n"])
        ax_top.annotate(f"{b['label']}  (n={n})",
                        xy=(x + n / 2, 1.12),
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


def _side_bars(ax, amp: np.ndarray, dele: np.ndarray, side: int,
               label: str) -> None:
    """One cohort's stacked amp/del frequencies growing *outward* from the
    centre (``side = -1`` leftwards, ``+1`` rightwards) with the combined
    frequency printed beyond the bar tip."""
    yy = np.arange(len(amp))
    h = 0.30
    if side > 0:
        ax.barh(yy, amp, height=h, color=_AMP, edgecolor="white", linewidth=0.3)
        ax.barh(yy, dele, left=amp, height=h, color=_DEL, edgecolor="white",
                linewidth=0.3)
        for y, f in zip(yy, amp + dele):
            ax.text(f + 0.03, y, f"{f:.2f}", fontsize=4.6, va="center",
                    ha="left", color="#454545")
        ax.set_xlim(0, 1.18)
        ax.set_xticks([0, 0.5, 1.0])
    else:
        ax.barh(yy, -amp, height=h, color=_AMP, edgecolor="white",
                linewidth=0.3)
        ax.barh(yy, -dele, left=-amp, height=h, color=_DEL, edgecolor="white",
                linewidth=0.3)
        for y, f in zip(yy, amp + dele):
            ax.text(-f - 0.03, y, f"{f:.2f}", fontsize=4.6, va="center",
                    ha="right", color="#454545")
        ax.set_xlim(-1.18, 0)
        ax.set_xticks([-1.0, -0.5, 0])
    ax.set_xticklabels(["1", ".5", "0"] if side < 0 else ["0", ".5", "1"],
                       fontsize=5.5)
    ax.set_yticks([])
    ax.set_ylim(len(amp) - 0.5, -0.5)
    ax.tick_params(length=2, pad=1)
    ax.set_xlabel(f"freq ({label})", fontsize=5.8, labelpad=1)
    ax.grid(axis="x", alpha=0.15)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)


def _waterfall_order(cell: pd.DataFrame, cols: list[str],
                     descending: bool) -> list[str]:
    """Samples ranked by their alteration pattern over the ranked rows.

    The key is the peak-by-peak pattern (peak 1's altered/unaltered split,
    then peak 2 within each group, …).  ``descending=True`` puts the
    most-altered patterns first (the v1 left-to-right cascade); the mirrored
    figure passes ``False`` for the cohort drawn left of the labels so its
    cascade falls away from the centre.
    """
    altered = cell.gt(0)
    calls = altered.sum(axis=0)
    pat = {c: tuple(bool(v) for v in altered[c]) for c in altered.columns}
    pool = [c for c in cols if c in altered.columns]

    def key(c: str):
        p = pat[c]
        if descending:                    # most altered first
            return (tuple(not v for v in p), -sum(p),
                    -float(calls.get(c, 0)), str(c))
        return (tuple(p), sum(p),         # least altered first (mirrored)
                float(calls.get(c, 0)), str(c))

    return sorted(pool, key=key)


def _draw_mirror_round(frame: pd.DataFrame, matrix: pd.DataFrame,
                       blocks: list[dict], cfg, name: str) -> None:
    """The v2 layout: freq | waterfall | labels | waterfall | freq.

    Left to right: cohort A's frequency bars, cohort A's sample waterfall,
    the peak labels, cohort B's waterfall, cohort B's frequency bars — the
    two waterfalls are the **same width** and mirror each other about the
    labels with their cascades falling *away* from the centre, and the two
    call strips share one y scale.
    """
    cell = pd.DataFrame(list(frame["cell"].values), columns=matrix.columns)
    n_samples = int(cell.shape[1])
    total = int((cell > 0).any(axis=0).sum())
    labels = [str(v) for v in frame["locus_label"]]
    # mirror: cohort A's cascade falls away from the labels (leftwards),
    # cohort B's away to the right -> most-altered samples hug the labels
    cols_a = _waterfall_order(cell, blocks[0]["columns"], descending=False)
    cols_b = _waterfall_order(cell, blocks[1]["columns"], descending=True)
    cell_a = cell.loc[:, cols_a]
    cell_b = cell.loc[:, cols_b]
    n_a, n_b = max(1, len(cols_a)), max(1, len(cols_b))

    strip_h, header_h, heat_h = 0.55, 0.34, 0.34 * len(frame) + 0.45
    fig_h = strip_h + header_h + heat_h
    lab_w, bar_w = 1.10, 0.36
    heat_w = 6.40 - lab_w - 2 * bar_w
    fig = plt.figure(figsize=(theme.capped_width(cfg, 7.0), fig_h))
    gs = fig.add_gridspec(
        2, 5, width_ratios=[bar_w, heat_w / 2.0, lab_w, heat_w / 2.0, bar_w],
        wspace=0.05, height_ratios=[strip_h, heat_h],
        top=1.0 - header_h / fig_h, bottom=0.22 / fig_h,
        left=0.085, right=0.985)
    ax_ba = fig.add_subplot(gs[1, 0])
    ax_ta = fig.add_subplot(gs[0, 1])
    ax_a = fig.add_subplot(gs[1, 1])
    ax_lab = fig.add_subplot(gs[1, 2])
    ax_b = fig.add_subplot(gs[1, 3])
    ax_tb = fig.add_subplot(gs[0, 3])
    ax_bb = fig.add_subplot(gs[1, 4])

    strip_max = max(float((cell_a > 0).sum(axis=0).max()),
                    float((cell_b > 0).sum(axis=0).max()))
    _draw_top_strip(ax_ta, cell_a, [dict(blocks[0], _n=n_a)],
                    vmax=strip_max, rot90=True)
    _draw_top_strip(ax_tb, cell_b, [dict(blocks[1], _n=n_b)],
                    vmax=strip_max, rot90=True)
    ax_tb.set_ylabel("")
    _draw_heat(ax_a, cell_a, [dict(blocks[0], _n=n_a)])
    _draw_heat(ax_b, cell_b, [dict(blocks[1], _n=n_b)])
    ax_a.set_yticks([])
    ax_b.set_yticks([])

    ax_lab.set_xticks([])
    ax_lab.set_yticks([])
    ax_lab.set_xlim(0, 1)
    ax_lab.set_ylim(len(frame) - 0.5, -0.5)
    for i, lab in enumerate(labels):
        # drawn centred inside the label column: tick labels would hang out
        # to the side and cover the neighbouring bar panel
        ax_lab.text(0.5, i, lab, transform=ax_lab.get_yaxis_transform(),
                    ha="center", va="center", fontsize=5, color="#333333")
    for s in ax_lab.spines.values():
        s.set_visible(False)

    _side_bars(ax_ba, frame["amp_0"].to_numpy(float),
               frame["del_0"].to_numpy(float), -1, blocks[0]["label"])
    _side_bars(ax_bb, frame["amp_1"].to_numpy(float),
               frame["del_1"].to_numpy(float), +1, blocks[1]["label"])

    ax_ta.annotate(f"{blocks[0]['label']}  (n={n_a})",
                   xy=(0.5, 1.05), xycoords="axes fraction", ha="center",
                   va="bottom", fontsize=6.5, fontweight="bold",
                   annotation_clip=False)
    ax_tb.annotate(f"{blocks[1]['label']}  (n={n_b})",
                   xy=(0.5, 1.05), xycoords="axes fraction", ha="center",
                   va="bottom", fontsize=6.5, fontweight="bold",
                   annotation_clip=False)
    pct = 100.0 * total / max(1, n_samples)
    fig.text(0.085, 1.0 - 0.08 / fig_h,
             f"Altered in {total} ({pct:.1f}%) of {n_samples} samples "
             f"(top-{len(frame)} peak loci)",
             fontsize=7.5, fontweight="bold", color="#111111",
             va="top", ha="left")
    save(fig, cfg, name)


def plot_oncoprint_v2(onco: dict[str, Any], cfg: dict[str, Any]) -> None:
    """Mirror-waterfall OncoPlot: freq | waterfall | labels | waterfall | freq.

    The two waterfalls are at equal width so the two cohorts are compared in
    the same space; each half's cascade falls away from the centre labels.
    """
    matrix: pd.DataFrame = onco["matrix"]
    meta: pd.DataFrame = onco["meta"]
    blocks: list[dict] = onco["blocks"]
    if matrix is None or matrix.empty or meta.empty or len(blocks) < 2:
        return
    style(cfg)
    comp = onco["comparison"]
    name = "P1_oncoplot_v2" if not comp else f"P1_oncoplot_v2_{comp}"
    blocks = [dict(b) for b in blocks]
    for b in blocks:
        b["_n"] = len(b["columns"])
    frame = _locus_table(meta, matrix, blocks)
    if frame.empty:
        return
    frame = frame.head(_TOP_ROWS).reset_index(drop=True)
    _draw_mirror_round(frame, matrix, blocks, cfg, name)


def figures_oncoplot(res: dict, cfg) -> None:
    """Entry point used by figures_p1; draws one figure per comparison."""
    for onco in res.get("onco", []) or []:
        try:
            plot_oncoplot(onco, cfg)
        except Exception:
            print(f"  [warn] oncoplot failed for {onco.get('comparison')}:")
            traceback.print_exc()
        try:
            plot_oncoprint_v2(onco, cfg)
        except Exception:
            print(f"  [warn] oncoplot v2 failed for "
                  f"{onco.get('comparison')}:")
            traceback.print_exc()
