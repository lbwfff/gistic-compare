"""Shared plotting style and per-module figures.

All figures follow the journal-ready conventions centralized in
:mod:`gistic_compare.theme` (dual-column width, Arial, Type-42 fonts,
Okabe-Ito colours, no in-figure titles). Descriptive text belongs in the
figure caption; only statistics/identifiers stay on the plot.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.ticker import FuncFormatter, MaxNLocator, NullFormatter
from scipy import stats

from .. import theme
from ..config import out_fig


def _palette(cfg: dict[str, Any]) -> dict[str, str]:
    """Cohort colours (config-driven Okabe-Ito; tab10 fallback)."""
    return theme.palette(cfg)


def _pc(ptype: str) -> str:
    return theme.ptype_color(ptype)


def _sides(cfg: dict[str, Any], comparison: str):
    """(label_a, label_b, colour_a, colour_b) for a comparison id."""
    return theme.side_meta(cfg, comparison)


def _kind(cfg: dict[str, Any], comparison: str) -> str:
    """``paired`` / ``unpaired`` for a comparison id ('' when unknown)."""
    for c in cfg.get("comparisons", []) or []:
        if f"{c.get('a')}_vs_{c.get('b')}" == str(comparison):
            return str(c.get("kind", "") or "")
    return ""


def _kind_note(cfg: dict[str, Any], comparison: str) -> str:
    k = _kind(cfg, comparison)
    return f"{k} comparison" if k else ""


def _n(cfg: dict[str, Any], cid: str) -> int | None:
    """Samples in a cohort, as recorded by the loader (None if unknown)."""
    return theme.n_samples(cfg, cid)


def _n_pair(cfg: dict[str, Any], comparison: str) -> str:
    """``n = 136 vs 90`` sample-size string for a comparison."""
    ca, cb = theme.comp_sides(cfg, comparison)
    na, nb = _n(cfg, ca), _n(cfg, cb)
    if na is None or nb is None:
        return ""
    return f"n = {na} vs {nb}"



def style(cfg: dict[str, Any] | None = None) -> None:
    theme.apply_theme(cfg)


def save(fig: plt.Figure, cfg: dict[str, Any], name: str) -> Path:
    for ext in ("pdf", "png"):
        p = out_fig(cfg, f"{name}.{ext}")
        fig.savefig(p)
    plt.close(fig)
    return out_fig(cfg, f"{name}.png")




# ---------------------------------------------------------------------------
# shared genome-track helpers (moved from the P1 section)
# ---------------------------------------------------------------------------

def _genome_x(df: pd.DataFrame) -> np.ndarray:
    """Linear genome x-coordinate from chrom + midpoint."""
    from ..config import CHROMS, HG19_CHROM_SIZES
    offsets: dict[str, int] = {}
    acc = 0
    for c in CHROMS:
        offsets[c] = acc
        acc += HG19_CHROM_SIZES[c]
    mid = (df["start"].astype(float) + df["end"].astype(float)) / 2.0
    return df["chrom"].map(offsets).to_numpy(float) + mid.to_numpy(float)

def _chrom_bounds(xdf: pd.DataFrame) -> list[tuple[float, str]]:
    from ..config import CHROMS, HG19_CHROM_SIZES
    bounds = []
    acc = 0.0
    for c in CHROMS:
        bounds.append((acc + HG19_CHROM_SIZES[c] / 2.0, c))
        acc += HG19_CHROM_SIZES[c]
    return bounds

def _draw_chrom_guides(ax, xdf: pd.DataFrame) -> None:
    """Alternating chromosome bands drawn at their true genome offsets.

    The offsets must match :func:`_genome_x`, i.e. accumulate over *all*
    autosomes. Deriving them from only the observed chromosomes silently
    shifts every band and its dividers whenever chr1 has no data point.
    """
    from ..config import CHROMS, HG19_CHROM_SIZES
    acc = 0.0
    for i, c in enumerate(CHROMS):
        size = HG19_CHROM_SIZES.get(c)
        if size is None:
            continue
        if i % 2 == 1:
            ax.axvspan(acc, acc + size, color="#F4F4F4", zorder=0)
        if i > 0:
            ax.axvline(acc, color="#CCCCCC", lw=0.5, zorder=1)
        acc += size

def _clean_label(v) -> str:
    """String form of a label cell, with pandas/NumPy missing values → ''."""
    if v is None:
        return ""
    try:
        if isinstance(v, float) and not np.isfinite(v):
            return ""
    except TypeError:
        pass
    s = str(v).strip()
    return "" if s.lower() in ("nan", "none", "<na>", "nat") else s

def _label_series(df: pd.DataFrame, cols: list[str],
                  fallback: str | None = None) -> list[str]:
    """First non-empty label among ``cols``, else the ``fallback`` column.

    Without the fallback a missing cytoband yields "" and the point silently
    loses its label instead of degrading to the peak id.
    """
    out = []
    for _, r in df.iterrows():
        lab = ""
        for c in cols:
            lab = _clean_label(r.get(c))
            if lab:
                break
        if not lab and fallback:
            lab = _clean_label(r.get(fallback))
        out.append(lab)
    return out

def _annotate_labeled_points(ax, xs, ys, labels,
                             fontsize: float = 4.8,
                             min_sep_pt: float = 16.0,
                             offsets: tuple[tuple[float, float], ...] = (
                                 (4, 4), (4, -8), (-4, 4), (-4, -8),
                                 (8, 0), (-8, 0), (0, 10), (0, -12))) -> None:
    """Place text labels near points with simple screen-space collision skip.

    Labels that would overlap an already-placed one are dropped, which keeps
    dense scatters readable. ``min_sep_pt`` is a *physical* separation in
    points: it is converted to display pixels with the live dpi so that the
    guard does not silently shrink (and stop working) at 600 dpi export.
    """
    if len(xs) == 0:
        return
    fig = ax.figure
    fig.canvas.draw()
    inv = ax.transData.inverted()
    px_per_pt = float(fig.dpi) / 72.0
    min_sep_px = min_sep_pt * px_per_pt
    placed: list[tuple[float, float]] = []
    # draw high-signal points first (larger y) so they win label space
    order = np.argsort(-np.asarray(ys, dtype=float))
    for oi in order:
        x, y = float(np.asarray(xs)[oi]), float(np.asarray(ys)[oi])
        lab = str(np.asarray(labels)[oi]).strip()
        if not lab or lab.lower() in ("nan", "none", "<na>", "nat"):
            continue
        if not (np.isfinite(x) and np.isfinite(y)):
            continue
        best = None
        for dx, dy in offsets:
            # offset in points → data via display transform
            disp = ax.transData.transform((x, y))
            txt_disp = (disp[0] + dx * px_per_pt, disp[1] + dy * px_per_pt)
            txt_data = inv.transform(txt_disp)
            # rough text extent in display px (fontsize pt → px)
            w_px = len(lab) * fontsize * 0.55 * px_per_pt
            corners = [
                txt_disp,
                (txt_disp[0] + w_px, txt_disp[1] + fontsize * px_per_pt),
            ]
            ok = True
            for cx, cy in corners:
                for px, py in placed:
                    if abs(cx - px) < min_sep_px and abs(cy - py) < min_sep_px:
                        ok = False
                        break
                if not ok:
                    break
            if ok:
                best = (txt_data, txt_disp, w_px)
                break
        if best is None:
            continue
        txt_data, txt_disp, w_px = best
        ax.annotate(
            lab, xy=(x, y), xytext=(txt_data[0], txt_data[1]),
            textcoords="data", fontsize=fontsize, color="#333333",
            ha="left" if txt_disp[0] >= ax.transData.transform((x, y))[0] else "right",
            va="bottom" if txt_disp[1] >= ax.transData.transform((x, y))[1] else "top",
            zorder=4,
        )
        placed.append(txt_disp)
        placed.append((txt_disp[0] + w_px, txt_disp[1]))


__all__ = [
    "Any",
    "FuncFormatter",
    "MaxNLocator",
    "NullFormatter",
    "Path",
    "_annotate_labeled_points",
    "_chrom_bounds",
    "_clean_label",
    "_draw_chrom_guides",
    "_genome_x",
    "_kind",
    "_kind_note",
    "_label_series",
    "_n",
    "_n_pair",
    "_palette",
    "_pc",
    "_sides",
    "matplotlib",
    "np",
    "out_fig",
    "pd",
    "plt",
    "save",
    "sns",
    "stats",
    "style",
    "theme",
]
