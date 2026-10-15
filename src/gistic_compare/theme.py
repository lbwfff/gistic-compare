"""Publication figure theme, palette, layout and annotation helpers.

Every figure in the pipeline funnels through here so that the output shares a
single journal-ready look: dual-column width, Arial, embedded Type-42 fonts,
colour-blind-safe (Okabe-Ito) colours, no in-figure titles, and a consistent
way of reporting statistics (r / rho, 95% CI, n, p).
"""

from __future__ import annotations

from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

# ---------------------------------------------------------------------------
# geometry
# ---------------------------------------------------------------------------

COL_WIDTH = 3.5      # single column (in)
FULL_WIDTH = 7.0     # double column (in)
PANEL_H = 2.6        # default panel height (in)


def figure_cfg(cfg: dict[str, Any] | None) -> dict[str, Any]:
    out = (cfg or {}).get("output") or {}
    return out.get("figure") or {}


def width(cfg: dict[str, Any] | None = None) -> float:
    return float(figure_cfg(cfg).get("width_in", FULL_WIDTH))


def figsize(cfg, nrows: float = 1.0,
            panel_h: float = PANEL_H, total_w: float | None = None):
    """Figure size in inches: default full (double-column) width."""
    w = total_w if total_w is not None else width(cfg)
    return (w, panel_h * nrows)


def capped_width(cfg, wanted: float) -> float:
    """Clamp a requested figure width to the journal column limit."""
    return min(float(wanted), width(cfg))


def grid_shape(n: int, max_cols: int = 2) -> tuple[int, int]:
    """Wrap ``n`` panels into a grid of at most ``max_cols`` columns.

    Keeps multi-comparison figures inside the page width instead of growing
    linearly with the number of panels.
    """
    n = max(1, int(n))
    max_cols = max(1, int(max_cols))
    ncols = min(max_cols, n)
    nrows = int(np.ceil(n / ncols))
    return nrows, ncols


# ---------------------------------------------------------------------------
# theme
# ---------------------------------------------------------------------------

def apply_theme(cfg: dict[str, Any] | None = None) -> None:
    fcfg = figure_cfg(cfg)
    dpi = int(fcfg.get("dpi", 600))
    plt.rcParams.update({
        # fonts (Arial first; DejaVu fallback keeps things working headless)
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8,
        "axes.titlesize": 8,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 6.5,
        "legend.title_fontsize": 7,
        # export
        "figure.dpi": dpi,
        "savefig.dpi": dpi,
        "pdf.fonttype": 42,     # TrueType — editable, journal-compatible
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        # axes / layout
        "axes.linewidth": 0.8,
        "axes.grid": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": "#333333",
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "xtick.color": "#333333",
        "ytick.color": "#333333",
        "axes.labelcolor": "#333333",
        "text.color": "#333333",
        "lines.linewidth": 1.0,
        "legend.frameon": False,
        "legend.handlelength": 1.2,
        "legend.handletextpad": 0.4,
        "legend.borderaxespad": 0.3,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    })


def finalize(fig: plt.Figure, layout: str = "tight", **kwargs) -> None:
    """Apply a layout engine without ``bbox_inches='tight'`` cropping."""
    if layout == "constrained":
        try:
            fig.set_layout_engine("constrained", **kwargs)
            return
        except Exception:
            pass
    if layout in ("tight", "constrained"):
        try:
            fig.tight_layout()
        except Exception:
            pass


# ---------------------------------------------------------------------------
# annotations
# ---------------------------------------------------------------------------

def annot_stat(ax: plt.Axes, text: str, loc: str = "upper left",
               fontsize: float = 6.5, **kw) -> None:
    """Place a statistics string at a fixed corner of the axes.

    ``loc`` is one of the four corners, ``"top center"`` (inside the axes,
    centred horizontally — used on grouped distribution panels where a
    corner annotation would sit on one cohort's points), ``"upper left"`` or
    ``"lower right"`` shortcuts.
    """
    anchors = {
        "upper left": (0.03, 0.97, "top", "left"),
        "upper right": (0.97, 0.97, "top", "right"),
        "lower left": (0.03, 0.03, "bottom", "left"),
        "lower right": (0.97, 0.03, "bottom", "right"),
        "top center": (0.5, 0.97, "top", "center"),
    }
    x, y, va, ha = anchors[loc]
    ax.text(x, y, text, transform=ax.transAxes, ha=ha, va=va,
            fontsize=fontsize, **kw)


# ---------------------------------------------------------------------------
# colour
# ---------------------------------------------------------------------------

# Okabe & Ito (2008) colour-blind-safe qualitative palette
OKABE = {
    "black": "#000000",
    "orange": "#E69F00",
    "sky_blue": "#56B4E9",
    "bluish_green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "reddish_purple": "#CC79A7",
}

PTYPE_COLORS = {"Amp": OKABE["vermillion"], "Del": OKABE["blue"]}


def palette(cfg: dict[str, Any] | None = None) -> dict[str, str]:
    """Cohort colours: per-run palette from config (tab10 fallback applied
    in config at load time); empty when no config is given."""
    return dict((cfg or {}).get("_palette") or {})


def _rgb(h: str) -> tuple[float, float, float]:
    h = h.lstrip("#")
    return tuple(int(h[i:i + 2], 16) / 255.0 for i in (0, 2, 4))  # type: ignore


def colour_distance(a: str, b: str) -> float:
    ra, ga, ba = _rgb(a)
    rb, gb, bb = _rgb(b)
    return ((ra - rb) ** 2 + (ga - gb) ** 2 + (ba - bb) ** 2) ** 0.5


def categorical(n: int, taken: list[str] | None = None) -> list[str]:
    """``n`` mutually distinguishable categorical colours.

    Picks greedily by "farthest from everything already used" (including
    ``taken``) instead of taking tab10 in order, which can hand two
    categories near-identical colours.
    """
    from .config import TAB10

    taken = list(taken or [])
    out: list[str] = []
    pool = list(TAB10)
    for _ in range(max(0, int(n))):
        if not pool:
            pool = list(TAB10)
        used = taken + out
        if used:
            best = max(pool, key=lambda c: min(colour_distance(c, u)
                                               for u in used))
        else:
            best = pool[0]
        out.append(best)
        pool.remove(best)
    return out


def ptype_color(ptype: str) -> str:
    return PTYPE_COLORS.get(ptype, "#555555")


# ---------------------------------------------------------------------------
# comparison context: cohort ids -> labels / colours
# ---------------------------------------------------------------------------

def comp_sides(cfg: dict[str, Any] | None, comparison: str) -> tuple[str, str]:
    """Split ``"cohortA_vs_cohortB"`` into ``("cohortA", "cohortB")``."""
    text = str(comparison)
    if "_vs_" in text:
        a, b = text.split("_vs_", 1)
        return a, b
    return text, ""


def cohort_label(cfg: dict[str, Any] | None, cid: str) -> str:
    """Display label for a cohort id (falls back to the id itself)."""
    return ((cfg or {}).get("cohorts", {}).get(cid) or {}).get("label", cid)


def n_samples(cfg: dict[str, Any] | None, cid: str) -> int | None:
    """Sample count per cohort, recorded by the loader as ``_n_samples``."""
    v = ((cfg or {}).get("_n_samples") or {}).get(cid)
    try:
        return int(v) if v is not None else None
    except (TypeError, ValueError):
        return None


def side_meta(cfg: dict[str, Any] | None, comparison: str
              ) -> tuple[str, str, str, str]:
    """``(label_a, label_b, colour_a, colour_b)`` for a comparison id."""
    ca, cb = comp_sides(cfg, comparison)
    pal = palette(cfg)
    return (cohort_label(cfg, ca), cohort_label(cfg, cb),
            pal.get(ca, "#555555"), pal.get(cb, "#777777"))


def annot_n(ax: plt.Axes, text: str, loc: str = "lower right",
            fontsize: float = 6.0, color: str = "#555555") -> None:
    """Sample-size / unit annotation, visually quieter than the statistics."""
    if not text:
        return
    annot_stat(ax, text, loc=loc, fontsize=fontsize, color=color)


# ---------------------------------------------------------------------------
# panel letters
# ---------------------------------------------------------------------------

def panel_labels_enabled(cfg: dict[str, Any] | None) -> bool:
    return bool(figure_cfg(cfg).get("panel_labels", False))


def panel_label(ax: plt.Axes, letter: str,
                cfg: dict[str, Any] | None = None,
                dx: float = -0.075, dy: float = 1.03) -> None:
    """Bold panel letter (a, b, c ...); no-op unless enabled in the config."""
    if cfg is not None and not panel_labels_enabled(cfg):
        return
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=9,
            fontweight="bold", ha="right", va="bottom", color="#111111",
            clip_on=False)


PANEL_LETTERS = "abcdefghijklmnopqrstuvwxyz"


def panel_letter(i: int) -> str:
    return PANEL_LETTERS[i] if 0 <= i < len(PANEL_LETTERS) else str(i + 1)


# ---------------------------------------------------------------------------
# statistic formatting
# ---------------------------------------------------------------------------

def fmt_p(p: float | None) -> str:
    if p is None or not np.isfinite(p):
        return "p = NA"
    if p < 0.001:
        return "p < 0.001"
    return f"p = {p:.3g}"


def fmt_ci(symbol: str, r: float, lo: float | None = None,
           hi: float | None = None, n: int | None = None,
           p: float | None = None, digits: int = 2) -> str:
    """e.g. ``r = 0.49 (95% CI 0.10-0.75, n = 13), p = 0.086``."""
    r_txt = f"{symbol} = {r:.{digits}f}" if np.isfinite(r) else f"{symbol} = NA"
    parts = [r_txt]
    if lo is not None and hi is not None and np.isfinite(lo) and np.isfinite(hi):
        parts.append(f"(95% CI {lo:.{digits}f}\u2013{hi:.{digits}f})")
    if n is not None:
        parts.append(f"n = {int(n)}")
    s = " ".join(parts)
    if p is not None:
        s += f", {fmt_p(p)}"
    return s


def sig_label(p: float) -> str:
    if p is None or not np.isfinite(p):
        return "ns"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def bh_adjust(pvals) -> np.ndarray:
    """Benjamini-Hochberg FDR, NaN-safe and monotone."""
    p = np.asarray(pvals, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order]
    m = ranked.size
    q = ranked * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    inv = np.empty(m, dtype=int)
    inv[order] = np.arange(m)
    out[ok] = q[inv]
    return out


def wilson_ci(k, n, level: float = 0.95):
    """Vectorised Wilson score interval for binomial proportions.

    ``k`` successes out of ``n`` trials; returns ``(lo, hi)`` with NaN where
    ``n <= 0``. Preferred over the normal approximation for the small
    denominators and extreme proportions common in peak-level figures.
    """
    k = np.asarray(k, dtype=float)
    n = np.asarray(n, dtype=float)
    k, n = np.broadcast_arrays(k, n)
    lo = np.full(k.shape, np.nan)
    hi = np.full(k.shape, np.nan)
    ok = np.isfinite(k) & np.isfinite(n) & (n > 0)
    if not ok.any():
        return lo, hi
    from scipy import stats as _stats

    z = float(_stats.norm.isf((1.0 - level) / 2.0))
    kk = np.clip(k[ok], 0.0, n[ok])
    nn = n[ok]
    p = kk / nn
    denom = 1.0 + z * z / nn
    centre = (p + z * z / (2.0 * nn)) / denom
    half = (z * np.sqrt(p * (1.0 - p) / nn + z * z / (4.0 * nn * nn))
            / denom)
    lo[ok] = np.clip(centre - half, 0.0, 1.0)
    hi[ok] = np.clip(centre + half, 0.0, 1.0)
    return lo, hi


def annotate_proportion(ax: plt.Axes, k, n, column: str = "ratio",
                        loc: str = "lower right") -> None:
    """Append a ``ratio = k/n`` line to the corner annotation."""
    k = np.asarray(k, dtype=float)
    n = np.asarray(n, dtype=float)
    tot_k, tot_n = float(np.nansum(k)), float(np.nansum(n))
    if not np.isfinite(tot_n) or tot_n <= 0:
        return
    annot_n(ax, f"{column} = {tot_k:.0f}/{tot_n:.0f}", loc=loc)


def errorbar_from_ci(ax: plt.Axes, x, y, lo, hi, color: str = "#333333",
                     capsize: float = 1.8, lw: float = 0.8,
                     zorder: float = 5) -> None:
    """Symmetric error bars from pre-computed CI bounds (NaN-aware)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    lo = np.asarray(lo, dtype=float)
    hi = np.asarray(hi, dtype=float)
    ok = np.isfinite(y) & np.isfinite(lo) & np.isfinite(hi)
    if not ok.any():
        return
    yerr = np.vstack([np.clip(y[ok] - lo[ok], 0, None),
                      np.clip(hi[ok] - y[ok], 0, None)])
    ax.errorbar(x[ok], y[ok], yerr=yerr, fmt="none", ecolor=color,
                elinewidth=lw, capsize=capsize, capthick=lw, zorder=zorder)


def bar_with_ci(ax: plt.Axes, x, height, lo=None, hi=None, width: float = 0.6,
                color: str = "#555555", label: str | None = None,
                hatch: str | None = None, alpha: float = 1.0,
                edgecolor: str = "white", linewidth: float = 0.4,
                err_color: str = "#333333", **kw):
    """Bar chart with optional CI whiskers; returns the BarContainer."""
    x = np.asarray(x, dtype=float)
    height = np.asarray(height, dtype=float)
    bars = ax.bar(x, height, width, color=color, label=label, hatch=hatch,
                  alpha=alpha, edgecolor=edgecolor, linewidth=linewidth, **kw)
    if lo is not None and hi is not None:
        errorbar_from_ci(ax, x, height, lo, hi, color=err_color)
    return bars


def strip_overlay(ax: plt.Axes, data, xcol: str, ycol: str,
                  order=None, hue: str | None = None, palette=None,
                  dodge: bool = False, size: float = 2.2,
                  alpha: float = 0.55) -> None:
    """Jittered raw points over a violin/box.

    Small-n distributions drawn as bare violins hide how few samples there
    are; the overlay makes n and the spread directly visible.
    """
    import seaborn as sns

    if data is None or len(data) == 0 or ycol not in getattr(data, "columns", []):
        return
    sns.stripplot(data=data, x=xcol, y=ycol, order=order, hue=hue,
                  palette=palette, dodge=dodge, size=size, alpha=alpha,
                  jitter=0.14, legend=False, ax=ax, zorder=4)


def group_sizes(data, xcol: str, order=None, hue: str | None = None,
                hue_order=None) -> str:
    """``n = 13/13/13/13`` (per x group) or per x×hue cell."""
    if data is None or len(data) == 0 or xcol not in data.columns:
        return ""
    groups = list(order) if order is not None else \
        list(dict.fromkeys(data[xcol]))
    parts = []
    for g in groups:
        sub = data[data[xcol] == g]
        if hue is not None and hue in sub.columns:
            cells = (list(hue_order) if hue_order is not None
                     else list(dict.fromkeys(sub[hue])))
            parts.append("/".join(str(int((sub[hue] == c).sum()))
                                  for c in cells))
        else:
            parts.append(str(int(len(sub))))
    return "n = " + "/".join(parts)


def sig_bracket(ax: plt.Axes, x1: float, x2: float, y: float, p: float,
                h: float | None = None, color: str = "#333333") -> None:
    """Significance bracket between two x positions, labelled by stars."""
    lo, hi = (x1, x2) if x1 <= x2 else (x2, x1)
    span = ax.get_ylim()[1] - ax.get_ylim()[0]
    if h is None:
        h = 0.02 * span
    ax.plot([lo, lo, hi, hi], [y, y + h, y + h, y], lw=0.7, c=color,
            clip_on=False)
    ax.text((lo + hi) / 2, y + h, sig_label(p), ha="center", va="bottom",
            fontsize=6.5, color=color, clip_on=False)
