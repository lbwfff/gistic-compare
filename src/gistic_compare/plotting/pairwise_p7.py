"""Figures for module P7 (amplitude differences among altered samples)."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ._core import *  # noqa: F401,F403

_CLASS_COLORS = {
    "both": "#111111",
    "freq_driven": "#E69F00",     # okabe orange
    "amp_driven": "#56B4E9",      # okabe sky blue
    "neither": "#CFCFCF",
}


def plot_p7_quadrant(df: pd.DataFrame, cfg, unit: str, prefix: str = "P7"
                     ) -> None:
    """Delta-frequency vs delta-amplitude, one panel per direction.

    Every point is one genomic unit (window/band/peak): x = how much more
    often it is altered in cohort A, y = how much stronger the alteration
    is among its altered samples.  Point size scales with frequency
    significance; text labels go only to units that are significant on at
    least one axis (strongest combined FDR first) — grey "neither" points
    are deliberately never labelled.
    """
    if df is None or df.empty:
        return
    style(cfg)
    fig, axes = plt.subplots(1, 2, figsize=theme.figsize(
        cfg, panel_h=2.9, total_w=7.0))
    for ax, direction in zip(axes, ("amp", "del")):
        d = df[np.isfinite(df[f"delta_amp_{direction}"])].copy()
        if d.empty:
            ax.set_visible(False)
            continue
        for cls, color in _CLASS_COLORS.items():
            sub = d[d[f"class_{direction}"] == cls]
            if sub.empty:
                continue
            size = -np.log10(np.clip(
                sub[f"fdr_freq_{direction}"].to_numpy(float), 1e-300, 1.0))
            ax.scatter(sub[f"delta_freq_{direction}"],
                       sub[f"delta_amp_{direction}"],
                       s=4 + 5 * size, color=color,
                       alpha=0.75 if cls != "neither" else 0.45,
                       edgecolor="none", label=cls, zorder=3)
        ax.axhline(0, color="#777777", lw=0.6, ls=":")
        ax.axvline(0, color="#777777", lw=0.6, ls=":")
        # label only units significant on at least one axis, strongest
        # combined FDR first; "neither" points stay unlabelled
        comb = np.minimum(d[f"fdr_freq_{direction}"].fillna(1.0),
                          d[f"fdr_amp_{direction}"].fillna(1.0))
        labels = np.where(d[f"class_{direction}"] != "neither",
                          d["unit"].astype(str).to_numpy(), "")
        preset_rings_and_top_labels(
            ax, d[f"delta_freq_{direction}"].to_numpy(dtype=float),
            d[f"delta_amp_{direction}"].to_numpy(dtype=float),
            labels, cfg,
            order_key=-comb.to_numpy(dtype=float), fontsize=5)
        la, lb, _, _ = _sides(cfg, d["comparison"].iloc[0])
        ax.set_xlabel(f"\u0394freq [{la} \u2212 {lb}]")
        ax.set_ylabel(f"\u0394amplitude among altered\n[{la} \u2212 {lb}]")
        ax.set_title(f"{direction.upper()}  (n={len(d)})", fontsize=7)
        ax.grid(alpha=0.2)
        ax.legend(fontsize=5.5, loc="lower left", handletextpad=0.2)
        theme.annot_n(ax, f"point = one {unit}; "
                          "size = -log10 FDR(freq)",
                      loc="upper right", fontsize=5.0)
    fig.suptitle(f"{df['comparison'].iloc[0]}: frequency \u00d7 amplitude "
                 f"by {unit}", fontsize=8)
    theme.finalize(fig, "tight")
    save(fig, cfg, f"{prefix}_{unit}_quadrant")


def plot_p7_volcano(df: pd.DataFrame, cfg, unit: str, prefix: str = "P7"
                    ) -> None:
    """Amplitude difference vs -log10 FDR (altered-only MW), per direction."""
    if df is None or df.empty:
        return
    style(cfg)
    fig, axes = plt.subplots(1, 2, figsize=theme.figsize(
        cfg, panel_h=2.7, total_w=7.0))
    for ax, direction in zip(axes, ("Amp", "Del")):
        d = direction.lower()
        fcol = f"fdr_amp_{d}"
        s = df[np.isfinite(df[fcol])].copy()
        if s.empty:
            ax.set_visible(False)
            continue
        s["y"] = -np.log10(np.clip(s[fcol].astype(float), 1e-300, 1.0))
        sig = s[s[fcol] < 0.05]
        col = theme.ptype_color(direction)
        ax.scatter(s[f"delta_amp_{d}"], s["y"], s=6, color="#CFCFCF",
                   edgecolor="none")
        ax.scatter(sig[f"delta_amp_{d}"], sig["y"], s=9, color=col,
                   edgecolor="none")
        ax.axvline(0, color="#777777", lw=0.6, ls=":")
        info = preset_rings_and_top_labels(
            ax, s[f"delta_amp_{d}"].to_numpy(dtype=float),
            s["y"].to_numpy(dtype=float),
            s["unit"].astype(str).to_numpy(), cfg,
            order_key=s["y"].to_numpy(dtype=float), fontsize=5)
        la, lb, _, _ = _sides(cfg, s["comparison"].iloc[0])
        ax.set_xlabel(f"\u0394 amplitude among altered [{la} \u2212 {lb}]")
        ax.set_ylabel("\u2212log10 FDR")
        ax.set_title(f"{direction}  (sig n={len(sig)})", fontsize=7)
        ax.grid(alpha=0.2)
    fig.suptitle(f"{df['comparison'].iloc[0]}: amplitude differential "
                 f"by {unit}", fontsize=8)
    if info["ringed"]:
        theme.annot_n(ax, "\u25cb = pre-specified loci/bands",
                      loc="upper left")
    theme.finalize(fig, "tight")
    save(fig, cfg, f"{prefix}_{unit}_volcano")


def _example_selection(df: pd.DataFrame, n_each: int = 3) -> pd.DataFrame:
    """Pick example units: the strongest by amplitude FDR and by frequency
    FDR (top ``n_each`` of each axis, de-duplicated, amplitude first).

    Ranking by a single combined score would fill the panel with adjacent
    windows of one story (e.g. a whole frequency-different chromosome
    arm); taking the best of *each axis* keeps both narratives visible.
    """
    d = df.copy()
    rank_freq = np.minimum(d["fdr_freq_amp"].fillna(1.0),
                           d["fdr_freq_del"].fillna(1.0))
    rank_amp = np.minimum(d["fdr_amp_amp"].fillna(1.0),
                          d["fdr_amp_del"].fillna(1.0))
    d["_rank_freq"], d["_rank_amp"] = rank_freq, rank_amp
    seen: list[str] = []
    rows = []
    for rank_col in ("_rank_amp", "_rank_freq"):
        take = d.sort_values(rank_col)
        for _, r in take.iterrows():
            u = str(r["unit"])
            if u in seen:
                continue
            if len(seen) >= 2 * n_each:
                break
            seen.append(u)
            rows.append(r)
        if len(seen) >= 2 * n_each:
            break
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    # amplitude stories first (the module's unique contribution)
    return out.sort_values("_rank_amp").reset_index(drop=True)


def plot_p7_examples(res: dict, cfg, unit: str = "window",
                     prefix: str = "P7", n_each: int = 3) -> None:
    """Per-sample copy value for the strongest units.

    Selection: the top units by amplitude FDR and by frequency FDR (up to
    three each).  Each point is one sample plotted at its cohort's fixed
    position (jittered); dashed/dotted lines mark GISTIC's alteration
    threshold (±0.1) and the high-level threshold (±0.9).
    """
    df = res.get(unit)
    if df is None or df.empty:
        return
    style(cfg)
    d = _example_selection(df, n_each=n_each)
    if d.empty:
        return
    comps = list(dict.fromkeys(d["comparison"]))
    if len(comps) != 1:
        return
    comp = comps[0]
    ca, cb = theme.comp_sides(cfg, comp)
    pal = _palette(cfg)
    acfg = res.get("acfg") or {}
    amp_thr = float(acfg.get("amp_thr", 0.1))
    high_amp = float(acfg.get("high_amp", 0.9))

    n = len(d)
    ncols = min(3, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=theme.figsize(
        cfg, panel_h=2.2, total_w=7.0), squeeze=False)
    axes = axes.ravel()
    first_drawn = None
    for k, (_, row) in enumerate(d.iterrows()):
        ax = axes[k]
        drawn = False
        for xi, cid in ((-0.55, ca), (0.55, cb)):
            tab = (res.get("matrices") or {}).get(f"{cid}_vals_{unit}")
            if tab is None or row["unit"] not in tab.index:
                continue
            vals = tab.loc[row["unit"]].to_numpy(float)
            vals = vals[np.isfinite(vals)]
            if len(vals) == 0:
                continue
            rng = np.random.default_rng(k)
            x = xi + rng.normal(0, 0.14, len(vals))
            ax.scatter(x, vals, s=7, color=pal.get(cid, "#888888"),
                       alpha=0.65, edgecolor="none")
            ax.text(xi, ax.get_ylim()[0], theme.cohort_label(cfg, cid),
                    fontsize=5.5, ha="center", va="bottom",
                    color=pal.get(cid, "#888888"))
            drawn = True
        if not drawn:
            ax.set_visible(False)
            continue
        if first_drawn is None:
            first_drawn = ax
        ax.axhline(amp_thr, color=theme.OKABE["vermillion"], lw=0.6, ls="--")
        ax.axhline(high_amp, color=theme.OKABE["vermillion"], lw=0.6, ls=":")
        ax.axhline(-amp_thr, color=theme.OKABE["blue"], lw=0.6, ls="--")
        ax.axhline(-high_amp, color=theme.OKABE["blue"], lw=0.6, ls=":")
        ax.axhline(0, color="#999999", lw=0.5)
        ax.set_xlim(-1.1, 1.1)
        ax.set_xticks([])
        ax.set_title(str(row["unit"]), fontsize=6.5)
        ax.tick_params(labelsize=6)
        ax.grid(axis="y", alpha=0.2)
    for k in range(n, len(axes)):
        axes[k].set_visible(False)
    if first_drawn is not None:
        theme.annot_n(first_drawn, "\u2013 \u2013 altered (0.1)   "
                      "\u00b7 \u00b7 high-level (0.9)",
                      loc="upper right", fontsize=5.0)
    fig.suptitle(f"{comp}: per-sample copy value \u2014 top units by "
                 f"amplitude / frequency FDR ({unit})", fontsize=8)
    theme.finalize(fig, "tight")
    save(fig, cfg, f"{prefix}_{unit}_examples")


def figures_p7(res: dict, cfg) -> None:
    for unit in ("window", "band", "peak"):
        df = res.get(unit)
        if df is None or df.empty:
            continue
        if unit == "peak":
            # peak rows have no 'unit' column; label by peak_key
            df = df.rename(columns={"peak_key": "unit"})
            plot_p7_quadrant(df, cfg, unit, prefix="P7")
            plot_p7_volcano(df, cfg, unit, prefix="P7")
            continue
        plot_p7_quadrant(df, cfg, unit)
        plot_p7_volcano(df, cfg, unit)
        plot_p7_examples(res, cfg, unit)
