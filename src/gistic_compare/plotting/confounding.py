"""Figures for module CF (confounder sensitivity).

(Split from the original single-module ``plotting``.)
"""

from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ._core import *  # noqa: F401,F403

# ---------------------------------------------------------------------------
# CF: confounder sensitivity (independent of P4/P5)
# ---------------------------------------------------------------------------

def _preset_bands(cfg) -> set:
    from ..presets import band_names
    return band_names(cfg)


def _cf_max_binary(cfg) -> float:
    return float((cfg.get("confounding") or {}).get(
        "plausible_max_binary", 0.20))


def plot_cf_required_effect(sens: pd.DataFrame, cfg: dict[str, Any]) -> None:
    """|observed Δfreq| vs effect a confounder would need to explain it."""
    if sens is None or sens.empty:
        return
    d = sens[sens["required_unit"] == "risk_diff"].copy()
    if d.empty:
        return
    style(cfg)
    maxb = _cf_max_binary(cfg)
    box = ["sex", "smoking", "mut_EGFR", "mut_KRAS", "mut_TP53"]
    order = [c for c in box if c in set(d["confounder"])]
    order += sorted(set(d["confounder"]) - set(order))
    cols = theme.categorical(len(order))
    cmap = {c: cols[i] for i, c in enumerate(order)}
    fig, ax = plt.subplots(figsize=theme.figsize(cfg, panel_h=3.6, total_w=5.2))
    ax.axhline(maxb, color=theme.OKABE["vermillion"], ls="--", lw=1.0,
               zorder=1)
    for c in order:
        s = d[d["confounder"] == c]
        ax.scatter(np.abs(s["diff"]), np.abs(s["required_effect"]),
                   s=14, alpha=0.7, color=cmap[c], label=c,
                   edgecolor="none", zorder=2)
    # one label per paper locus, at its largest |Dfreq|: iterating every row
    # would re-plot the same 4 names thousands of times on top of each other
    pap = d[d["band"].isin(_preset_bands(cfg))].copy()
    if not pap.empty:
        pap["_ad"] = pap["diff"].abs()
        pap = (pap.sort_values("_ad", ascending=False)
                  .drop_duplicates("band"))
        for k, (_, r) in enumerate(pap.iterrows()):
            dx, dy = (5, 4) if k % 2 == 0 else (5, -9)
            ax.annotate(str(r["band"]),
                        (abs(float(r["diff"])), abs(float(r["required_effect"]))),
                        fontsize=5.5, color="#333333",
                        xytext=(dx, dy), textcoords="offset points", zorder=3)
    comp0 = str(d["comparison"].iloc[0]) if "comparison" in d.columns else ""
    la, lb, _, _ = _sides(cfg, comp0)
    ax.set_yscale("log")
    ax.set_xlabel(f"observed |Δfreq|  ({la} − {lb})" if comp0
                  else "observed |Δfreq|")
    ax.set_ylabel("required confounder effect  δ* = |Δfreq| / |Δp|")
    ax.set_title("required confounder effect", fontsize=8)
    ax.legend(frameon=False, fontsize=6, loc="upper left")
    fig.text(0.5, 0.005,
             f"binary confounders; dashed line = plausible max δ* = {maxb:g}; "
             "above it the finding is robust",
             fontsize=5.5, color="#555555", ha="center", va="bottom")
    theme.annot_n(ax, f"n = {len(d)} findings", loc="lower right")
    fig.tight_layout(rect=(0.0, 0.035, 1.0, 1.0))
    save(fig, cfg, "CF_required_effect")


def plot_cf_evalue_volcano(ev_rr: pd.DataFrame, cfg: dict[str, Any]) -> None:
    """Tier-1 generic E-value vs effect size for every significant finding."""
    if ev_rr is None or ev_rr.empty:
        return
    style(cfg)
    d = ev_rr.copy()
    d["y"] = -np.log10(np.clip(d["mw_fdr"].astype(float), 1e-300, 1.0))
    d = d[np.isfinite(d["evalue_rr"]) & np.isfinite(d["y"])]
    if d.empty:
        return
    fig, ax = plt.subplots(figsize=theme.figsize(cfg, panel_h=3.6, total_w=5.4))
    sc = ax.scatter(d["rr"], d["y"], c=d["evalue_rr"].astype(float),
                    cmap="viridis", s=16, alpha=0.8, edgecolor="none")
    ax.axvline(1.0, color="#888888", lw=0.8, ls=":")
    for _, r in d[d["band"].isin(_preset_bands(cfg))].iterrows():
        ax.annotate(str(r["band"]), (float(r["rr"]), float(r["y"])),
                    fontsize=5.5, xytext=(3, 3), textcoords="offset points")
    cbar = fig.colorbar(sc, ax=ax, shrink=0.6)
    cbar.set_label("E-value (RR scale)", fontsize=7)
    ax.set_xscale("log")
    comp0 = str(d["comparison"].iloc[0]) if "comparison" in d.columns else ""
    la, lb, _, _ = _sides(cfg, comp0)
    ax.set_xlabel(f"risk ratio  ({la} vs {lb})" if comp0
                  else "risk ratio")
    ax.set_ylabel("−log10  MW FDR")
    ax.set_title("generic E-value per finding (Tier 1)", fontsize=8)
    theme.annot_n(ax, f"n = {len(d)} findings", loc="lower right")
    theme.finalize(fig, "tight")
    save(fig, cfg, "CF_evalue_volcano")


def plot_cf_bias_tornado(sens: pd.DataFrame, cfg: dict[str, Any]) -> None:
    """Required δ* per confounder for the four paper loci (log x-axis).

    A lollipop (not bars) so the log scale is honest: each point is a required
    effect δ*; confounders are ordered by cohort exposure gap |Δp| (largest =
    most leverage first). Points right of the dashed line are robust.
    """
    if sens is None or sens.empty:
        return
    d = sens[(sens["required_unit"] == "risk_diff")
             & (sens["band"].isin(_preset_bands(cfg)))
             & (sens["level"] == "band")].copy()
    if d.empty:
        return
    style(cfg)
    gap = {c: abs(float(g["delta_p"].iloc[0]))
           for c, g in d.groupby("confounder")}
    confs = sorted(gap, key=lambda c: -gap[c])          # most leverage first
    d["label"] = d["band"].astype(str) + " (" + d["direction"].str[0] + ")"
    labels = list(dict.fromkeys(d["label"]))
    maxb = _cf_max_binary(cfg)
    ypos = np.arange(len(labels))
    h = 0.8 / max(1, len(confs))
    cols = theme.categorical(len(confs))
    fig, ax = plt.subplots(figsize=theme.figsize(
        cfg, panel_h=0.5 * len(labels) + 1.8, total_w=6.2))
    allx = []
    for i, c in enumerate(confs):
        xs, ys = [], []
        for lab in labels:
            s = d[(d["label"] == lab) & (d["confounder"] == c)]
            if not len(s):
                continue
            x = abs(float(s["required_effect"].iloc[0]))
            if not np.isfinite(x) or x <= 0:
                continue
            xs.append(x)
            ys.append(ypos[labels.index(lab)]
                      + (i - (len(confs) - 1) / 2) * h)
        allx += xs
        ax.scatter(xs, ys, s=24, color=cols[i], label=f"{c} (Δp={gap[c]:.2f})",
                   zorder=3, edgecolor="white", linewidth=0.4)
    ax.axvline(maxb, color=theme.OKABE["vermillion"], ls="--", lw=1.0,
               label=f"plausible max ({maxb:g})", zorder=2)
    ax.set_xscale("log")
    lo = min([v for v in allx if v > 0] + [maxb]) * 0.6
    ax.set_xlim(lo, max(allx) * 1.8)
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=6.5)
    ax.invert_yaxis()
    ax.set_xlabel("required confounder effect  \u03b4*   (log scale)")
    ax.set_title("paper loci \u2014 effect needed per confounder\n"
                 "(right of dashed line => finding is robust)", fontsize=8)
    ax.legend(frameon=False, fontsize=5.5, ncol=2, loc="lower left")
    ax.grid(axis="x", alpha=0.25, which="both")
    theme.finalize(fig, "tight")
    save(fig, cfg, "CF_bias_tornado")


def figures_confounding(res: dict, cfg) -> None:
    plot_cf_required_effect(res.get("sensitivity"), cfg)
    plot_cf_evalue_volcano(res.get("evalue_rr"), cfg)
    plot_cf_bias_tornado(res.get("sensitivity"), cfg)


