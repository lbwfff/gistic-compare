"""Figures for module P2 (call concordance).

(Split from the original single-module ``plotting``.)
"""

from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ._core import *  # noqa: F401,F403


def plot_p2_peak_agreement(peak_agreement: pd.DataFrame,
                           cfg: dict[str, Any]) -> None:
    if peak_agreement is None or peak_agreement.empty:
        return
    style(cfg)
    fig, axes = plt.subplots(1, 2, figsize=theme.figsize(
        cfg, panel_h=3.6), sharey=True)
    order = [p for p in ("Amp", "Del") if p in set(peak_agreement["ptype"])]
    for k, (ax, (metric, ylab)) in enumerate(zip(axes, (
            ("exact_agreement", "exact agreement"),
            ("weighted_kappa", "weighted κ")))):
        sns.violinplot(data=peak_agreement, x="ptype", y=metric, hue="ptype",
                       order=order, hue_order=order,
                       palette=theme.PTYPE_COLORS,
                       cut=0, inner="box", ax=ax, saturation=0.85,
                       legend=False)
        theme.strip_overlay(ax, peak_agreement, "ptype", metric, order=order,
                            hue="ptype", palette=theme.PTYPE_COLORS,
                            dodge=True, size=2.0)
        # metric as a panel title: with sharey the 2nd panel's ylabel is not
        # reliably rendered, which would leave the panel unlabelled
        ax.set_title(ylab, fontsize=8)
        ax.set_xlabel("")
        ax.set_ylabel("agreement" if k == 0 else "")
        ax.set_ylim(-0.05, 1.20)
        # annotate medians and n above each violin
        for i, pt in enumerate(order):
            g = peak_agreement.loc[peak_agreement["ptype"] == pt, metric]
            if g.empty:
                continue
            ax.text(i, 1.04, f"med {g.median():.2f}\n(n = {len(g)})",
                    ha="center", va="bottom", fontsize=5.5)
        comps_pa = list(dict.fromkeys(peak_agreement["comparison"]))
        if comps_pa:
            kn = _kind_note(cfg, comps_pa[0])
            if kn:
                theme.annot_n(ax, kn, loc="lower right", fontsize=5.5)
        theme.panel_label(ax, theme.panel_letter(k), cfg)
    theme.finalize(fig, "tight")
    save(fig, cfg, "P2_peak_agreement")


def plot_p2_sample_agreement(sample_agreement: pd.DataFrame,
                             cfg: dict[str, Any]) -> None:
    if sample_agreement is None or sample_agreement.empty:
        return
    style(cfg)
    comps = sorted(sample_agreement["comparison"].unique())
    nrows, ncols = theme.grid_shape(len(comps), max_cols=2)
    fig, axes = plt.subplots(nrows, ncols, squeeze=False, sharey=True,
                             figsize=(theme.width(cfg), 3.6 * nrows))
    flat = axes.ravel()
    for k, comp in enumerate(comps):
        ax = flat[k]
        g = sample_agreement[sample_agreement["comparison"] == comp]
        order = [p for p in ("Amp", "Del") if p in set(g["ptype"])]
        sns.violinplot(data=g, x="ptype", y="exact_agreement", hue="ptype",
                       order=order, hue_order=order,
                       palette=theme.PTYPE_COLORS,
                       cut=0, inner="box", ax=ax, saturation=0.85,
                       legend=False)
        theme.strip_overlay(ax, g, "ptype", "exact_agreement", order=order,
                            hue="ptype", palette=theme.PTYPE_COLORS,
                            dodge=True, size=2.0)
        ax.set_title(comp, fontsize=8)
        ax.set_xlabel("")
        ax.set_ylabel("exact agreement")
        ax.set_ylim(-0.05, 1.20)
        theme.annot_n(ax, theme.group_sizes(g, "ptype", order) +
                      (f"\n{_n_pair(cfg, comp)}" if _n_pair(cfg, comp) else ""),
                      loc="upper left")
        theme.panel_label(ax, theme.panel_letter(k), cfg)
    for ax in flat[len(comps):]:
        ax.set_visible(False)
    theme.finalize(fig, "constrained")
    save(fig, cfg, "P2_sample_agreement")


def plot_p2_mcnemar_volcano(mcnemar: pd.DataFrame,
                            peak_agreement: pd.DataFrame,
                            cfg: dict[str, Any]) -> None:
    if mcnemar is None or mcnemar.empty:
        return
    style(cfg)
    m = mcnemar.copy()
    if peak_agreement is not None and not peak_agreement.empty:
        m = m.merge(
            peak_agreement[["comparison", "ptype", "a_peak_id",
                            "a_freq", "b_freq"]],
            on=["comparison", "ptype", "a_peak_id"], how="left")
        m["delta_freq"] = m["a_freq"] - m["b_freq"]
    else:
        m["delta_freq"] = np.nan
    # single scatter → single-column width (avoid flat full-width panel)
    fig, ax = plt.subplots(figsize=theme.figsize(
        cfg, panel_h=4.0, total_w=5.2))
    xs_all, ys_all, labs_all = [], [], []
    for ptype, col in (("Amp", _pc("Amp")), ("Del", _pc("Del"))):
        g = m[m["ptype"] == ptype]
        if g.empty or g["delta_freq"].isna().all():
            continue
        y = -np.log10(np.clip(g["mcnemar_fdr"].to_numpy(float), 1e-300, None))
        ax.scatter(g["delta_freq"], y, s=18, alpha=0.75, color=col,
                   label=ptype, edgecolors="none")
        # peak label: a_cytoband (append b if boundaries differ); fall back to
        # the peak id so a missing cytoband never drops the label entirely
        for (_, r), yi in zip(g.iterrows(), y):
            la = _clean_label(r.get("a_cytoband"))
            lb = _clean_label(r.get("b_cytoband"))
            if la and lb and la != lb:
                lab = f"{la}/{lb}"
            else:
                lab = la or lb or _clean_label(r.get("a_peak_id"))
            xs_all.append(float(r["delta_freq"]))
            ys_all.append(float(yi))
            labs_all.append(lab)
    _annotate_labeled_points(ax, xs_all, ys_all, labs_all,
                             fontsize=5.0, min_sep_pt=17)
    if not m["delta_freq"].isna().all():
        ax.axvline(0, color="#888888", ls="--", lw=0.7)
    ax.axhline(-np.log10(0.05), color="#888888", ls=":", lw=0.7)
    if not m.empty:
        comp0 = str(m["comparison"].iloc[0])
        la, lb, _, _ = _sides(cfg, comp0)
        npair = _n_pair(cfg, comp0)
        xl = f"Δ frequency ({la} − {lb})"
        if npair:
            xl += f"   (n = {len(m)} peaks, {npair})"
        ax.set_xlabel(xl)
    from matplotlib.lines import Line2D
    h, _lab = ax.get_legend_handles_labels()
    h += [Line2D([], [], ls=":", color="#888888", lw=0.7,
                 label="FDR = 0.05")]
    ax.legend(handles=h, frameon=False, fontsize=6.5, loc="upper left")
    theme.panel_label(ax, theme.panel_letter(0), cfg)
    theme.finalize(fig, "tight")
    save(fig, cfg, "P2_mcnemar_volcano")


def plot_p2_peak_kappa_sorted(peak_agreement: pd.DataFrame,
                              cfg: dict[str, Any]) -> None:
    """Per-peak weighted κ sorted ascending (find worst-agreement peaks)."""
    if peak_agreement is None or peak_agreement.empty:
        return
    style(cfg)
    comps = sorted(peak_agreement["comparison"].unique())
    nrows, ncols = theme.grid_shape(len(comps), max_cols=2)
    fig, axes = plt.subplots(nrows, ncols, squeeze=False, sharey=True,
                             figsize=(theme.width(cfg), 4.4 * nrows))
    flat = axes.ravel()
    max_len = 0
    for k, comp in enumerate(comps):
        ax = flat[k]
        g = peak_agreement[peak_agreement["comparison"] == comp].copy()
        if g.empty:
            ax.set_visible(False)
            continue
        g["label"] = _label_series(g, ["a_cytoband", "b_cytoband"],
                                   "a_peak_id")
        g = g.sort_values("weighted_kappa").reset_index(drop=True)
        max_len = max(max_len, len(g))
        y = np.arange(len(g))
        for ptype in ("Amp", "Del"):
            mask = (g["ptype"] == ptype).to_numpy()
            if not mask.any():
                continue
            ax.scatter(g.loc[mask, "weighted_kappa"], y[mask], s=16,
                       color=_pc(ptype), label=ptype,
                       edgecolors="none", zorder=3)
        # label 3 worst overall
        for yi in range(min(3, len(g))):
            ax.text(g["weighted_kappa"].iloc[yi], yi,
                    f" {g['label'].iloc[yi]}",
                    fontsize=5, va="center", color="#333333")
        ax.axvline(0.8, color="#AAAAAA", ls=":", lw=0.7)
        ax.set_title(comp, fontsize=8)
        ax.set_xlabel("weighted κ")
        ax.set_xlim(0, 1.02)
        ax.set_ylim(-0.8, len(g) - 0.2)
        ax.tick_params(labelsize=6)
        n_note = f"n = {len(g)} peaks"
        n_pairs = _n_pair(cfg, comp)
        theme.annot_n(ax, n_note + ("\n" + n_pairs if n_pairs else ""),
                      loc="upper left")
        theme.panel_label(ax, theme.panel_letter(k), cfg)
    for ax in flat[len(comps):]:
        ax.set_visible(False)
    for ax in flat[:len(comps)]:
        ax.set_ylim(-0.8, max(max_len - 0.2, 1.0))
    flat[0].set_ylabel("peaks (ranked by κ)")
    if flat[0].get_legend_handles_labels()[0]:
        flat[0].legend(frameon=False, fontsize=6, loc="lower right")
    theme.finalize(fig, "constrained")
    save(fig, cfg, "P2_peak_kappa_sorted")


def plot_p2_discordant_stacked(mcnemar: pd.DataFrame,
                               cfg: dict[str, Any]) -> None:
    """Top discordant peaks: samples altered only in A vs only in B.

    Bars are coloured by cohort (A / B) so the legend matches what is drawn;
    the peak type stays readable from the y tick label.
    """
    if mcnemar is None or mcnemar.empty:
        return
    style(cfg)
    m = mcnemar.copy()
    m["total_disc"] = m["only_a_altered"] + m["only_b_altered"]
    comps = sorted(m["comparison"].unique())
    nrows, ncols = theme.grid_shape(len(comps), max_cols=2)
    fig, axes = plt.subplots(nrows, ncols, squeeze=False, sharex=True,
                             figsize=(theme.width(cfg), 4.4 * nrows))
    flat = axes.ravel()
    from matplotlib.patches import Patch
    for k, comp in enumerate(comps):
        ax = flat[k]
        la, lb, cola, colb = _sides(cfg, comp)
        g = m[m["comparison"] == comp]
        if g.empty:
            ax.set_visible(False)
            continue
        rows = []
        for ptype in ("Amp", "Del"):
            s = g[g["ptype"] == ptype].nlargest(8, "total_disc")
            for _, r in s.iterrows():
                rows.append((ptype, r))
        rows = rows[::-1]  # highest on top after barh
        ypos = np.arange(len(rows))
        a_vals = np.array([float(r["only_a_altered"]) for _, r in rows])
        b_vals = np.array([float(r["only_b_altered"]) for _, r in rows])
        ax.barh(ypos, a_vals, color=cola, edgecolor="white", linewidth=0.4,
                height=0.7)
        ax.barh(ypos, b_vals, left=a_vals, color=colb, edgecolor="white",
                linewidth=0.4, height=0.7, hatch="//")
        labels = [f"{_clean_label(r.get('a_cytoband')) or _clean_label(r.get('a_peak_id'))}"
                  f" ({pt})" for pt, r in rows]
        for yi, (pt, r) in enumerate(rows):
            fdr = r.get("mcnemar_fdr", np.nan)
            star = "*" if np.isfinite(fdr) and fdr < 0.05 else ""
            if star:
                ax.text(a_vals[yi] + b_vals[yi], yi, star, va="center",
                        ha="left", fontsize=7, color="#333333")
        ax.set_yticks(ypos)
        ax.set_yticklabels(labels, fontsize=5.5)
        ax.set_title(comp, fontsize=8)
        n_pair = _n_pair(cfg, comp)
        ax.set_xlabel("discordant samples"
                      + (f"   ({n_pair} paired)" if n_pair else ""))
        ax.tick_params(labelsize=6)
        ax.legend(handles=[
            Patch(facecolor=cola, edgecolor="white", label=f"only {la}"),
            Patch(facecolor=colb, edgecolor="white", hatch="//",
                  label=f"only {lb}"),
        ], frameon=False, fontsize=6, loc="lower right")
        theme.panel_label(ax, theme.panel_letter(k), cfg)
    for ax in flat[len(comps):]:
        ax.set_visible(False)
    theme.finalize(fig, "constrained", wspace=0.30)
    save(fig, cfg, "P2_discordant_stacked")


def plot_p2_sample_ranked(sample_agreement: pd.DataFrame,
                          cfg: dict[str, Any]) -> None:
    """Samples ranked by agreement (identify outlier samples)."""
    if sample_agreement is None or sample_agreement.empty:
        return
    style(cfg)
    comps = sorted(sample_agreement["comparison"].unique())
    nrows, ncols = theme.grid_shape(len(comps), max_cols=2)
    fig, axes = plt.subplots(nrows, ncols, squeeze=False, sharey=True,
                             figsize=(theme.width(cfg), 4.0 * nrows))
    flat = axes.ravel()
    max_len = 0
    for k, comp in enumerate(comps):
        ax = flat[k]
        g = sample_agreement[sample_agreement["comparison"] == comp].copy()
        if g.empty:
            ax.set_visible(False)
            continue
        # mean across ptypes per sample
        agg = (g.groupby("sample", as_index=False)
               .agg(acc=("exact_agreement", "mean"),
                    kappa=("weighted_kappa", "mean")))
        agg = agg.sort_values("acc").reset_index(drop=True)
        max_len = max(max_len, len(agg))
        y = np.arange(len(agg))
        ax.scatter(agg["acc"], y, s=8, alpha=0.7, color=theme.OKABE["blue"],
                   edgecolors="none")
        # mark bottom 5
        for i in range(min(5, len(agg))):
            ax.text(agg["acc"].iloc[i], i, f" {agg['sample'].iloc[i]}",
                    fontsize=4.5, va="center", color="#333333")
        med = float(agg["acc"].median())
        ax.axvline(med, color="#888888", ls="--", lw=0.7)
        ax.text(med, len(agg) * 0.02, f" med {med:.2f}",
                fontsize=5.5, color="#555555")
        ax.set_xlim(-0.02, 1.02)
        ax.set_title(comp, fontsize=8)
        ax.set_xlabel("exact agreement")
        ax.set_ylim(-1, len(agg))
        ax.tick_params(labelsize=6)
        theme.annot_n(ax, f"n = {len(agg)} samples", loc="upper left")
        theme.panel_label(ax, theme.panel_letter(k), cfg)
    for ax in flat[len(comps):]:
        ax.set_visible(False)
    for ax in flat[:len(comps)]:
        ax.set_ylim(-1, max(max_len, 1))
    flat[0].set_ylabel("samples (ranked)")
    theme.finalize(fig, "constrained")
    save(fig, cfg, "P2_sample_ranked")


def figures_p2(res: dict, cfg) -> None:
    plot_p2_peak_agreement(res.get("peak_agreement"), cfg)
    plot_p2_sample_agreement(res.get("sample_agreement"), cfg)
    plot_p2_mcnemar_volcano(res.get("mcnemar"), res.get("peak_agreement"),
                            cfg)
    plot_p2_peak_kappa_sorted(res.get("peak_agreement"), cfg)
    plot_p2_discordant_stacked(res.get("mcnemar"), cfg)
    plot_p2_sample_ranked(res.get("sample_agreement"), cfg)


