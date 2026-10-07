"""Figures for module A (peak overlap).

(Split from the original single-module ``plotting``.)
"""

from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ._core import *  # noqa: F401,F403

# ---------------------------------------------------------------------------
# A: peak overlap
# ---------------------------------------------------------------------------

def plot_jaccard(jaccard_df: pd.DataFrame, cfg: dict[str, Any]) -> None:
    if jaccard_df.empty:
        return
    style(cfg)
    q_main = cfg["analysis"]["q_thresholds"][0]
    show = jaccard_df[
        jaccard_df["method"].str.startswith("interval")
        & (jaccard_df["q_thresh"] == q_main)
    ].copy()
    if show.empty:
        return
    show["tier"] = show["method"].str.replace("interval_", "", regex=False)
    tiers = [t for t in ("moderate", "strict", "loose")
             if t in set(show["tier"])]
    comps = list(dict.fromkeys(show["comparison"]))
    ptypes = ["Amp", "Del"]
    n = len(tiers)
    comp_labels = []
    for comp in comps:
        s = show[(show["comparison"] == comp) & (show["tier"] == tiers[0])]
        if len(s):
            r = s.iloc[0]
            comp_labels.append(
                f"{comp}\n({int(r['n_a'])}/{int(r['n_b'])}, "
                f"{int(r['n_pairs'])})")
        else:
            comp_labels.append(comp)
    fig, axes = plt.subplots(
        n, 1, sharex=True, squeeze=False,
        figsize=theme.figsize(cfg, nrows=n, panel_h=1.9))
    x = np.arange(len(comps))
    w = 0.36
    ymax = float(show["jaccard"].max()) if not show.empty else 0.5
    for ti, tier in enumerate(tiers):
        ax = axes[ti][0]
        for pi, pt in enumerate(ptypes):
            vals = []
            for comp in comps:
                s = show[(show["tier"] == tier)
                         & (show["comparison"] == comp)
                         & (show["ptype"] == pt)]
                vals.append(float(s["jaccard"].iloc[0]) if len(s) else np.nan)
            bars = ax.bar(x + (pi - 0.5) * w, vals, w, color=_pc(pt),
                          label=pt, edgecolor="white", linewidth=0.4)
            ax.bar_label(bars, fmt="%.2f", fontsize=5.5, padding=1,
                         color="#444444")
        ax.set_ylim(0, max(0.5, ymax * 1.25))
        ax.set_ylabel("Jaccard")
        ax.text(0.01, 0.96, tier, transform=ax.transAxes, fontsize=7,
                fontstyle="italic", ha="left", va="top")
    axes[-1][0].set_xticks(x)
    axes[-1][0].set_xticklabels(comp_labels, rotation=15, ha="right")
    axes[-1][0].set_xlabel("comparison  (n = peaks A / peaks B, matched)")
    axes[0][0].legend(loc="upper right", ncol=2)
    for ti in range(n):
        theme.panel_label(axes[ti][0], theme.panel_letter(ti), cfg)
    theme.finalize(fig, "tight")
    save(fig, cfg, "A_jaccard_interval")


def plot_q_scatter(matched: pd.DataFrame, cfg: dict[str, Any],
                   corr: pd.DataFrame | None = None) -> None:
    """Scatter of rank-normalized -log10(q) (z) with Pearson bootstrap CI."""
    if matched.empty:
        return
    style(cfg)
    q_main = cfg["analysis"]["q_thresholds"][0]
    m = matched[matched["q_thresh"] == q_main].copy()
    if m.empty or "a_z" not in m.columns:
        return
    z_ci: dict[tuple[str, str], tuple[float, float, float, int]] = {}
    if corr is not None and not corr.empty:
        zsub = corr[(corr["q_thresh"] == q_main)
                    & (corr["metric"] == "ranknorm_neg_log10_q")]
        for _, row in zsub.iterrows():
            z_ci[(row["comparison"], row["ptype"])] = (
                row["corr"], row["ci_lo"], row["ci_hi"], int(row["n_pairs"]))

    comps = list(dict.fromkeys(m["comparison"]))
    ptypes = ["Amp", "Del"]
    allz = np.concatenate([m["a_z"].to_numpy(float), m["b_z"].to_numpy(float)])
    glo, ghi = float(np.nanmin(allz)), float(np.nanmax(allz))
    gpad = (ghi - glo) * 0.06 + 0.15
    glo, ghi = glo - gpad, ghi + gpad
    mid = (glo + ghi) / 2
    fig, axes = plt.subplots(
        len(comps), len(ptypes), squeeze=False, sharex=True, sharey=True,
        figsize=theme.figsize(cfg, nrows=len(comps), panel_h=2.6))
    for i, comp in enumerate(comps):
        la, lb, _, _ = _sides(cfg, comp)
        for j, ptype in enumerate(ptypes):
            ax = axes[i][j]
            sub = m[(m["comparison"] == comp) & (m["a_type"] == ptype)]
            if sub.empty:
                ax.set_visible(False)
                continue
            az = sub["a_z"].to_numpy(float)
            bz = sub["b_z"].to_numpy(float)
            ax.scatter(az, bz, s=16, alpha=0.8, c=_pc(ptype),
                       edgecolors="none", zorder=3)
            ax.plot([glo, ghi], [glo, ghi], "k--", lw=0.7, alpha=0.5,
                    zorder=1)
            if np.isfinite(az).sum() >= 3 and np.std(az) > 0:
                k, b = np.polyfit(az, bz, 1)
                xs = np.linspace(np.nanmin(az), np.nanmax(az), 50)
                ax.plot(xs, k * xs + b, color="#888888", lw=0.9, zorder=2)
            if (comp, ptype) in z_ci:
                r, ci_lo, ci_hi, n = z_ci[(comp, ptype)]
            else:
                from ..peak_overlap import bootstrap_pearson
                boot = cfg["analysis"].get("bootstrap", {})
                bp = bootstrap_pearson(
                    az, bz, n_boot=int(boot.get("n", 2000)),
                    seed=int(boot.get("seed", 42)),
                    ci=float(boot.get("ci", 0.95)))
                r, ci_lo, ci_hi, n = (bp["r"], bp["ci_lo"], bp["ci_hi"], bp["n"])
            sub2 = sub.copy()
            sub2["d"] = (sub2["a_z"] - sub2["b_z"]).abs()
            top = sub2.nlargest(1, "d")
            stat_loc = ("lower right"
                        if (len(top) and top.iloc[0]["b_z"] >= mid)
                        else "upper left")
            theme.annot_stat(ax, theme.fmt_ci("r", r, ci_lo, ci_hi, n),
                             loc=stat_loc)
            for _, row in top.iterrows():
                if np.isfinite(row["a_z"]) and np.isfinite(row["b_z"]):
                    dx = 5 if row["a_z"] <= mid else -6
                    dy = -8 if row["b_z"] >= mid else 5
                    ax.annotate(str(row.get("a_cytoband", "")),
                                (row["a_z"], row["b_z"]), fontsize=5.5,
                                alpha=0.85, textcoords="offset points",
                                xytext=(dx, dy))
            ax.set_xlim(glo, ghi)
            ax.set_ylim(glo, ghi)
            ax.tick_params(axis="x", pad=5)
            ax.tick_params(axis="y", pad=2)
            ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
            ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
            if i == 0:
                ax.set_title(ptype)
                theme.panel_label(ax, theme.panel_letter(j), cfg)
            if j == 0:
                ax.set_ylabel(f"{comp}\n{lb}: z(-log10 q)")
            if i == len(comps) - 1:
                ax.set_xlabel(f"{la}: z(-log10 q)")
    theme.finalize(fig, "constrained", wspace=0.35, hspace=0.30)
    save(fig, cfg, "A_matched_q_scatter")


def plot_matched_frequency(matched: pd.DataFrame, cfg: dict[str, Any],
                           corr: pd.DataFrame | None = None) -> None:
    """Frequency-frequency scatter with bootstrap Spearman CI."""
    if matched.empty or "a_freq_altered" not in matched.columns:
        return
    style(cfg)
    q_main = cfg["analysis"]["q_thresholds"][0]
    m = matched[matched["q_thresh"] == q_main].copy()
    if m.empty:
        return
    f_ci: dict[tuple[str, str], tuple[float, float, float, float, int]] = {}
    if corr is not None and not corr.empty:
        fsub = corr[(corr["q_thresh"] == q_main)
                    & (corr["metric"] == "freq_altered")
                    & (corr["corr_type"] == "spearman_bootstrap")]
        for _, row in fsub.iterrows():
            f_ci[(row["comparison"], row["ptype"])] = (
                row["corr"], row["ci_lo"], row["ci_hi"], row["pvalue"],
                int(row["n_pairs"]))

    comps = list(dict.fromkeys(m["comparison"]))
    ptypes = ["Amp", "Del"]
    fig, axes = plt.subplots(
        len(comps), len(ptypes), squeeze=False, sharex=True, sharey=True,
        figsize=theme.figsize(cfg, nrows=len(comps), panel_h=2.6))
    for i, comp in enumerate(comps):
        la, lb, _, _ = _sides(cfg, comp)
        for j, ptype in enumerate(ptypes):
            ax = axes[i][j]
            sub = m[(m["comparison"] == comp) & (m["a_type"] == ptype)]
            if sub.empty:
                ax.set_visible(False)
                continue
            af = sub["a_freq_altered"].to_numpy(float)
            bf = sub["b_freq_altered"].to_numpy(float)
            ax.scatter(af, bf, s=16, alpha=0.8, c=_pc(ptype),
                       edgecolors="none", zorder=3)
            ax.plot([0, 1], [0, 1], "k--", lw=0.7, alpha=0.5, zorder=1)
            if np.isfinite(af).sum() >= 3 and np.std(af) > 0:
                k, b = np.polyfit(af, bf, 1)
                xs = np.linspace(0, 1, 50)
                ax.plot(xs, k * xs + b, color="#888888", lw=0.9, zorder=2)
            if (comp, ptype) in f_ci:
                r, ci_lo, ci_hi, p, n = f_ci[(comp, ptype)]
            else:
                from ..peak_overlap import bootstrap_spearman
                boot = cfg["analysis"].get("bootstrap", {})
                fs = bootstrap_spearman(
                    af, bf, n_boot=int(boot.get("n", 2000)),
                    seed=int(boot.get("seed", 42)),
                    ci=float(boot.get("ci", 0.95)))
                r, ci_lo, ci_hi, p, n = (fs["r"], fs["ci_lo"], fs["ci_hi"],
                                         fs["p"], fs["n"])
            theme.annot_stat(ax, theme.fmt_ci("\u03c1", r, ci_lo, ci_hi, n, p),
                             loc="lower right")
            ax.set_xlim(-0.02, 1.02)
            ax.set_ylim(-0.02, 1.02)
            ticks = [0.0, 0.25, 0.5, 0.75, 1.0]
            ax.set_xticks(ticks)
            ax.set_yticks(ticks)
            ax.tick_params(pad=2)
            if i == 0:
                ax.set_title(ptype)
                theme.panel_label(ax, theme.panel_letter(j), cfg)
            if j == 0:
                ax.set_ylabel(f"{comp}\n{lb}: frequency altered")
            if i == len(comps) - 1:
                ax.set_xlabel(f"{la}: frequency altered")
    theme.finalize(fig, "constrained", wspace=0.30, hspace=0.26)
    save(fig, cfg, "A_matched_frequency_scatter")


def plot_permutation(perm: pd.DataFrame, cfg: dict[str, Any],
                     null_samples: list | None = None) -> None:
    """Histogram of permutation null Jaccard with observed line."""
    if perm.empty:
        return
    style(cfg)
    q_main = cfg["analysis"]["q_thresholds"][0]
    sub = perm[perm["q_thresh"] == q_main].copy()
    if sub.empty:
        return
    null_map: dict[tuple[str, str], dict] = {}
    for item in null_samples or []:
        if item.get("q_thresh") == q_main:
            null_map[(item["comparison"], item["ptype"])] = item

    comps = list(dict.fromkeys(sub["comparison"]))
    ptypes = ["Amp", "Del"]
    fig, axes = plt.subplots(
        len(comps), len(ptypes), squeeze=False, sharex=False,
        figsize=theme.figsize(cfg, nrows=len(comps), panel_h=2.1))
    for i, comp in enumerate(comps):
        for j, ptype in enumerate(ptypes):
            ax = axes[i][j]
            row_df = sub[(sub["comparison"] == comp) & (sub["ptype"] == ptype)]
            if row_df.empty:
                ax.set_visible(False)
                continue
            row = row_df.iloc[0]
            key = (comp, ptype)
            if key in null_map:
                null_j = np.asarray(null_map[key]["null_jaccards"], dtype=float)
                obs = float(null_map[key]["observed_jaccard"])
            else:
                null_j = np.array([])
                obs = float(row["observed_jaccard"])
            if len(null_j):
                ax.hist(null_j, bins=30, color="#B8B8B8", edgecolor="#8A8A8A",
                        density=True, alpha=0.9, zorder=1)
                ax.axvline(np.mean(null_j), color="#333333", lw=0.9, ls=":",
                           label="null mean", zorder=2)
            ax.axvline(obs, color=_pc(ptype), lw=1.6, label="observed",
                       zorder=3)
            # x-range covers the null bulk plus observed, with a little margin
            if len(null_j):
                nmin, nmax = np.quantile(null_j, [0.01, 0.99])
                nmin, nmax = float(nmin), float(nmax)
            else:
                nmin = nmax = obs
            lo = min(nmin, obs)
            hi = max(nmax, obs)
            pad = (hi - lo) * 0.06 + 0.002
            ax.set_xlim(lo - pad, hi + pad)
            theme.annot_stat(
                ax, f"obs = {obs:.3f}\n{theme.fmt_p(row.pvalue)}\n"
                    f"n_perm = {int(row.n_perm)}", loc="upper left")
            ax.tick_params(pad=2)
            if i == 0:
                ax.set_title(ptype)
                theme.panel_label(ax, theme.panel_letter(j), cfg)
            if i == len(comps) - 1:
                ax.set_xlabel("Jaccard")
            if j == 0:
                ax.set_ylabel(f"{comp}\ndensity")
    axes[0][0].legend(loc="upper right", fontsize=6)
    theme.finalize(fig, "constrained", wspace=0.30, hspace=0.35)
    save(fig, cfg, "A_permutation_null")


def plot_three_way(upset_input: pd.DataFrame, cfg: dict[str, Any],
                   *, prefix: str = "A",
                   name_pattern: str | None = None) -> None:
    """UpSet plot (bar + matrix + set sizes) of peak groups by cohort membership.

    ``name_pattern`` format fields: {prefix}, {ptype}, {ptype_l}.
    Default: ``{prefix}_three_way_membership_{ptype}``.
    """
    if upset_input.empty:
        return
    style(cfg)
    pat = name_pattern or "{prefix}_three_way_membership_{ptype}"
    cohorts = list(cfg["cohorts"].keys())
    label_map = {c: cfg["cohorts"][c].get("label", c) for c in cohorts}
    for ptype in ("Amp", "Del"):
        sub = upset_input[upset_input["type"] == ptype].copy()
        if sub.empty:
            continue
        membership = sub["cohorts"].str.split(",")
        combos = membership.apply(lambda s: tuple(sorted(s)))
        vc = combos.value_counts()
        order = sorted(vc.index, key=lambda t: (-len(t), -vc[t]))
        sizes = [vc[c] for c in order]
        set_totals = {c: int(membership.apply(lambda s, c=c: c in s).sum())
                      for c in cohorts}

        n_sub = len(order)
        n_set = len(cohorts)
        fig = plt.figure(figsize=(theme.width(cfg), 4.2))
        gs = fig.add_gridspec(2, 2, width_ratios=[3, 1], height_ratios=[1, 2],
                              hspace=0.08, wspace=0.08)
        ax_bar = fig.add_subplot(gs[0, 0])
        ax_mat = fig.add_subplot(gs[1, 0], sharex=ax_bar)
        ax_set = fig.add_subplot(gs[1, 1], sharey=ax_mat)

        xs = np.arange(n_sub)
        bar_c = theme.OKABE["blue"]
        set_c = theme.OKABE["orange"]
        ax_bar.bar(xs, sizes, color=bar_c, width=0.7)
        ax_bar.bar_label(ax_bar.containers[0], fmt="%d", fontsize=5.5,
                         padding=1, color="#444444")
        ax_bar.set_ylim(0, max(sizes) * 1.20)
        ax_bar.set_ylabel("n groups")
        ax_bar.tick_params(axis="x", labelbottom=False)
        ax_bar.text(0.01, 1.02, f"{ptype}  (n = {len(sub)} groups)",
                    transform=ax_bar.transAxes,
                    fontsize=8, ha="left", va="bottom")

        for i, combo in enumerate(order):
            for j, c in enumerate(cohorts):
                on = c in combo
                ax_mat.scatter(i, j, s=80 if on else 34,
                               c=bar_c if on else "#DDDDDD",
                               edgecolors="#333333" if on else "#BBBBBB",
                               linewidths=0.6, zorder=3)
            if len(order) > 1:
                on_rows = [j for j, c in enumerate(cohorts) if c in combo]
                if len(on_rows) > 1:
                    ax_mat.plot([i, i], [min(on_rows), max(on_rows)],
                                color="#888888", lw=1.0, zorder=2)
        ax_mat.set_yticks(range(n_set))
        ax_mat.set_yticklabels([label_map[c] for c in cohorts], fontsize=6.5)
        ax_mat.set_ylim(-0.6, n_set - 0.4)
        ax_mat.set_xlim(-0.6, n_sub - 0.4)
        ax_mat.set_xticks([])
        ax_mat.set_xlabel("")
        ax_mat.grid(False)
        ax_bar.grid(False)

        ax_set.barh(range(n_set), [set_totals[c] for c in cohorts],
                    color=set_c, height=0.6)
        ax_set.bar_label(ax_set.containers[0], fmt="%d", fontsize=5.5,
                         padding=1, color="#444444")
        ax_set.set_xlim(0, max(set_totals.values()) * 1.25)
        ax_set.set_xlabel("n groups")
        ax_set.tick_params(axis="y", labelleft=False)
        ax_set.grid(False)

        theme.finalize(fig, "constrained")
        save(fig, cfg, pat.format(prefix=prefix, ptype=ptype,
                                  ptype_l=ptype.lower()))


