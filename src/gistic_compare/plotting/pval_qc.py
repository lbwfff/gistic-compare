"""Figures for module PQ (p-value QC: Q-Q plots + inflation)."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from ._core import *  # noqa: F401,F403


def plot_pq_qq(qq: list[dict], inflation=None, cfg=None, prefix: str = "PQ",
               per_fig: int = 6) -> None:
    """One panel per p-value set: observed vs uniform-expected -log10(p).

    Points hugging the diagonal = healthy; a systematic lift above it =
    inflation (lambda annotated per panel)."""
    if not qq:
        return
    style(cfg)
    lam_map = {}
    if inflation is not None and len(inflation):
        pooled = (inflation[inflation["comparison"] == "all"]
                  if "comparison" in inflation.columns else inflation)
        if pooled.empty:
            pooled = inflation
        for r in pooled.itertuples():
            lam_map[r.set_id] = r.lambda_infl
    # wrap sets into pages of `per_fig` panels
    for page in range(0, len(qq), per_fig):
        chunk = qq[page:page + per_fig]
        ncols = 3
        nrows = int(np.ceil(len(chunk) / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=theme.figsize(
            cfg, panel_h=2.6, total_w=7.0), squeeze=False)
        axes = axes.ravel()
        for k, entry in enumerate(chunk):
            ax = axes[k]
            exp, obs = entry["expected"], entry["observed"]
            ax.scatter(exp, obs, s=3, color="#4C72B0", alpha=0.5,
                       edgecolor="none")
            lim = max(float(obs.max()), float(exp.max()), 0.5) * 1.05
            ax.plot([0, lim], [0, lim], color="#999999", lw=0.7, ls="--")
            lam = lam_map.get(entry["set_id"], np.nan)
            ax.set_title(f"{entry['set_id']}  "
                         f"(\u03bb = {lam:.2f})"
                         if np.isfinite(lam) else entry["set_id"],
                         fontsize=6.2)
            ax.set_xlabel("expected \u2212log10 p", fontsize=6.5)
            ax.set_ylabel("observed \u2212log10 p", fontsize=6.5)
            ax.tick_params(labelsize=6)
            ax.grid(alpha=0.2)
        for k in range(len(chunk), len(axes)):
            axes[k].set_visible(False)
        fig.suptitle("Module PQ \u2014 p-value Q-Q (\u03bb = median(z\u00b2)"
                     "/0.4549; >1.1 = systematic bias)", fontsize=8)
        theme.finalize(fig, "tight")
        name = (f"{prefix}_qq" if len(qq) <= per_fig
                else f"{prefix}_qq_{page // per_fig + 1}")
        save(fig, cfg, name)


def figures_pq(res: dict, cfg) -> None:
    plot_pq_qq(res.get("qq"), res.get("inflation"), cfg)
