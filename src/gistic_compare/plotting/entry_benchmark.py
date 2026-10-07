"""figures_a..figures_g entry points.

(Split from the original single-module ``plotting``.)
"""

from __future__ import annotations

from ._core import *  # noqa: F401,F403
from .benchmark_a import (
    plot_jaccard,
    plot_matched_frequency,
    plot_permutation,
    plot_q_scatter,
    plot_three_way,
)
from .benchmark_b import plot_genome_tracks, plot_track_correlations
from .benchmark_c import plot_agreement, plot_frequency, plot_sample_agreement
from .benchmark_d import plot_driver_heatmap, plot_gene_jaccard
from .benchmark_e import plot_fga, plot_fga_vs_peaks
from .benchmark_f import plot_novelty_bars, plot_novelty_facet

# ---------------------------------------------------------------------------
# entry points
# ---------------------------------------------------------------------------

def figures_a(res: dict, cfg) -> None:
    plot_jaccard(res["jaccard"], cfg)
    plot_q_scatter(res["matched_pairs"], cfg, corr=res.get("correlations"))
    plot_matched_frequency(res["matched_pairs"], cfg, corr=res.get("correlations"))
    plot_permutation(res["permutation"], cfg,
                     null_samples=res.get("permutation_nulls"))
    plot_three_way(res["three_way"], cfg)


def figures_b(res: dict, cfg) -> None:
    plot_track_correlations(res["correlations"], cfg, method="pearson")
    plot_track_correlations(res["correlations"], cfg, method="spearman")
    plot_genome_tracks(res["grid"], cfg, ptype="Amp")
    plot_genome_tracks(res["grid"], cfg, ptype="Del")


def figures_c(res: dict, cfg) -> None:
    plot_agreement(res["peak_agreement"], cfg)
    plot_sample_agreement(res["sample_agreement"], cfg)
    plot_frequency(res["frequency"], cfg)


def figures_d(res: dict, cfg) -> None:
    plot_driver_heatmap(res["drivers"], cfg)
    plot_gene_jaccard(res["jaccard"], cfg)


def figures_e(res: dict, cfg) -> None:
    plot_fga(res["qc"], cfg, stats_df=res.get("stats"))
    plot_fga_vs_peaks(res["qc"], cfg, corr=res.get("fga_corr"))


def figures_f(res: dict, cfg) -> None:
    plot_novelty_facet(res["summary"], cfg)
    plot_novelty_bars(res["summary"], cfg)


