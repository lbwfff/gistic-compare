"""figures_g entry point.

(Split from the original single-module ``plotting``.)
"""

from __future__ import annotations

from ._core import *  # noqa: F401,F403
from .benchmark_g import (
    plot_paper_cytoband_heatmap,
    plot_paper_driver_heatmap,
    plot_paper_genome,
    plot_paper_top_freq,
)


def figures_g(res: dict, cfg) -> None:
    plot_paper_genome(res["genome"], cfg)
    plot_paper_cytoband_heatmap(res["cytoband"], cfg)
    plot_paper_top_freq(res["top_freq"], cfg)
    if res.get("driver") is not None and not res["driver"].empty:
        plot_paper_driver_heatmap(res["driver"], cfg)


