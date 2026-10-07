"""Module P5 (pairwise workflow): regional differential on GISTIC's marker matrix.

This is the *parallel* version of P4 that, instead of resampling the post-CBS
`tumors.seg`, uses GISTIC's own per-marker / per-sample D matrix
(`<prefix>.D.cap1.5.mat`, a MATLAB v7.3 `SegArray`). Same window grid and the
same statistics as P4 (Fisher on amp/del tails + Mann-Whitney continuous +
cytoband tests + paper loci + random-split null), so results are directly
comparable with P4.

Outputs are prefixed `P5_`.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from .marker_matrix import marker_window_matrix
from .pairwise_regional import (
    _load_band_genes,
    _regional_cfg,
    analyze_windows,
    write_p4_like,
)
from .score_track import make_grid


def run_module_p5(cfg: dict[str, Any], data: dict[str, dict]) -> dict[str, Any]:
    rcfg = _regional_cfg(cfg)
    grid = make_grid(rcfg["bin_size"])
    cap_file = rcfg.get("marker_cap_file", "D.cap1.5.mat")

    vals: dict[str, np.ndarray] = {}
    for cid, d in data.items():
        c = cfg["cohorts"][cid]
        mat = c["_abs"] / f"{c['prefix']}.{cap_file}"
        if not mat.exists():
            raise FileNotFoundError(f"[P5] marker matrix not found: {mat}")
        v, samples = marker_window_matrix(mat, grid)
        vals[cid] = v
        print(f"  [P5] {cid}: {v.shape[0]} windows × {len(samples)} samples "
              f"from {mat.name}")

    bands, band_genes = _load_band_genes(cfg, rcfg, data)
    res = analyze_windows(cfg, vals, bands, band_genes, grid, rcfg, tag="P5")
    write_p4_like(cfg, res, "P5")
    return res
