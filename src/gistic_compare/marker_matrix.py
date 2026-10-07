"""Read GISTIC's marker-level D matrix (MATLAB v7.3 SegArray) as window values.

The `D.cap1.5.mat` file stores `D.dat` as a MATLAB `SegArray` (also produced by
GISTIC's own code): for each sample a list of breakpoints over the marker axis
(`#refs#/p` -> ir/jc) plus the per-segment value (`#refs#/o`). Marker
coordinates are `D/pos` (bp) and `D/chrn`; sample ids are the chars in
`D/sdesc`.  This is the per-marker, per-sample copy value GISTIC scores on.

We resample it onto the same fixed window grid used by module P4, returning a
`windows × samples` matrix so the P4 statistics can be reused unchanged.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .config import CHROMS


def _decode_sdesc(sd: np.ndarray) -> list[str]:
    return ["".join(chr(int(c)) for c in sd[:, j] if int(c) > 0).strip()
            for j in range(sd.shape[1])]


def load_marker_matrix(mat_path) -> dict:
    """Return dict(pos, chrom, samples, ir, jc, segval) from D.cap1.5.mat."""
    import h5py
    p = Path(mat_path)
    with h5py.File(p, "r") as h:
        pos = h["D/pos"][:].ravel().astype(np.int64)
        chrn = h["D/chrn"][:].ravel()
        sd = h["D/sdesc"][:]
        samples = _decode_sdesc(sd)
        ir = h["#refs#/p/ir"][:].astype(np.int64)
        jc = h["#refs#/p/jc"][:].astype(np.int64)
        segval = h["#refs#/o"][:].ravel().astype(float)
    chrom = np.array(["chr" + str(int(c)) for c in chrn], dtype=object)
    return dict(pos=pos, chrom=chrom, samples=samples, ir=ir, jc=jc,
                segval=segval)


def marker_window_matrix(mat_path, grid: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
    """Resample the markers × samples matrix onto `grid` (windows × samples)."""
    mm = load_marker_matrix(mat_path)
    pos, chrom, samples = mm["pos"], mm["chrom"], mm["samples"]
    M = len(pos)
    W = len(grid)

    # marker -> window index (by nearest / overlapping bin; use midpoint binning)
    win_of = np.full(M, -1, dtype=np.int64)
    grid_by_chrom = {c: grid.index[(grid["chrom"] == c)].to_numpy() for c in CHROMS}
    for c in CHROMS:
        gi = grid_by_chrom[c]
        if len(gi) == 0:
            continue
        starts = grid.loc[gi, "start"].to_numpy()
        ends = grid.loc[gi, "end"].to_numpy()
        sel = np.flatnonzero(chrom == c)
        if len(sel) == 0:
            continue
        ps = pos[sel]
        idx = np.searchsorted(starts, ps, side="right") - 1
        ok = (idx >= 0) & (ps < ends[np.clip(idx, 0, len(ends) - 1)])
        win_of[sel[ok]] = gi[idx[ok]]

    cnt = np.bincount(win_of[win_of >= 0], minlength=W).astype(float)
    out = np.full((W, len(samples)), np.nan)
    ir, jc, segval = mm["ir"], mm["jc"], mm["segval"]
    for j in range(len(samples)):
        a, b = int(jc[j]), int(jc[j + 1])
        if b <= a:
            continue
        rows = ir[a:b]
        vv = segval[a:b]
        markerval = np.empty(M)
        for k in range(len(rows)):
            r0 = rows[k]
            r1 = rows[k + 1] if k + 1 < len(rows) else M
            markerval[r0:r1] = vv[k]
        w = win_of
        s = np.bincount(w[w >= 0], weights=markerval[w >= 0], minlength=W)
        out[:, j] = np.where(cnt > 0, s / cnt, np.nan)
    return out, samples
