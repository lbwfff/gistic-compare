"""Module GS (pairwise): G-score differential with a permutation test.

Unlike P4/P5 (which test the per-sample copy-number *distribution*), GS works
in GISTIC's own language. For every genomic unit (1 Mb window or cytoband) and
direction it forms the GISTIC-style amplitude of each sample, aggregates it
into a G-score, and asks whether the two cohorts' G-scores differ.

Amplitude is applied per window (mirroring how GISTIC scores per marker), then
pooled over the unit:

    amp window amplitude  a = min(cap, max(0,  v - amp_thr))
    del window amplitude  d = min(cap, max(0,  del_thr - v))
    sample unit score     A_i = mean over the unit's windows of a
    G_cohort              = mean over samples of A_i        (size-normalised)
    dG                    = G_A - G_B

`v` is the per-sample 1 Mb window value from GISTIC's own per-marker D matrix
(the data GISTIC scores on).

The p-value uses the exact permutation variance of dG (random cohort labels)
with a normal tail -- a permutation test in normal approximation. This avoids
the 1/(B+1) Monte-Carlo floor, which cannot survive BH over hundreds of units.
`analysis.gscore.method: mc` switches to Monte-Carlo (large B).

Each cohort's own GISTIC q at the unit is attached for reference only: the
difference of two cohort-specific q-values has no null and is NOT tested.

Outputs: GS_window_gscore.tsv, GS_band_gscore.tsv, GS_summary.tsv
"""

from __future__ import annotations

import warnings
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from .config import out_table
from .marker_matrix import marker_window_matrix
from .pairwise_regional import _load_band_genes, _locus_flag, _regional_cfg
from .score_track import make_grid


def _gscore_cfg(cfg: dict[str, Any], rcfg: dict) -> dict:
    g = (cfg.get("analysis") or {}).get("gscore") or {}
    return dict(
        cap=float(g.get("cap", 1.5)),
        perm=int(g.get("perm", 2000)),
        method=str(g.get("method", "normal")),   # normal | mc
        seed=int(g.get("seed", rcfg.get("seed", 42))),
    )


def _window_groups(grid: pd.DataFrame) -> list[tuple[str, str, int, int, np.ndarray]]:
    out = []
    for i, r in enumerate(grid.itertuples()):
        name = f"{r.chrom}:{int(r.start/1e6)}-{int(r.end/1e6)}Mb"
        out.append((name, r.chrom, int(r.start), int(r.end),
                    np.array([i], dtype=int)))
    return out


def _band_groups(grid: pd.DataFrame,
                 bands: pd.DataFrame) -> list[tuple[str, str, int, int, np.ndarray]]:
    """Window indices assigned to each band (overlap > half of the shorter)."""
    gstart = grid["start"].to_numpy()
    gend = grid["end"].to_numpy()
    glen = (gend - gstart).astype(float)
    gchrom = grid["chrom"].to_numpy()
    out = []
    for b in bands.itertuples():
        band_len = float(b.end - b.start)
        ov = (np.minimum(gend, b.end) - np.maximum(gstart, b.start)).astype(float)
        thr = 0.5 * np.minimum(glen, band_len)
        idx = np.flatnonzero((gchrom == b.chrom) & (ov > thr))
        if len(idx) == 0:
            continue
        out.append((str(b.band), b.chrom, int(b.start), int(b.end), idx))
    return out


def _perm_sd(amp: np.ndarray, nA: int, nB: int) -> float:
    """Exact permutation SD of dG = mean_A - mean_B under random labels.

    Var(dG) = s2_pop * n^2 / (nA * nB * (n - 1)), s2_pop the population
    variance (denominator n). Exact for sampling without replacement.
    """
    n = nA + nB
    if n < 3 or nA < 1 or nB < 1:
        return np.nan
    return float(np.sqrt(amp.var() * n * n / (nA * nB * (n - 1))))


def _dg_mc(amp: np.ndarray, nA: int, nB: int, obs: float,
           rng: np.random.Generator, B: int) -> float:
    """Monte-Carlo permutation p for |dG| (floor 1/(B+1); needs large B)."""
    n = nA + nB
    r = rng.random((B, n))
    sel = np.argpartition(r, nA, axis=1)[:, :nA]
    sum_a = amp[sel].sum(axis=1)
    dp = sum_a / nA - (amp.sum() - sum_a) / nB
    return (1.0 + float(np.sum(np.abs(dp) >= abs(obs)))) / (B + 1.0)


def _unit_q(scores: pd.DataFrame, chrom: str, start: int, end: int,
            ptype: str) -> float:
    """Max -log10(q) of the cohort's GISTIC regions (ptype) overlapping a unit."""
    if scores is None or scores.empty:
        return np.nan
    s = scores[(scores["type"] == ptype) & (scores["chrom"] == chrom)
               & (scores["start"] < end) & (scores["end"] > start)]
    if s.empty:
        return np.nan
    return float(s["neg_log10_q"].max())


def _fdr(p: np.ndarray) -> np.ndarray:
    from statsmodels.stats.multitest import multipletests
    out = np.full_like(p, np.nan, dtype=float)
    ok = np.isfinite(p)
    if ok.any():
        out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out


def _test_groups(cfg: dict[str, Any], data: dict[str, dict],
                 amp_w: dict[str, np.ndarray], del_w: dict[str, np.ndarray],
                 groups: list, gcfg: dict, fdr_thr: float,
                 rng: np.random.Generator) -> pd.DataFrame:
    cap, B = gcfg["cap"], gcfg["perm"]
    rows = []
    for cmp_ in cfg["comparisons"]:
        ca, cb = cmp_["a"], cmp_["b"]
        comp = f"{ca}_vs_{cb}"
        if cmp_.get("kind") == "paired":
            print(f"  [GS] skipped paired comparison {comp}")
            continue
        la = cfg.get("cohorts", {}).get(ca, {}).get("label", ca)
        lb = cfg.get("cohorts", {}).get(cb, {}).get("label", cb)
        sca, scb = data[ca].get("scores"), data[cb].get("scores")
        for name, chrom, start, end, idx in groups:
            rec = dict(comparison=comp, unit=name, chrom=chrom,
                       start=int(start), end=int(end), n_windows=len(idx))
            for dirn in ("Amp", "Del"):
                wa = amp_w[ca] if dirn == "Amp" else del_w[ca]
                wb = amp_w[cb] if dirn == "Amp" else del_w[cb]
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    xa = np.nanmean(wa[idx], axis=0)
                    xb = np.nanmean(wb[idx], axis=0)
                d = dirn.lower()
                nA = int(np.isfinite(xa).sum())
                nB = int(np.isfinite(xb).sum())
                if nA < 3 or nB < 3:
                    rec[f"p_{d}"] = np.nan
                    continue
                amp = np.concatenate([xa[np.isfinite(xa)],
                                      xb[np.isfinite(xb)]])
                GA = float(amp[:nA].mean())
                GB = float(amp[nA:].mean())
                dG = GA - GB
                sd = _perm_sd(amp, nA, nB)
                if gcfg["method"] == "mc":
                    p = _dg_mc(amp, nA, nB, dG, rng, B)
                else:
                    p = (float(2.0 * stats.norm.sf(abs(dG / sd)))
                         if np.isfinite(sd) and sd > 0 else np.nan)
                rec.update({
                    f"G_{d}_a": GA, f"G_{d}_b": GB, f"dG_{d}": dG,
                    f"sd_{d}": sd,
                    f"z_{d}": (dG / sd if sd and sd > 0 else np.nan),
                    f"p_{d}": p, f"n_{d}_a": nA, f"n_{d}_b": nB,
                    f"dir_{d}": (la if dG > 0 else lb),
                    f"neglogq_{d}_a": _unit_q(sca, chrom, start, end, dirn),
                    f"neglogq_{d}_b": _unit_q(scb, chrom, start, end, dirn),
                })
            rows.append(rec)
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    for d in ("amp", "del"):
        df[f"fdr_{d}"] = _fdr(df[f"p_{d}"].to_numpy(float))
        df[f"sig_{d}"] = df[f"fdr_{d}"] < fdr_thr
    return df


def _tag_paper(df: pd.DataFrame, rcfg: dict) -> pd.DataFrame:
    """Tag rows with preset (study-specific) locus / band names, if any."""
    loci = [n for n, *_ in rcfg.get("loci", [])]
    band_names = [n for n, _s in rcfg.get("band_marks", [])]
    lbl = []
    for r in df.itertuples():
        if str(r.unit) in band_names:
            lbl.append(str(r.unit))
        else:
            f = _locus_flag(rcfg.get("loci", []), r.chrom, r.start, r.end)
            lbl.append(f if f in loci else "")
    df["paper_locus"] = lbl
    return df


def run_module_gs(cfg: dict[str, Any], data: dict[str, dict]) -> dict[str, Any]:
    rcfg = _regional_cfg(cfg)
    gcfg = _gscore_cfg(cfg, rcfg)
    grid = make_grid(rcfg["bin_size"])
    cap = gcfg["cap"]
    amp_thr, del_thr = rcfg["amp_thr"], rcfg["del_thr"]
    fdr_thr = rcfg["fdr"]

    # per-sample window values from GISTIC's D matrix, then window-level
    # amplitudes (GISTIC applies the cap per marker).
    amp_w: dict[str, np.ndarray] = {}
    del_w: dict[str, np.ndarray] = {}
    cap_file = rcfg.get("marker_cap_file", "D.cap1.5.mat")
    for cid, c in cfg["cohorts"].items():
        mat = c["_abs"] / f"{c['prefix']}.{cap_file}"
        if not mat.exists():
            raise FileNotFoundError(f"[GS] marker matrix not found: {mat}")
        v, _ = marker_window_matrix(mat, grid)
        amp_w[cid] = np.clip(v - amp_thr, 0, cap)
        del_w[cid] = np.clip(del_thr - v, 0, cap)
        print(f"  [GS] {cid}: {v.shape[0]} windows x {v.shape[1]} samples "
              f"from {mat.name}")

    bands, _ = _load_band_genes(cfg, rcfg, data)
    rng = np.random.default_rng(gcfg["seed"])

    win = _tag_paper(_test_groups(cfg, data, amp_w, del_w,
                                  _window_groups(grid), gcfg, fdr_thr, rng),
                     rcfg)
    win.to_csv(out_table(cfg, "GS_window_gscore.tsv"), sep="\t", index=False)

    if bands.empty:
        print("  [GS] no cytobands; wrote window table only")
        return dict(window=win, band=pd.DataFrame(), summary=pd.DataFrame())

    band = _tag_paper(_test_groups(cfg, data, amp_w, del_w,
                                   _band_groups(grid, bands), gcfg, fdr_thr,
                                   rng), rcfg)
    band.to_csv(out_table(cfg, "GS_band_gscore.tsv"), sep="\t", index=False)

    summary = pd.DataFrame([
        dict(unit="window", n_units=len(win),
             n_sig_amp=int(win["sig_amp"].sum()),
             n_sig_del=int(win["sig_del"].sum()),
             min_fdr_amp=float(np.nanmin(win["fdr_amp"])),
             min_fdr_del=float(np.nanmin(win["fdr_del"]))),
        dict(unit="band", n_units=len(band),
             n_sig_amp=int(band["sig_amp"].sum()),
             n_sig_del=int(band["sig_del"].sum()),
             min_fdr_amp=float(np.nanmin(band["fdr_amp"])),
             min_fdr_del=float(np.nanmin(band["fdr_del"]))),
    ])
    summary.insert(0, "comparison",
                   win["comparison"].iloc[0] if len(win) else "")
    summary.insert(6, "cap", cap)
    summary.insert(7, "amp_thr", amp_thr)
    summary.insert(8, "del_thr", del_thr)
    summary.insert(9, "method", gcfg["method"])
    summary.to_csv(out_table(cfg, "GS_summary.tsv"), sep="\t", index=False)

    print(f"  [GS] window: sig amp={int(win['sig_amp'].sum())} "
          f"del={int(win['sig_del'].sum())}; "
          f"band: sig amp={int(band['sig_amp'].sum())} "
          f"del={int(band['sig_del'].sum())}  (method={gcfg['method']})")
    return dict(window=win, band=band, summary=summary, rcfg=rcfg)
