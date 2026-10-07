"""Module P4 (pairwise workflow): regional between-cohort differential analysis.

Answers the group-difference question directly:
    "Does the copy-number distribution differ between the two cohorts in this
     genomic region?"  (group difference, not recurrence)

P1–P3 are *detection / consistency* analyses on GISTIC peaks. P4 is independent
of GISTIC peak calling and works from per-sample tumors.seg.

Design (config `analysis.regional`):
  * unit  : fixed windows (default 1 Mb) over chr1–22
  * data  : per-sample segment means, length-weighted onto the window grid
  * calls : window value >= amp_thr -> amplified; <= del_thr -> deleted
  * per window × cohort: alteration frequency (frequency rule)
  * tests : two-sided Fisher on binary calls (amp/del) + Mann-Whitney on the
            continuous window value; BH-FDR across windows
  * regions: consecutive FDR-significant windows in the same direction are
            merged into differential regions (fewer, interpretable units)
  * exclusivity (paper wording): freq >= excl_freq_min in one cohort AND
            freq <= excl_freq_max in the other, per direction
  * null calibration: random within-cohort splits (equal size) run through the
            same binary pipeline -> false-positive background

Outputs: P4_regional_diff.tsv, P4_differential_regions.tsv,
         P4_exclusive_regions.tsv, P4_summary.tsv, P4_null_calibration.tsv
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from .config import CHROMS, out_table
from .presets import _resolve
from .presets import from_cfg as _preset_from_cfg
from .score_track import make_grid


def _mde_two_props(n_a: float, n_b: float, p_base: float,
                   alpha: float = 0.05, power: float = 0.8
                   ) -> tuple[float, float, float]:
    """MDE between two proportions from a fixed baseline frequency.

    Arcsine (Cohen's h) approximation: the test on ``2*asin(sqrt(p))`` has
    variance ``1/n_a + 1/n_b``, so an 80%-power (two-sided ``alpha``)
    contrast needs an arcsine half-gap of
    ``h/2 = (z_{1-a/2} + z_{power}) / 2 * sqrt(1/n_a + 1/n_b)``.
    With the baseline frequency held fixed we return the detectable
    upward shift, the detectable downward shift, and their maximum
    (conservative single-number MDE).
    """
    from math import asin, pi, sin, sqrt

    from scipy.stats import norm
    za = float(norm.ppf(1 - alpha / 2))
    zb = float(norm.ppf(power))
    half_gap = 0.5 * (za + zb) * sqrt(1.0 / n_a + 1.0 / n_b)
    phi = asin(sqrt(min(max(p_base, 1e-12), 1 - 1e-12)))
    p_up = sin(min(pi / 2, phi + half_gap)) ** 2
    p_dn = sin(max(0.0, phi - half_gap)) ** 2
    mde_up = max(0.0, p_up - p_base)
    mde_dn = max(0.0, p_base - p_dn)
    return mde_up, mde_dn, max(mde_up, mde_dn)


_POWER_BASELINES = (0.05, 0.10, 0.20, 0.30, 0.50)


def _mde_window(na: np.ndarray, nb: np.ndarray, p_bar: np.ndarray,
                alpha: float, power: float) -> np.ndarray:
    """Window-wise conservative MDE (vectorised over windows)."""
    from scipy.stats import norm
    za = float(norm.ppf(1 - alpha / 2))
    zb = float(norm.ppf(power))
    with np.errstate(divide="ignore", invalid="ignore"):
        half_gap = 0.5 * (za + zb) * np.sqrt(
            1.0 / np.maximum(na, 1) + 1.0 / np.maximum(nb, 1))
        phi = np.arcsin(np.sqrt(np.clip(p_bar, 0.0, 1.0)))
    p_up = np.sin(np.clip(phi + half_gap, 0, np.pi / 2)) ** 2
    p_dn = np.sin(np.clip(phi - half_gap, 0, np.pi / 2)) ** 2
    return np.maximum(p_up - p_bar, p_bar - p_dn)


def _power_summary(df: pd.DataFrame, rcfg: dict,
                   alpha: float = 0.05, power: float = 0.8) -> pd.DataFrame:
    """Per-comparison MDE table + observed-vs-MDE contrasts."""
    rows = []
    if df.empty:
        return pd.DataFrame()
    for comp, g in df.groupby("comparison"):
        na = float(np.nanmedian(g["n_a"]))
        nb = float(np.nanmedian(g["n_b"]))
        if not na or not nb:
            continue
        for direction, pcol in (("Amp", "freq_amp"), ("Del", "freq_del")):
            pa = g[f"{pcol}_a"].to_numpy(float)
            pb = g[f"{pcol}_b"].to_numpy(float)
            pooled = float(np.nanmean(np.concatenate([pa, pb])))
            for base in _POWER_BASELINES:
                mde_up, mde_dn, mde = _mde_two_props(na, nb, base, alpha,
                                                     power)
                rows.append(dict(
                    comparison=comp, direction=direction,
                    n_a=round(na), n_b=round(nb),
                    alpha=alpha, power_target=power,
                    baseline_freq=base,
                    mde_up=round(mde_up, 4), mde_down=round(mde_dn, 4),
                    mde=round(mde, 4),
                    observed_mean_pooled_freq=round(pooled, 4),
                ))
    return pd.DataFrame(rows)


def _window_mde_contrast(df: pd.DataFrame, rcfg: dict,
                         alpha: float = 0.05,
                         power: float = 0.8) -> dict[str, float]:
    """How many windows show a difference larger than the per-window MDE."""
    out: dict[str, float] = {}
    for direction in ("amp", "del"):
        na = df["n_a"].to_numpy(float)
        nb = df["n_b"].to_numpy(float)
        pa = df[f"freq_{direction}_a"].to_numpy(float)
        pb = df[f"freq_{direction}_b"].to_numpy(float)
        with np.errstate(invalid="ignore"):
            pooled = ((na * np.nan_to_num(pa) + nb * np.nan_to_num(pb))
                      / (na + nb))
            mde = _mde_window(na, nb, pooled, alpha, power)
            diff = np.abs(pa - pb)
        ok = np.isfinite(mde) & np.isfinite(diff)
        out[f"n_windows_diff_gt_mde_{direction}"] = int(
            (diff[ok] > mde[ok]).sum()) if ok.any() else 0
        out[f"max_abs_diff_{direction}"] = (round(float(np.nanmax(diff)), 4)
                                            if ok.any() else np.nan)
    return out


def _regional_cfg(cfg: dict[str, Any]) -> dict[str, Any]:
    r = (cfg.get("analysis") or {}).get("regional") or {}
    preset = _preset_from_cfg(cfg)
    # reference tables resolve relative to the config file's directory
    cfg_dir = (Path(cfg["_config_path"]).parent
               if cfg.get("_config_path") else Path.cwd())

    def _rp(key: str) -> str | None:
        v = r.get(key)
        return str(_resolve(v, cfg_dir)) if v else None

    return dict(
        bin_size=int(r.get("bin_size", 1_000_000)),
        amp_thr=float(r.get("amp_thr", 0.10)),
        del_thr=float(r.get("del_thr", -0.10)),
        excl_freq_min=float(r.get("excl_freq_min", 0.10)),
        excl_freq_max=float(r.get("excl_freq_max", 0.05)),
        fdr=float(r.get("fdr", 0.05)),
        null_splits=int(r.get("null_splits", 20)),
        seed=int(r.get("seed", 42)),
        cytoband_file=_rp("cytoband_file"),
        refgene_file=_rp("refgene_file"),
        marker_cap_file=r.get("marker_cap_file", "D.cap1.5.mat"),
        power_alpha=float(r.get("power_alpha", 0.05)),
        power_target=float(r.get("power_target", 0.8)),
        # preset (study-specific) loci / bookmarked bands, all empty by default
        loci=list(preset["loci"]),
        band_marks=list(preset["bands"]),
    )


def _window_values(seg: pd.DataFrame, grid: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
    """Length-weighted per-sample segment mean on the window grid."""
    samples = sorted(seg["sample"].unique())
    sidx = {s: i for i, s in enumerate(samples)}
    W, S = len(grid), len(samples)
    vsum = np.zeros((W, S))
    wsum = np.zeros((W, S))
    gstart = grid["start"].to_numpy()
    gend = grid["end"].to_numpy()
    gix = {c: np.flatnonzero((grid["chrom"] == c).to_numpy()) for c in CHROMS}
    for c in CHROMS:
        gi = gix[c]
        if len(gi) == 0:
            continue
        sub = seg[seg["chrom"] == c]
        if sub.empty:
            continue
        gs, ge = gstart[gi], gend[gi]
        for r in sub.itertuples():
            ov = np.minimum(ge, r.end) - np.maximum(gs, r.start)
            m = ov > 0
            if not m.any():
                continue
            j = sidx[r.sample]
            idx = gi[m]
            w = ov[m].astype(float)
            vsum[idx, j] += w * r.seg_mean
            wsum[idx, j] += w
    with np.errstate(invalid="ignore", divide="ignore"):
        vals = np.where(wsum > 0, vsum / wsum, np.nan)
    return vals, samples


def _counts(vals: np.ndarray, thr: float, greater: bool):
    finite = np.isfinite(vals)
    n = finite.sum(axis=1)
    if greater:
        k = ((vals >= thr) & finite).sum(axis=1)
    else:
        k = ((vals <= thr) & finite).sum(axis=1)
    return k.astype(int), n.astype(int)


def _fisher_p(ka: int, na: int, kb: int, nb: int) -> float:
    if na == 0 or nb == 0:
        return np.nan
    return float(stats.fisher_exact([[ka, na - ka], [kb, nb - kb]],
                                    alternative="two-sided")[1])


def _ztest_p_vec(ka, na, kb, nb) -> np.ndarray:
    with np.errstate(invalid="ignore", divide="ignore"):
        pa, pb = ka / na, kb / nb
        p = (ka + kb) / (na + nb)
        se = np.sqrt(p * (1 - p) * (1 / na + 1 / nb))
        return 2 * stats.norm.sf(np.abs((pa - pb) / se))


def _fdr(p: np.ndarray) -> np.ndarray:
    from statsmodels.stats.multitest import multipletests
    out = np.full_like(p, np.nan, dtype=float)
    ok = np.isfinite(p)
    if ok.any():
        out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out


def _locus_flag(loci, chrom: str, start: int, end: int) -> str:
    from .presets import loci_flag
    return loci_flag(loci, chrom, start, end)


def _differential(vals_a: np.ndarray, vals_b: np.ndarray,
                  grid: pd.DataFrame, rcfg: dict) -> pd.DataFrame:
    amp, dthr = rcfg["amp_thr"], rcfg["del_thr"]
    ka_amp, na = _counts(vals_a, amp, True)
    kb_amp, nb = _counts(vals_b, amp, True)
    ka_del, _ = _counts(vals_a, dthr, False)
    kb_del, _ = _counts(vals_b, dthr, False)

    rows = []
    for i in range(len(grid)):
        if na[i] == 0 or nb[i] == 0:
            continue
        fa_amp, fb_amp = ka_amp[i] / na[i], kb_amp[i] / nb[i]
        fa_del, fb_del = ka_del[i] / na[i], kb_del[i] / nb[i]
        va, vb = vals_a[i], vals_b[i]
        ma, mb = np.isfinite(va), np.isfinite(vb)
        if ma.sum() >= 3 and mb.sum() >= 3:
            try:
                p_mw = float(stats.mannwhitneyu(va[ma], vb[mb],
                                                alternative="two-sided")[1])
            except ValueError:
                p_mw = np.nan
        else:
            p_mw = np.nan
        rows.append(dict(
            win=int(i),
            chrom=grid["chrom"].iloc[i], start=int(grid["start"].iloc[i]),
            end=int(grid["end"].iloc[i]), n_a=int(na[i]), n_b=int(nb[i]),
            freq_amp_a=fa_amp, freq_amp_b=fb_amp,
            freq_del_a=fa_del, freq_del_b=fb_del,
            n_amp_a=int(ka_amp[i]), n_amp_b=int(kb_amp[i]),
            n_del_a=int(ka_del[i]), n_del_b=int(kb_del[i]),
            median_a=float(np.nanmedian(va)), median_b=float(np.nanmedian(vb)),
            p_amp=_fisher_p(int(ka_amp[i]), int(na[i]), int(kb_amp[i]), int(nb[i])),
            p_del=_fisher_p(int(ka_del[i]), int(na[i]), int(kb_del[i]), int(nb[i])),
            p_mw=p_mw,
        ))
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df["p_amp_fdr"] = _fdr(df["p_amp"].to_numpy(float))
    df["p_del_fdr"] = _fdr(df["p_del"].to_numpy(float))
    df["p_mw_fdr"] = _fdr(df["p_mw"].to_numpy(float))
    df["diff_amp"] = df["freq_amp_a"] - df["freq_amp_b"]
    df["diff_del"] = df["freq_del_a"] - df["freq_del_b"]
    df["paper_locus"] = [_locus_flag(rcfg["loci"], r.chrom, r.start, r.end)
                         for r in df.itertuples()]
    return df


def _region_tests(vals_a: np.ndarray, vals_b: np.ndarray,
                  df: pd.DataFrame, rcfg: dict,
                  raw_alpha: float = 0.10) -> pd.DataFrame:
    """Pool direction-consistent candidate windows into regions, then test once.

    Candidate region = maximal run of consecutive windows (same chrom) with
    MW p < raw_alpha and the same sign of median difference (A - B). For each
    region the per-sample window values are averaged across the region and a
    single between-cohort test is run (MW + Fisher on region-level calls).
    """
    amp, dthr = rcfg["amp_thr"], rcfg["del_thr"]
    s = df[df["p_mw"] < raw_alpha].copy()
    if s.empty:
        return pd.DataFrame()
    s["dir"] = np.sign(s["median_a"] - s["median_b"]).astype(float)
    blocks, cur = [], None
    for _, r in s.iterrows():
        if (cur is not None and r["chrom"] == cur["chrom"]
                and r["dir"] == cur["dir"] and r["win"] == cur["last"] + 1):
            cur["wins"].append(int(r["win"])); cur["last"] = int(r["win"])
            cur["end"] = int(r["end"])
        else:
            if cur is not None:
                blocks.append(cur)
            cur = dict(chrom=r["chrom"], start=int(r["start"]),
                       end=int(r["end"]), dir=float(r["dir"]),
                       wins=[int(r["win"])], last=int(r["win"]))
    if cur is not None:
        blocks.append(cur)

    rows = []
    for b in blocks:
        idx = np.array(b["wins"])
        ra = np.nanmean(vals_a[idx], axis=0)
        rb = np.nanmean(vals_b[idx], axis=0)
        ka, na = int(np.nansum(ra >= amp)), int(np.isfinite(ra).sum())
        kb, nb = int(np.nansum(rb >= amp)), int(np.isfinite(rb).sum())
        da = int(np.nansum(ra <= dthr))
        db = int(np.nansum(rb <= dthr))
        p_amp = _fisher_p(ka, na, kb, nb)
        p_del = _fisher_p(da, na, db, nb)
        try:
            p_mw = float(stats.mannwhitneyu(ra[np.isfinite(ra)],
                                            rb[np.isfinite(rb)],
                                            alternative="two-sided")[1])
        except ValueError:
            p_mw = np.nan
        rows.append(dict(
            chrom=b["chrom"], start=b["start"], end=b["end"],
            n_windows=len(idx),
            direction="A_higher" if b["dir"] > 0 else "B_higher",
            freq_amp_a=round(ka / na, 3) if na else np.nan,
            freq_amp_b=round(kb / nb, 3) if nb else np.nan,
            freq_del_a=round(da / na, 3) if na else np.nan,
            freq_del_b=round(db / nb, 3) if nb else np.nan,
            p_amp_region=p_amp, p_del_region=p_del, p_mw_region=p_mw,
            paper_loci=",".join(sorted({_locus_flag(
                rcfg["loci"], b["chrom"], b["start"], b["end"])} - {""})),
        ))
    reg = pd.DataFrame(rows)
    if reg.empty:
        return reg
    reg["p_min_region"] = reg[["p_amp_region", "p_del_region",
                               "p_mw_region"]].min(axis=1)
    reg["region_fdr"] = _fdr(reg["p_min_region"].to_numpy(float))
    return reg.sort_values("p_min_region").reset_index(drop=True)


def _null_calibration(seg: pd.DataFrame, grid: pd.DataFrame,
                      rcfg: dict) -> pd.DataFrame:
    amp, dthr = rcfg["amp_thr"], rcfg["del_thr"]
    vals, _ = _window_values(seg, grid)
    rng = np.random.default_rng(rcfg["seed"])
    n = vals.shape[1]
    half = n // 2
    rows = []
    for _ in range(rcfg["null_splits"]):
        perm = rng.permutation(n)
        A = vals[:, perm[:half]]; B = vals[:, perm[half:2 * half]]
        ka_amp, na = _counts(A, amp, True); kb_amp, nb = _counts(B, amp, True)
        ka_del, _ = _counts(A, dthr, False); kb_del, _ = _counts(B, dthr, False)
        pa = _ztest_p_vec(ka_amp, na, kb_amp, nb)
        pd_ = _ztest_p_vec(ka_del, na, kb_del, nb)
        fdr = _fdr(np.concatenate([pa, pd_]))
        rows.append(dict(
            n_fdr05=int(np.nansum(fdr < 0.05)),
            n_raw01=int(np.nansum(np.concatenate([pa, pd_]) < 0.01)),
            n_raw001=int(np.nansum(np.concatenate([pa, pd_]) < 0.001)),
        ))
    r = pd.DataFrame(rows)
    return pd.DataFrame([dict(
        n_samples=n, half=half, n_windows=len(grid),
        null_splits=rcfg["null_splits"],
        sig_windows_mean=float(r["n_fdr05"].mean()),
        sig_windows_p95=float(np.percentile(r["n_fdr05"], 95)),
        raw01_mean=float(r["n_raw01"].mean()),
        raw001_mean=float(r["n_raw001"].mean()),
    )])


def bundled_refs_dir() -> Path:
    """Directory with the tables shipped inside the package."""
    return Path(__file__).resolve().parent / "data" / "refs"


def _read_bundled_table(name: str) -> pd.DataFrame:
    p = bundled_refs_dir() / name
    if not p.exists():
        return pd.DataFrame()
    return pd.read_csv(p, sep="\t")


def _load_refgene(path) -> pd.DataFrame:
    """Read hg19 refgene (gene symbol + coordinates) from a GISTIC .mat file."""
    import scipy.io as sio
    p = Path(path)
    if not p.exists():
        return pd.DataFrame(columns=["chrom", "start", "end", "gene"])
    m = sio.loadmat(str(p), squeeze_me=True, struct_as_record=False)
    rg = m.get("rg")
    if rg is None:
        return pd.DataFrame(columns=["chrom", "start", "end", "gene"])
    rows = []
    for e in rg:
        d = {f: getattr(e, f) for f in e._fieldnames}
        chrom = str(d.get("chr", ""))
        if chrom not in CHROMS:
            continue
        sym = str(d.get("symb", "")).strip()
        if not sym:
            continue
        try:
            s, en = int(d["start"]), int(d["end"])
        except (TypeError, ValueError):
            continue
        rows.append((chrom, s, en, sym))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "gene"])


def _load_file_or_bundled(path, name: str,
                          mat_loader, columns: list[str]) -> pd.DataFrame:
    """Cytoband/refgene table: user file (.mat or .tsv), else bundled TSV."""
    if path is not None:
        p = Path(path)
        if p.exists():
            if p.suffix == ".mat":
                return mat_loader(p)
            return pd.read_csv(p, sep="\t")
        print(f"  [refs] warning: file not found, using bundled {name}: {p}")
    df = _read_bundled_table(name)
    return df if not df.empty else pd.DataFrame(columns=columns)


def _genes_by_band(refgene: pd.DataFrame, bands: pd.DataFrame,
                   which: str = "mid") -> dict[str, set]:
    """Assign each gene to a cytoband (by 'mid'/'start' position)."""
    out: dict[str, set] = {}
    if refgene.empty or bands.empty:
        return out
    for chrom, gg in refgene.groupby("chrom"):
        bb = bands[bands["chrom"] == chrom].sort_values("start")
        if bb.empty:
            continue
        starts = bb["start"].to_numpy()
        ends = bb["end"].to_numpy()
        names = bb["band"].to_numpy()
        pos = ((gg["start"].to_numpy() + gg["end"].to_numpy()) // 2
               if which == "mid" else gg["start"].to_numpy())
        for p, sym in zip(pos, gg["gene"].to_numpy()):
            idx = int(np.searchsorted(starts, p, side="right")) - 1
            if idx >= 0 and p < ends[idx]:
                out.setdefault(str(names[idx]), set()).add(str(sym))
    return out


def _load_cytobands(path) -> pd.DataFrame:
    import scipy.io as sio
    p = Path(path)
    if not p.exists():
        return pd.DataFrame()
    m = sio.loadmat(str(p), squeeze_me=True, struct_as_record=False)
    cyto = m.get("cyto")
    if cyto is None:
        return pd.DataFrame()
    rows = []
    for c in cyto:
        d = {f: getattr(c, f) for f in c._fieldnames}
        chrom = str(d["chr"])
        if not chrom.startswith("chr"):
            chrom = "chr" + chrom
        if chrom not in CHROMS:
            continue
        rows.append(dict(chrom=chrom, start=int(d["start"]), end=int(d["end"]),
                         band=str(d["name"])))
    return pd.DataFrame(rows)


def _pooled_stats(va: np.ndarray, vb: np.ndarray, idx: np.ndarray,
                  rcfg: dict) -> dict | None:
    """Pool per-sample window means over `idx`, then one between-cohort test."""
    amp, dthr = rcfg["amp_thr"], rcfg["del_thr"]
    if len(idx) == 0:
        return None
    ra = np.nanmean(va[idx], axis=0)
    rb = np.nanmean(vb[idx], axis=0)
    ra, rb = ra[np.isfinite(ra)], rb[np.isfinite(rb)]
    if len(ra) < 3 or len(rb) < 3:
        return None
    ka, kb = int((ra >= amp).sum()), int((rb >= amp).sum())
    da, db = int((ra <= dthr).sum()), int((rb <= dthr).sum())
    try:
        p_mw = float(stats.mannwhitneyu(ra, rb, alternative="two-sided")[1])
    except ValueError:
        p_mw = np.nan
    return dict(
        n_a=len(ra), n_b=len(rb),
        freq_amp_a=round(ka / len(ra), 3), freq_amp_b=round(kb / len(rb), 3),
        freq_del_a=round(da / len(ra), 3), freq_del_b=round(db / len(rb), 3),
        p_amp=_fisher_p(ka, len(ra), kb, len(rb)),
        p_del=_fisher_p(da, len(ra), db, len(rb)),
        p_mw=p_mw,
        median_a=float(np.median(ra)), median_b=float(np.median(rb)),
    )


def _band_tests(va: np.ndarray, vb: np.ndarray, grid: pd.DataFrame,
                bands: pd.DataFrame, rcfg: dict) -> pd.DataFrame:
    """Pre-specified cytoband-level differential (no data-driven selection).

    A window belongs to a band when their overlap exceeds half of the shorter
    of the two (so sub-window bands, e.g. 100 kb 19p13.13, are still tested).
    """
    if bands is None or bands.empty:
        return pd.DataFrame()
    gstart = grid["start"].to_numpy()
    gend = grid["end"].to_numpy()
    glen = (gend - gstart).astype(float)
    gchrom = grid["chrom"].to_numpy()
    rows = []
    for b in bands.itertuples():
        band_len = float(b.end - b.start)
        ov = (np.minimum(gend, b.end) - np.maximum(gstart, b.start)).astype(float)
        thr = 0.5 * np.minimum(glen, band_len)
        idx = np.flatnonzero((gchrom == b.chrom) & (ov > thr))
        st = _pooled_stats(va, vb, idx, rcfg)
        if st is None:
            continue
        st.update(chrom=b.chrom, start=int(b.start), end=int(b.end),
                  band=b.band, n_windows=int(len(idx)))
        rows.append(st)
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    for c in ("p_amp", "p_del", "p_mw"):
        df[c + "_fdr"] = _fdr(df[c].to_numpy(float))
    df["diff_amp"] = df["freq_amp_a"] - df["freq_amp_b"]
    df["diff_del"] = df["freq_del_a"] - df["freq_del_b"]
    return df.sort_values("p_mw").reset_index(drop=True)


def _paper_loci_summary(va: np.ndarray, vb: np.ndarray, grid: pd.DataFrame,
                        rcfg: dict) -> pd.DataFrame:
    if not rcfg["loci"]:
        return pd.DataFrame()
    mid = ((grid["start"] + grid["end"]) // 2).to_numpy()
    gchrom = grid["chrom"].to_numpy()
    rows = []
    for name, chrom, lo, hi, dirn, who in rcfg["loci"]:
        idx = np.flatnonzero((gchrom == chrom) & (mid >= lo) & (mid < hi))
        st = _pooled_stats(va, vb, idx, rcfg)
        if st is None:
            continue
        st.update(locus=name, paper_direction=dirn, expected_cohort=who,
                  chrom=chrom, start=lo, end=hi)
        rows.append(st)
    return pd.DataFrame(rows)


def _load_band_genes(cfg: dict[str, Any], rcfg: dict,
                     data: dict[str, dict]) -> tuple[pd.DataFrame, dict]:
    """Cytoband table + gene-sets per band (shared by P4 and P5).

    Cytoband/refgene tables come from the user's config when provided
    (`.mat` GISTIC assets or plain TSVs), otherwise from the small tables
    shipped inside the package (`gistic_compare/data/refs/`).
    """
    bands = _load_file_or_bundled(
        rcfg.get("cytoband_file"), "hg19_cytobands.tsv",
        _load_cytobands, ["chrom", "start", "end", "band"])
    band_genes: dict[str, set] = {}
    if not bands.empty:
        refgene = _load_file_or_bundled(
            rcfg.get("refgene_file"), "hg19_refgene.tsv",
            _load_refgene, ["chrom", "start", "end", "gene"])
        if not refgene.empty:
            band_genes = _genes_by_band(refgene, bands)
    if not band_genes:  # fallback: GISTIC per-peak gene lists
        for cid, d in data.items():
            gg = d.get("genes")
            if gg is None or gg.empty or "cytoband" not in gg.columns:
                continue
            for cb2, bg in gg.groupby("cytoband"):
                band_genes.setdefault(str(cb2).strip(), set()).update(
                    str(x).strip().strip("[]")
                    for x in bg["gene"] if str(x).strip())
    return bands, band_genes


def _null_from_vals(vals: np.ndarray, grid: pd.DataFrame,
                    rcfg: dict) -> pd.DataFrame:
    amp, dthr = rcfg["amp_thr"], rcfg["del_thr"]
    rng = np.random.default_rng(rcfg["seed"])
    n = vals.shape[1]
    half = n // 2
    rows = []
    for _ in range(rcfg["null_splits"]):
        perm = rng.permutation(n)
        A = vals[:, perm[:half]]
        B = vals[:, perm[half:2 * half]]
        ka_amp, na = _counts(A, amp, True)
        kb_amp, nb = _counts(B, amp, True)
        ka_del, _ = _counts(A, dthr, False)
        kb_del, _ = _counts(B, dthr, False)
        pa = _ztest_p_vec(ka_amp, na, kb_amp, nb)
        pdl = _ztest_p_vec(ka_del, na, kb_del, nb)
        fdr = _fdr(np.concatenate([pa, pdl]))
        both = np.concatenate([pa, pdl])
        rows.append(dict(n_fdr05=int(np.nansum(fdr < 0.05)),
                         n_raw01=int(np.nansum(both < 0.01)),
                         n_raw001=int(np.nansum(both < 0.001))))
    r = pd.DataFrame(rows)
    return pd.DataFrame([dict(
        n_samples=n, half=half, n_windows=len(grid),
        null_splits=rcfg["null_splits"],
        sig_windows_mean=float(r["n_fdr05"].mean()),
        sig_windows_p95=float(np.percentile(r["n_fdr05"], 95)),
        raw01_mean=float(r["n_raw01"].mean()),
        raw001_mean=float(r["n_raw001"].mean()),
    )])


def analyze_windows(cfg: dict[str, Any], vals: dict[str, np.ndarray],
                    bands: pd.DataFrame, band_genes: dict,
                    grid: pd.DataFrame, rcfg: dict,
                    tag: str = "P4") -> dict[str, Any]:
    """Shared analysis on a per-cohort (windows × samples) matrix."""
    fdr_thr = rcfg["fdr"]
    diff_frames, region_frames, band_frames, paper_frames = [], [], [], []
    exc_rows, summary_rows = [], []
    for cmp_ in cfg["comparisons"]:
        ca, cb = cmp_["a"], cmp_["b"]
        comp = f"{ca}_vs_{cb}"
        kind = cmp_.get("kind", "unpaired")
        if kind == "paired":
            print(f"  [{tag}] skipped paired comparison {comp} "
                  f"(needs a paired design)")
            continue
        df = _differential(vals[ca], vals[cb], grid, rcfg)
        if df.empty:
            continue
        df.insert(0, "comparison", comp)
        df.insert(1, "kind", kind)
        diff_frames.append(df)

        reg = _region_tests(vals[ca], vals[cb], df, rcfg)
        if not reg.empty:
            reg.insert(0, "comparison", comp)
            reg.insert(1, "kind", kind)
            region_frames.append(reg)

        if not bands.empty:
            bt = _band_tests(vals[ca], vals[cb], grid, bands, rcfg)
            if not bt.empty:
                bt.insert(0, "comparison", comp)
                bt.insert(1, "kind", kind)
                band_frames.append(bt)
        pt = _paper_loci_summary(vals[ca], vals[cb], grid, rcfg)
        if not pt.empty:
            pt.insert(0, "comparison", comp)
            pt.insert(1, "kind", kind)
            paper_frames.append(pt)

        for dirn in ("Amp", "Del"):
            fa_all = df["freq_amp_a"] if dirn == "Amp" else df["freq_del_a"]
            fb_all = df["freq_amp_b"] if dirn == "Amp" else df["freq_del_b"]
            pv_all = df["p_amp"] if dirn == "Amp" else df["p_del"]
            fdr_all = df["p_amp_fdr"] if dirn == "Amp" else df["p_del_fdr"]
            for side in ("a", "b"):
                if side == "a":
                    freq_ok = (fa_all >= rcfg["excl_freq_min"]) & \
                              (fb_all <= rcfg["excl_freq_max"])
                else:
                    freq_ok = (fb_all >= rcfg["excl_freq_min"]) & \
                              (fa_all <= rcfg["excl_freq_max"])
                sig_ok = (fdr_all < fdr_thr) | (df["p_mw_fdr"] < fdr_thr)
                sel = df[freq_ok | sig_ok]
                for _, r in sel.iterrows():
                    fa = r["freq_amp_a"] if dirn == "Amp" else r["freq_del_a"]
                    fb = r["freq_amp_b"] if dirn == "Amp" else r["freq_del_b"]
                    exc_rows.append(dict(
                        comparison=comp, kind=kind, direction=dirn, side=side,
                        chrom=r["chrom"], start=int(r["start"]),
                        end=int(r["end"]), paper_locus=r["paper_locus"],
                        freq_a=round(float(fa), 3), freq_b=round(float(fb), 3),
                        fisher_p=(float(pv_all.loc[r.name])),
                        fdr=(float(fdr_all.loc[r.name])),
                        mw_fdr=(float(r["p_mw_fdr"])),
                        freq_exclusive=bool(freq_ok.loc[r.name]),
                        significant=bool((fdr_all.loc[r.name] < fdr_thr)
                                         or (r["p_mw_fdr"] < fdr_thr)),
                    ))

        sig_win = int(((df["p_amp_fdr"] < fdr_thr) | (df["p_del_fdr"] < fdr_thr)
                       | (df["p_mw_fdr"] < fdr_thr)).sum())
        freq_exc = int((((df["freq_amp_a"] >= rcfg["excl_freq_min"]) & (df["freq_amp_b"] <= rcfg["excl_freq_max"]))
                        | ((df["freq_amp_b"] >= rcfg["excl_freq_min"]) & (df["freq_amp_a"] <= rcfg["excl_freq_max"]))
                        | ((df["freq_del_a"] >= rcfg["excl_freq_min"]) & (df["freq_del_b"] <= rcfg["excl_freq_max"]))
                        | ((df["freq_del_b"] >= rcfg["excl_freq_min"]) & (df["freq_del_a"] <= rcfg["excl_freq_max"]))).sum())
        paper = df[df["paper_locus"] != ""]
        mde = _window_mde_contrast(df, rcfg,
                                   alpha=rcfg.get("power_alpha", 0.05),
                                   power=rcfg.get("power_target", 0.8))
        summary_rows.append(dict(
            comparison=comp, kind=kind, bin_size=rcfg["bin_size"],
            amp_thr=rcfg["amp_thr"], del_thr=rcfg["del_thr"], fdr_thr=fdr_thr,
            n_windows=len(df), n_sig_windows=sig_win,
            n_regions=(len(reg) if not reg.empty else 0),
            n_regions_fdr05=(int((reg["region_fdr"] < fdr_thr).sum())
                             if not reg.empty else 0),
            n_freq_exclusive_windows=freq_exc, n_paper_windows=len(paper),
            paper_loci_hit=",".join(sorted(set(paper["paper_locus"]))),
            max_abs_diff_amp=mde["max_abs_diff_amp"],
            max_abs_diff_del=mde["max_abs_diff_del"],
            n_windows_diff_gt_mde_amp=mde["n_windows_diff_gt_mde_amp"],
            n_windows_diff_gt_mde_del=mde["n_windows_diff_gt_mde_del"],
        ))

    regional = (pd.concat(diff_frames, ignore_index=True) if diff_frames
                else pd.DataFrame())
    regions = (pd.concat(region_frames, ignore_index=True) if region_frames
               else pd.DataFrame())
    band_df = (pd.concat(band_frames, ignore_index=True) if band_frames
               else pd.DataFrame())
    sig_rows = []
    for comp in (list(dict.fromkeys(band_df["comparison"]))
                 if not band_df.empty else []):
        ca, cb = (comp.split("_vs_") + ["A", "B"])[:2]
        la = cfg.get("cohorts", {}).get(ca, {}).get("label", ca)
        lb = cfg.get("cohorts", {}).get(cb, {}).get("label", cb)
        for _, r in band_df[band_df["comparison"] == comp].iterrows():
            sig_amp = np.isfinite(r["p_amp"]) and r["p_amp"] < 0.05
            sig_del = np.isfinite(r["p_del"]) and r["p_del"] < 0.05
            sig_mw = np.isfinite(r["p_mw_fdr"]) and r["p_mw_fdr"] < 0.05
            if not (sig_amp or sig_del or sig_mw):
                continue
            genes = sorted(band_genes.get(str(r["band"]).strip(), set()))
            sig_rows.append(dict(
                comparison=comp, band=r["band"], chrom=r["chrom"],
                start=int(r["start"]), end=int(r["end"]),
                sig_amp=sig_amp, sig_del=sig_del, sig_mw_fdr=sig_mw,
                direction_amp=(la if r["diff_amp"] > 0 else lb),
                direction_del=(la if r["diff_del"] > 0 else lb),
                freq_amp_a=round(r["freq_amp_a"], 3),
                freq_amp_b=round(r["freq_amp_b"], 3), diff_amp=r["diff_amp"],
                freq_del_a=round(r["freq_del_a"], 3),
                freq_del_b=round(r["freq_del_b"], 3), diff_del=r["diff_del"],
                p_amp=r["p_amp"], p_del=r["p_del"], p_mw=r["p_mw"],
                p_mw_fdr=r["p_mw_fdr"], n_genes=len(genes),
                genes=",".join(genes),
            ))
    sig_df = pd.DataFrame(sig_rows)
    paper_df = (pd.concat(paper_frames, ignore_index=True) if paper_frames
                else pd.DataFrame())
    exclusive = pd.DataFrame(exc_rows)
    summary = pd.DataFrame(summary_rows)
    power = _power_summary(regional, rcfg,
                           alpha=rcfg.get("power_alpha", 0.05),
                           power=rcfg.get("power_target", 0.8))

    null_rows = []
    for cid, v in vals.items():
        nc = _null_from_vals(v, grid, rcfg)
        nc.insert(0, "cohort", cid)
        null_rows.append(nc)
    null = pd.concat(null_rows, ignore_index=True) if null_rows else pd.DataFrame()

    return dict(regional=regional, regions=regions, bands=band_df,
                significant=sig_df, paper=paper_df, exclusive=exclusive,
                summary=summary, power=power, null=null, grid=grid, rcfg=rcfg,
                band_genes=band_genes)


def write_p4_like(cfg: dict[str, Any], res: dict, prefix: str) -> None:
    for name, fn in (("regional_diff", "regional"),
                     ("differential_regions", "regions"),
                     ("cytoband_diff", "bands"),
                     ("significant_bands", "significant"),
                     ("paper_loci", "paper"),
                     ("exclusive_regions", "exclusive"),
                     ("summary", "summary"),
                     ("power", "power"),
                     ("null_calibration", "null")):
        res[fn].to_csv(out_table(cfg, f"{prefix}_{name}.tsv"), sep="\t",
                       index=False)


def run_module_p4(cfg: dict[str, Any], data: dict[str, dict]) -> dict[str, Any]:
    rcfg = _regional_cfg(cfg)
    grid = make_grid(rcfg["bin_size"])
    vals: dict[str, np.ndarray] = {}
    for cid, d in data.items():
        vals[cid], _ = _window_values(d["seg"], grid)
    bands, band_genes = _load_band_genes(cfg, rcfg, data)
    res = analyze_windows(cfg, vals, bands, band_genes, grid, rcfg, tag="P4")
    write_p4_like(cfg, res, "P4")
    return res
