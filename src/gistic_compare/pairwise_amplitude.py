"""Module P7 (pairwise): copy-number *amplitude* differences among altered
samples.

P4/P5 ask whether the frequency of alteration differs between cohorts; GS
asks whether the G-score (frequency-weighted amplitude) differs.  None of
them answers a question that comes up constantly when two cohorts are put
side by side: **"this locus is altered at the same frequency in both
cohorts — but is it altered more strongly in one of them?"**  P7 splits the
per-locus signal into the two axes GISTIC itself uses:

  frequency  — fraction of samples altered in the unit
  amplitude  — the per-sample copy-change value, compared **among the
               altered samples only** (alteredness is defined per sample
               from that sample's own value, so the altered subset does not
               depend on cohort labels; a plain two-sample test on the
               subset is therefore exact — no extra permutation needed)

Units (three layers, all optional):
  window / band — 1 Mb windows and cytobands from GISTIC's own marker
                  matrix (``D.cap1.5.mat``, the data GISTIC scores on)
  peak          — union peaks (P1-style moderate matching) with the
                  per-sample ``- CN values`` rows of all_lesions.conf_90

GISTIC's own amplitude thresholds (0.1 / 0.9, as printed in the
``Amplitude Threshold`` column) are the defaults: altered = value beyond
+/- 0.1, high-level = value beyond +/- 0.9.

Per unit x direction x comparison the table reports: alteration frequency
(Fisher + BH), the mean amplitude among altered samples (Mann-Whitney +
BH, rank-biserial effect size), the high-level fraction among altered
samples (Fisher among altered + BH), and a four-way class combining the
frequency and amplitude axes (freq_driven / amp_driven / both / neither).

Outputs: P7_window_amplitude.tsv, P7_band_amplitude.tsv,
         P7_peak_amplitude.tsv, P7_summary.tsv
"""

from __future__ import annotations

import warnings
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from .config import out_table
from .marker_matrix import marker_window_matrix
from .pairwise_model import build_union_keys
from .pairwise_regional import (
    _fdr,
    _fisher_p,
    _load_band_genes,
    _regional_cfg,
)
from .score_track import make_grid


def _amp_cfg(cfg: dict[str, Any], rcfg: dict) -> dict:
    a = (cfg.get("analysis") or {}).get("amplitude") or {}
    return dict(
        amp_thr=float(a.get("amp_thr", 0.10)),
        del_thr=float(a.get("del_thr", -0.10)),
        high_amp=float(a.get("high_amp", 0.90)),
        high_del=float(a.get("high_del", -0.90)),
        fdr=float(a.get("fdr", rcfg["fdr"])),
    )


# ---------------------------------------------------------------------------
# shared per-unit statistics
# ---------------------------------------------------------------------------

def _rank_biserial(xa: np.ndarray, xb: np.ndarray) -> float:
    try:
        u1 = stats.mannwhitneyu(xa, xb, alternative="two-sided")[0]
    except ValueError:
        return np.nan
    return 2.0 * u1 / (len(xa) * len(xb)) - 1.0


def _direction_stats(va: np.ndarray, vb: np.ndarray, thr: float,
                     high: float, greater: bool, d: str) -> dict:
    """Frequency + altered-only amplitude stats for one unit/direction.

    ``greater`` selects the amplification direction (value >= thr) vs the
    deletion direction (value <= thr).  ``va``/``vb`` are the per-sample
    unit values of the two cohorts (NaN where a sample has no value).
    Column names carry the direction suffix ``d`` ('amp'/'del').
    """
    fin_a, fin_b = np.isfinite(va), np.isfinite(vb)
    na, nb = int(fin_a.sum()), int(fin_b.sum())
    if na == 0 or nb == 0:
        return {f"n_{d}_a": na, f"n_{d}_b": nb,
                f"freq_{d}_a": np.nan, f"freq_{d}_b": np.nan,
                f"p_freq_{d}": np.nan,
                f"n_alt_{d}_a": 0, f"n_alt_{d}_b": 0,
                f"amp_alt_{d}_a": np.nan, f"amp_alt_{d}_b": np.nan,
                f"delta_amp_{d}": np.nan,
                f"p_amp_{d}": np.nan, f"r_rb_{d}": np.nan,
                f"frac_high_{d}_a": np.nan, f"frac_high_{d}_b": np.nan,
                f"p_high_{d}": np.nan}
    if greater:
        alt_a, alt_b = va[fin_a] >= thr, vb[fin_b] >= thr
        hi_a, hi_b = va[fin_a][alt_a] >= high, vb[fin_b][alt_b] >= high
    else:
        alt_a, alt_b = va[fin_a] <= thr, vb[fin_b] <= thr
        hi_a, hi_b = va[fin_a][alt_a] <= high, vb[fin_b][alt_b] <= high
    ka, kb = int(alt_a.sum()), int(alt_b.sum())
    out = {
        f"n_{d}_a": na, f"n_{d}_b": nb,
        f"freq_{d}_a": ka / na, f"freq_{d}_b": kb / nb,
        f"p_freq_{d}": _fisher_p(ka, na, kb, nb),
        f"n_alt_{d}_a": ka, f"n_alt_{d}_b": kb,
    }
    if ka >= 2 and kb >= 2:
        xa, xb = va[fin_a][alt_a], vb[fin_b][alt_b]
        try:
            p_amp = float(stats.mannwhitneyu(xa, xb,
                                             alternative="two-sided")[1])
        except ValueError:
            p_amp = np.nan
        out.update({f"amp_alt_{d}_a": float(np.mean(xa)),
                    f"amp_alt_{d}_b": float(np.mean(xb)),
                    f"delta_amp_{d}": float(np.mean(xa) - np.mean(xb)),
                    f"p_amp_{d}": p_amp,
                    f"r_rb_{d}": _rank_biserial(xa, xb)})
    if ka >= 1 and kb >= 1:
        out.update({f"frac_high_{d}_a": float(np.mean(hi_a)),
                    f"frac_high_{d}_b": float(np.mean(hi_b)),
                    f"p_high_{d}": _fisher_p(int(hi_a.sum()), ka,
                                             int(hi_b.sum()), kb)})
    return out


def _classify(df: pd.DataFrame, fdr_thr: float) -> pd.DataFrame:
    """Four-way class per direction from the frequency/amplitude axes."""
    for d in ("amp", "del"):
        sf = df[f"fdr_freq_{d}"] < fdr_thr
        sa = df[f"fdr_amp_{d}"] < fdr_thr
        df[f"class_{d}"] = np.select(
            [sf & sa, sf & ~sa, ~sf & sa], ["both", "freq_driven",
                                            "amp_driven"], "neither")
    return df


def _add_fdr(df: pd.DataFrame) -> pd.DataFrame:
    for d in ("amp", "del"):
        for base, fdr in ((f"p_freq_{d}", f"fdr_freq_{d}"),
                          (f"p_amp_{d}", f"fdr_amp_{d}"),
                          (f"p_high_{d}", f"fdr_high_{d}")):
            if base in df.columns:
                df[fdr] = _fdr(df[base].to_numpy(float))
    df["delta_freq_amp"] = df["freq_amp_a"] - df["freq_amp_b"]
    df["delta_freq_del"] = df["freq_del_a"] - df["freq_del_b"]
    return df


# ---------------------------------------------------------------------------
# window / band layer (GISTIC marker matrix)
# ---------------------------------------------------------------------------

def _matrix_source(cfg: dict[str, Any], data: dict[str, dict],
                   cids: set[str], rcfg: dict) -> str:
    """One matrix source for all cohorts in the run: GISTIC's D matrix when
    every cohort has it, else seg-derived window values (P4 machinery).
    Mixing the two within a run would make cohorts incomparable."""
    cap = rcfg.get("marker_cap_file", "D.cap1.5.mat")
    have_d = all(
        (cfg["cohorts"][cid]["_abs"]
         / f"{cfg['cohorts'][cid]['prefix']}.{cap}").exists()
        for cid in cids)
    return "d_matrix" if have_d else "seg"


def _sample_ids_for(cfg: dict[str, Any], data: dict[str, dict], cid: str,
                    rcfg: dict, src: str) -> list[str]:
    """Sample order of the cohort's window matrix for the given source."""
    if src == "d_matrix":
        c = cfg["cohorts"][cid]
        cap = rcfg.get("marker_cap_file", "D.cap1.5.mat")
        from .marker_matrix import load_marker_matrix
        return list(load_marker_matrix(c["_abs"] / f"{c['prefix']}.{cap}")
                    ["samples"])
    return sorted(data[cid]["seg"]["sample"].unique())


def _matrix_for_cohort(cfg: dict[str, Any], data: dict[str, dict], cid: str,
                       grid: pd.DataFrame, rcfg: dict,
                       source: str) -> np.ndarray:
    """Windows x samples unit values from the chosen source."""
    if source == "d_matrix":
        c = cfg["cohorts"][cid]
        cap = rcfg.get("marker_cap_file", "D.cap1.5.mat")
        return marker_window_matrix(c["_abs"] / f"{c['prefix']}.{cap}",
                                    grid)[0]
    from .pairwise_regional import _window_values
    v, _ = _window_values(data[cid]["seg"], grid)
    return v


def _window_groups(grid: pd.DataFrame) -> list[tuple[str, str, int, int,
                                                 np.ndarray]]:
    out = []
    for i, r in enumerate(grid.itertuples()):
        out.append((f"{r.chrom}:{int(r.start/1e6)}-{int(r.end/1e6)}Mb",
                    r.chrom, int(r.start), int(r.end),
                    np.array([i], dtype=int)))
    return out


def _band_groups(grid: pd.DataFrame,
                 bands: pd.DataFrame) -> list[tuple[str, str, int, int,
                                                    np.ndarray]]:
    """Window indices per band (same overlap rule as P4/GS)."""
    gstart = grid["start"].to_numpy()
    gend = grid["end"].to_numpy()
    glen = (gend - gstart).astype(float)
    gchrom = grid["chrom"].to_numpy()
    out = []
    for b in bands.itertuples():
        band_len = float(b.end - b.start)
        ov = (np.minimum(gend, b.end)
              - np.maximum(gstart, b.start)).astype(float)
        thr = 0.5 * np.minimum(glen, band_len)
        idx = np.flatnonzero((gchrom == b.chrom) & (ov > thr))
        if len(idx):
            out.append((str(b.band), b.chrom, int(b.start), int(b.end), idx))
    return out


def _layer_stats(cfg: dict[str, Any], data: dict[str, dict], groups: list,
                 grid: pd.DataFrame, acfg: dict
                 ) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    """One table + per-unit per-sample value tables for a window/band layer."""
    from .gs_score import _unit_q  # per-unit cohort-q helper (display only)
    rows = []
    unit_names = []
    vals: dict[str, list[np.ndarray]] = {}
    for cmp_ in cfg["comparisons"]:
        ca, cb = cmp_["a"], cmp_["b"]
        if cmp_.get("kind") == "paired":
            print(f"  [P7] {ca}_vs_{cb}: paired comparison skipped — P7 "
                  "implements unpaired tests only (paired inference is on "
                  "the roadmap)")
            continue
        comp = f"{ca}_vs_{cb}"
        mat_a = data[ca]["_p7_matrix"]
        mat_b = data[cb]["_p7_matrix"]
        for name, chrom, start, end, idx in groups:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                va = np.nanmean(mat_a[idx], axis=0)
                vb = np.nanmean(mat_b[idx], axis=0)
            unit_names.append(name)
            vals.setdefault(ca, []).append(va)
            vals.setdefault(cb, []).append(vb)
            rec = dict(comparison=comp, unit=name, chrom=chrom,
                       start=int(start), end=int(end), n_windows=len(idx))
            for dirn in ("Amp", "Del"):
                greater = dirn == "Amp"
                thr = acfg["amp_thr"] if greater else acfg["del_thr"]
                high = acfg["high_amp"] if greater else acfg["high_del"]
                rec.update(_direction_stats(va, vb, thr, high, greater,
                                            dirn.lower()))
                rec[f"neglogq_{dirn.lower()}_a"] = _unit_q(
                    data[ca]["scores"], chrom, start, end, dirn)
                rec[f"neglogq_{dirn.lower()}_b"] = _unit_q(
                    data[cb]["scores"], chrom, start, end, dirn)
            rows.append(rec)
    df = _add_fdr(_unit_table(rows))
    df = _classify(df, acfg["fdr"])
    val_dfs = {}
    samples = {cid: list(data[cid]["_p7_sample_ids"]) for cid in vals}
    for cid, arrs in vals.items():
        val_dfs[cid] = pd.DataFrame(np.vstack(arrs), index=unit_names,
                                    columns=samples[cid])
    return df, val_dfs


def _unit_table(rows: list[dict]) -> pd.DataFrame:
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# peak layer (union peaks + all_lesions CN values rows)
# ---------------------------------------------------------------------------

def _peak_values(cfg: dict[str, Any], data: dict[str, dict], ca: str, cb: str,
                 keys: pd.DataFrame, acfg: dict) -> list[dict]:
    """Per-sample CN per union peak key from each cohort's own CN rows."""
    rows = []
    cn = {cid: data[cid]["peak_cn"] for cid in (ca, cb)}
    ksub = keys[keys["cohort"].isin([ca, cb])]
    for key, g in ksub.groupby("peak_key", sort=False):
        first = g.sort_values(["chrom", "start"]).iloc[0]
        rec = dict(peak_key=key, ptype=str(g["ptype"].iloc[0]),
                   cytoband=first["cytoband"], chrom=first["chrom"],
                   start=int(first["start"]), end=int(first["end"]))
        vals = {}
        for cid in (ca, cb):
            gc = g[g["cohort"] == cid]
            if gc.empty:
                vals[cid] = (None, False)
                continue
            sub = cn[cid][cn[cid]["peak_id"].isin(set(gc["peak_id"]))]
            v = (sub.groupby("sample")["cn"].mean() if not sub.empty
                 else pd.Series(dtype=float))
            vals[cid] = (v, True)
        rec["present_a"] = bool(vals[ca][1])
        rec["present_b"] = bool(vals[cb][1])
        for dirn in ("Amp", "Del"):
            greater = dirn == "Amp"
            thr = acfg["amp_thr"] if greater else acfg["del_thr"]
            high = acfg["high_amp"] if greater else acfg["high_del"]
            d = dirn.lower()
            if not (rec["present_a"] and rec["present_b"]):
                rec[f"present_{d}"] = ("a_only" if rec["present_a"] else
                                       "b_only" if rec["present_b"] else
                                       "none")
                continue
            rec[f"present_{d}"] = "both"
            rec.update(_direction_stats(vals[ca][0].to_numpy(float),
                                        vals[cb][0].to_numpy(float),
                                        thr, high, greater, d))
        rows.append(rec)
    return rows


# ---------------------------------------------------------------------------
# module entry
# ---------------------------------------------------------------------------

def run_module_p7(cfg: dict[str, Any], data: dict[str, dict]) -> dict[str, Any]:
    rcfg = _regional_cfg(cfg)
    acfg = _amp_cfg(cfg, rcfg)
    fdr_thr = acfg["fdr"]

    res: dict[str, Any] = dict(window=pd.DataFrame(), band=pd.DataFrame(),
                               peak=pd.DataFrame(), summary=pd.DataFrame(),
                               acfg=acfg)
    written = []
    unpaired_cohorts = {cid for cmp_ in cfg["comparisons"]
                        if cmp_.get("kind") != "paired"
                        for cid in (cmp_["a"], cmp_["b"])}

    # ---- window / band layers ----
    grid = make_grid(rcfg["bin_size"])
    res["matrices"] = {}
    if unpaired_cohorts:
        src = _matrix_source(cfg, data, unpaired_cohorts, rcfg)
        if src == "seg":
            print("  [P7] marker matrix (D.cap1.5.mat) missing for some "
                  "cohort — using seg-derived window values (P4 machinery)")
        for cid in unpaired_cohorts:
            data[cid]["_p7_matrix"] = _matrix_for_cohort(
                cfg, data, cid, grid, rcfg, src)
            data[cid]["_p7_sample_ids"] = _sample_ids_for(
                cfg, data, cid, rcfg, src)
        bands, _ = _load_band_genes(cfg, rcfg, data)
        win, win_vals = _layer_stats(cfg, data, _window_groups(grid), grid,
                                     acfg)
        for cid, vdf in win_vals.items():
            res["matrices"][f"{cid}_vals_window"] = vdf
        win.to_csv(out_table(cfg, "P7_window_amplitude.tsv"), sep="\t",
                   index=False)
        res["window"] = win
        written.append(f"window n={len(win)}")
        if not bands.empty:
            band, band_vals = _layer_stats(cfg, data, _band_groups(grid,
                                                                  bands),
                                           grid, acfg)
            for cid, vdf in band_vals.items():
                res["matrices"][f"{cid}_vals_band"] = vdf
            band.to_csv(out_table(cfg, "P7_band_amplitude.tsv"), sep="\t",
                        index=False)
            res["band"] = band
            written.append(f"band n={len(band)}")
        for cid in unpaired_cohorts:
            data[cid].pop("_p7_matrix", None)
            data[cid].pop("_p7_sample_ids", None)

    # ---- peak layer (union peaks + CN values rows) ----
    q_main = cfg["analysis"]["q_thresholds"][0]
    keys = build_union_keys(cfg, data, q_main)
    if not keys.empty:
        rows = []
        for cmp_ in cfg["comparisons"]:
            ca, cb = cmp_["a"], cmp_["b"]
            if cmp_.get("kind") == "paired":
                print(f"  [P7] {ca}_vs_{cb}: paired comparison skipped — P7 "
                      "implements unpaired tests only (paired inference is "
                      "on the roadmap)")
                continue
            comp = f"{ca}_vs_{cb}"
            for rec in _peak_values(cfg, data, ca, cb, keys, acfg):
                rec["comparison"] = comp
                rec["kind"] = cmp_.get("kind", "unpaired")
                rows.append(rec)
        peak = pd.DataFrame(rows)
        if not peak.empty:
            # single-side rows lack the stats columns; normalise the schema
            for d in ("amp", "del"):
                for col in (f"n_{d}_a", f"n_{d}_b", f"freq_{d}_a",
                            f"freq_{d}_b", f"p_freq_{d}", f"n_alt_{d}_a",
                            f"n_alt_{d}_b", f"amp_alt_{d}_a",
                            f"amp_alt_{d}_b", f"delta_amp_{d}", f"p_amp_{d}",
                            f"r_rb_{d}", f"frac_high_{d}_a",
                            f"frac_high_{d}_b", f"p_high_{d}"):
                    if col not in peak.columns:
                        peak[col] = np.nan
            peak = _classify(_add_fdr(peak), fdr_thr)
        peak.to_csv(out_table(cfg, "P7_peak_amplitude.tsv"), sep="\t",
                    index=False)
        res["peak"] = peak
        written.append(f"peak n={len(peak)}")

    # ---- summary ----
    srows = []
    for layer in ("window", "band", "peak"):
        df = res[layer]
        if df.empty:
            continue
        for comp, g in df.groupby("comparison"):
            row = dict(layer=layer, comparison=comp, n_units=len(g),
                       fdr_thr=fdr_thr, amp_thr=acfg["amp_thr"],
                       del_thr=acfg["del_thr"], high_amp=acfg["high_amp"],
                       high_del=acfg["high_del"])
            for d in ("amp", "del"):
                row[f"n_sig_freq_{d}"] = int((g[f"fdr_freq_{d}"] <
                                              fdr_thr).sum())
                row[f"n_sig_amp_{d}"] = int((g[f"fdr_amp_{d}"] <
                                             fdr_thr).sum())
                for c in ("both", "freq_driven", "amp_driven", "neither"):
                    row[f"n_class_{c}_{d}"] = int((g[f"class_{d}"] == c).sum())
            srows.append(row)
    summary = pd.DataFrame(srows)
    summary.to_csv(out_table(cfg, "P7_summary.tsv"), sep="\t", index=False)
    res["summary"] = summary

    print(f"  [P7] layers: {', '.join(written) if written else 'none'}"
          f"  (amp_thr={acfg['amp_thr']}, high_amp={acfg['high_amp']})")
    return res
