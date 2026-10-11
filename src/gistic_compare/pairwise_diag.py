"""Module PD (pairwise): baseline/batch diagnostic for the lambda inflation.

Motivated by module PQ: on mixed-source comparisons the MW and amp-side
p-value sets can show strong genomic inflation (lambda >> 1) even when the
event-count burden is balanced (module P9).  PD answers **why** and **what
survives**:

1. baseline check   — per-sample baseline = median ``seg_mean`` over the
   *neutral* segments (|seg_mean| < ``neutral_band``, default 0.25);
   cohort-level baseline distributions are compared, and the pooled
   ``seg_mean`` distributions are overlaid.  A shifted baseline is the
   classic batch signature.
2. shift structure  — the per-window median difference (CUE - CHB) along
   the genome, raw and after baseline removal.  A uniform offset = pure
   batch/plane effect; regional structure = biology or CN difference.
3. centering sensitivity (the core experiment) — the P4-style window
   tests (Fisher amp/del + MW) are run twice, on the raw and on the
   **within-sample centered** segments; lambda is reported for both and
   every window is assigned a trust layer:
     robust         significant in raw AND centered  -> quote it
     baseline_only  significant only before centering -> baseline-shift
                    driven candidate, downgrade
     centered_only  significant only after centering (rare; inspect)
     ns             not significant in either

Caveat (stats.md): within-sample centering also removes true whole-sample
scale differences (e.g. real aneuploidy shifts).  Read it together with
P9, which showed whether event-count burden differs.

Outputs: PD_lambda.tsv, PD_shift_regions.tsv, PD_summary.tsv
Figures: PD_baseline_dists, PD_shift_track, PD_lambda_bars
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from .config import out_table
from .pairwise_regional import (
    _counts,
    _fdr,
    _fisher_p,
    _regional_cfg,
    _window_values,
)
from .pval_qc import _lambda
from .score_track import make_grid


def _diag_cfg(cfg: dict[str, Any], rcfg: dict) -> dict:
    d = (cfg.get("analysis") or {}).get("diag") or {}
    return dict(
        neutral_band=float(d.get("neutral_band", 0.25)),
        smooth=int(d.get("smooth_windows", 20)),
        fdr=float(d.get("fdr", rcfg["fdr"])),
    )


def _baseline_per_sample(seg: pd.DataFrame, band: float) -> pd.DataFrame:
    """Per-sample baseline: median seg_mean over neutral segments."""
    neutral = seg[seg["seg_mean"].abs() < band]
    base = neutral.groupby("sample")["seg_mean"].median()
    return base.rename("baseline").reset_index()


def _center_seg(seg: pd.DataFrame, band: float) -> pd.DataFrame:
    base = _baseline_per_sample(seg, band).set_index("sample")["baseline"]
    out = seg.copy()
    out["seg_mean"] = out["seg_mean"] - out["sample"].map(base)
    return out


def _window_stats(vals_a: np.ndarray, vals_b: np.ndarray,
                  grid: pd.DataFrame, rcfg: dict) -> pd.DataFrame:
    """P4-style per-window tests (Fisher amp/del + MW) on one variant."""
    amp, dthr = rcfg["amp_thr"], rcfg["del_thr"]
    ka_amp, na = _counts(vals_a, amp, True)
    kb_amp, nb = _counts(vals_b, amp, True)
    ka_del, _ = _counts(vals_a, dthr, False)
    kb_del, _ = _counts(vals_b, dthr, False)
    rows = []
    for i in range(len(grid)):
        va, vb = vals_a[i], vals_b[i]
        ma, mb = np.isfinite(va), np.isfinite(vb)
        if ma.sum() < 3 or mb.sum() < 3:
            continue
        try:
            p_mw = float(stats.mannwhitneyu(va[ma], vb[mb],
                                            alternative="two-sided")[1])
        except ValueError:
            p_mw = np.nan
        rows.append(dict(
            win=int(i), chrom=grid["chrom"].iloc[i],
            start=int(grid["start"].iloc[i]), end=int(grid["end"].iloc[i]),
            freq_amp_a=ka_amp[i] / na[i] if na[i] else np.nan,
            freq_amp_b=kb_amp[i] / nb[i] if nb[i] else np.nan,
            freq_del_a=ka_del[i] / na[i] if na[i] else np.nan,
            freq_del_b=kb_del[i] / nb[i] if nb[i] else np.nan,
            p_amp=_fisher_p(int(ka_amp[i]), int(na[i]),
                            int(kb_amp[i]), int(nb[i])),
            p_del=_fisher_p(int(ka_del[i]), int(na[i]),
                            int(kb_del[i]), int(nb[i])),
            p_mw=p_mw,
            median_a=float(np.nanmedian(va)),
            median_b=float(np.nanmedian(vb)),
        ))
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    for c in ("p_amp", "p_del", "p_mw"):
        df[f"{c}_fdr"] = _fdr(df[c].to_numpy(float))
    df["delta_median"] = df["median_a"] - df["median_b"]
    return df


def _lambda_rows(stats_df: pd.DataFrame, variant: str, comp: str
                 ) -> list[dict]:
    rows = []
    for col, direction in (("p_amp", "Amp"), ("p_del", "Del"),
                           ("p_mw", "MW")):
        p = stats_df[col].to_numpy(float)
        p = p[np.isfinite(p)]
        if len(p) == 0:
            continue
        rows.append(dict(
            comparison=comp, variant=variant, set_col=col, direction=direction,
            n=len(p), n_p_lt_005=int((p < 0.05).sum()),
            median_p=round(float(np.median(p)), 5),
            lambda_infl=(round(_lambda(p), 3) if np.isfinite(_lambda(p))
                         else np.nan),
        ))
    return rows


def run_module_pd(cfg: dict[str, Any], data: dict[str, dict]) -> dict[str, Any]:
    rcfg = _regional_cfg(cfg)
    dcfg = _diag_cfg(cfg, rcfg)
    grid = make_grid(rcfg["bin_size"])

    # per-cohort window values, raw and centered
    raw_vals: dict[str, np.ndarray] = {}
    cen_vals: dict[str, np.ndarray] = {}
    baselines = []
    for cid, d in data.items():
        seg = d["seg"]
        b = _baseline_per_sample(seg, dcfg["neutral_band"])
        b.insert(0, "cohort", cid)
        baselines.append(b)
        raw_vals[cid], _ = _window_values(seg, grid)
        seg_c = _center_seg(seg, dcfg["neutral_band"])
        cen_vals[cid], _ = _window_values(seg_c, grid)
    baseline_df = pd.concat(baselines, ignore_index=True)

    lam_rows, layer_frames, sum_rows = [], [], []
    track_frames = []
    for cmp_ in cfg["comparisons"]:
        ca, cb = cmp_["a"], cmp_["b"]
        if cmp_.get("kind") == "paired":
            print(f"  [PD] skipped paired comparison {ca}_vs_{cb}")
            continue
        comp = f"{ca}_vs_{cb}"
        raw = _window_stats(raw_vals[ca], raw_vals[cb], grid, rcfg)
        cen = _window_stats(cen_vals[ca], cen_vals[cb], grid, rcfg)
        if raw.empty or cen.empty:
            continue
        for variant, st in (("raw", raw), ("centered", cen)):
            lam_rows.extend(_lambda_rows(st, variant, comp))

        # trust layers from the MW statistic (the inflated one)
        layer = pd.DataFrame(dict(
            win=raw["win"], chrom=raw["chrom"], start=raw["start"],
            end=raw["end"],
            delta_median_raw=raw["delta_median"],
            delta_median_cen=cen["delta_median"],
            p_mw_raw=raw["p_mw"], p_mw_fdr_raw=raw["p_mw_fdr"],
            p_mw_cen=cen["p_mw"], p_mw_fdr_cen=cen["p_mw_fdr"],
            p_amp_fdr_raw=raw["p_amp_fdr"], p_amp_fdr_cen=cen["p_amp_fdr"],
        ))
        sig_r = layer["p_mw_fdr_raw"] < dcfg["fdr"]
        sig_c = layer["p_mw_fdr_cen"] < dcfg["fdr"]
        layer["layer"] = np.select(
            [sig_r & sig_c, sig_r & ~sig_c, ~sig_r & sig_c],
            ["robust", "baseline_only", "centered_only"], "ns")
        layer.insert(0, "comparison", comp)
        layer_frames.append(layer)

        # smoothed shift tracks for the figure
        tr = layer[["comparison", "chrom", "start", "end",
                    "delta_median_raw", "delta_median_cen"]].copy()
        tr["raw_smooth"] = tr["delta_median_raw"].rolling(
            dcfg["smooth"], min_periods=3, center=True).mean()
        tr["cen_smooth"] = tr["delta_median_cen"].rolling(
            dcfg["smooth"], min_periods=3, center=True).mean()
        track_frames.append(tr)

        n_robust = int((layer["layer"] == "robust").sum())
        n_base = int((layer["layer"] == "baseline_only").sum())
        sum_rows.append(dict(
            comparison=comp, neutral_band=dcfg["neutral_band"],
            fdr=dcfg["fdr"], n_windows=len(layer),
            n_robust=n_robust, n_baseline_only=n_base,
            n_centered_only=int((layer["layer"] == "centered_only").sum()),
            median_baseline_a=round(float(baseline_df.loc[
                baseline_df["cohort"] == ca, "baseline"].median()), 4),
            median_baseline_b=round(float(baseline_df.loc[
                baseline_df["cohort"] == cb, "baseline"].median()), 4),
        ))
        lam_this = [r for r in lam_rows if r["comparison"] == comp]
        fmt = lambda v: (f"{v['lambda_infl']:.2f}"          # noqa: E731
                         if np.isfinite(v["lambda_infl"]) else "NA")
        print(f"  [PD] {comp}: lambda raw    "
              + " ".join(f"{r['set_col']}={fmt(r)}" for r in lam_this
                         if r["variant"] == "raw"))
        print(f"  [PD] {comp}: lambda centered "
              + " ".join(f"{r['set_col']}={fmt(r)}" for r in lam_this
                         if r["variant"] == "centered"))
        print(f"  [PD] {comp}: MW layers -> robust={n_robust} "
              f"baseline_only={n_base}")

    lam_df = pd.DataFrame(lam_rows)
    layers = (pd.concat(layer_frames, ignore_index=True) if layer_frames
              else pd.DataFrame())
    summary = pd.DataFrame(sum_rows)
    track = (pd.concat(track_frames, ignore_index=True) if track_frames
             else pd.DataFrame())

    lam_df.to_csv(out_table(cfg, "PD_lambda.tsv"), sep="\t", index=False)
    layers.to_csv(out_table(cfg, "PD_shift_regions.tsv"), sep="\t",
                  index=False)
    summary.to_csv(out_table(cfg, "PD_summary.tsv"), sep="\t", index=False)
    track.to_csv(out_table(cfg, "PD_shift_track.tsv"), sep="\t", index=False)
    baseline_df.to_csv(out_table(cfg, "PD_baselines.tsv"), sep="\t",
                       index=False)
    seg_values = {cid: d["seg"]["seg_mean"].clip(-1.5, 1.5).to_numpy(float)
                  for cid, d in data.items()}
    return dict(lambda_df=lam_df, layers=layers, summary=summary,
                track=track, baselines=baseline_df, seg_values=seg_values,
                dcfg=dcfg, skipped=False)
