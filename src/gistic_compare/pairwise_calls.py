"""Module P2 (pairwise workflow): per-sample call concordance for paired pairs.

Reuses call-stat helpers from sample_calls; only comparisons with
kind == 'paired' and a non-empty shared sample set are analysed.
Unpaired pairs are skipped (tables keep headers for parseability).
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from . import intervals as iv
from .config import out_table
from .sample_calls import (
    _HAS_STATSMODELS,
    _exact_mcnemar_p,
    _pivot_calls,
    _weighted_kappa,
)

_COLS = {
    "peak_agreement": ["comparison", "ptype", "a_peak_id", "b_peak_id",
                       "a_cytoband", "b_cytoband", "n_samples",
                       "exact_agreement", "weighted_kappa", "a_freq", "b_freq"],
    "sample_agreement": ["comparison", "ptype", "sample", "n_peaks",
                         "exact_agreement", "weighted_kappa"],
    "mcnemar": ["comparison", "ptype", "a_peak_id", "b_peak_id",
                "a_cytoband", "b_cytoband", "only_a_altered",
                "only_b_altered", "mcnemar_p", "mcnemar_fdr"],
    "summary": ["comparison", "kind", "ptype", "n_matched_peaks",
                "n_paired_samples", "mean_exact_agreement",
                "median_exact_agreement", "mean_weighted_kappa",
                "median_weighted_kappa"],
}


def _fdr_bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    out = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    if _HAS_STATSMODELS:
        from statsmodels.stats.multitest import multipletests
        _, fdr, _, _ = multipletests(p[ok], method="fdr_bh")
        out[ok] = fdr
        return out
    # fallback BH
    order = np.argsort(p[ok])
    ranked = p[ok][order]
    fdr_sorted = ranked * len(ranked) / (np.arange(len(ranked)) + 1)
    fdr_sorted = np.minimum.accumulate(fdr_sorted[::-1])[::-1]
    fdr_sorted = np.clip(fdr_sorted, 0, 1)
    fdr = np.empty_like(fdr_sorted)
    fdr[order] = fdr_sorted
    out[ok] = fdr
    return out


def run_module_p2(cfg: dict[str, Any], data: dict[str, dict]) -> dict[str, Any]:
    ana = cfg["analysis"]
    q_main = ana["q_thresholds"][0]
    ov = ana["overlap"]

    peak_rows: list[dict[str, Any]] = []
    sample_rows: list[dict[str, Any]] = []
    mcnemar_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    skipped: list[str] = []

    for cmp_ in cfg["comparisons"]:
        ca, cb = cmp_["a"], cmp_["b"]
        kind = cmp_.get("kind", "unknown")
        comp = f"{ca}_vs_{cb}"
        if kind != "paired":
            skipped.append(f"{comp}({kind})")
            continue

        calls_a = _pivot_calls(data[ca]["calls"])
        calls_b = _pivot_calls(data[cb]["calls"])
        common = sorted(set(calls_a.columns) & set(calls_b.columns))
        if not common:
            skipped.append(f"{comp}(no shared samples)")
            continue
        calls_a = calls_a[common]
        calls_b = calls_b[common]

        for ptype in ("Amp", "Del"):
            fa = iv.filter_peaks(data[ca]["peaks"], q_main, ptype)
            fb = iv.filter_peaks(data[cb]["peaks"], q_main, ptype)
            res = iv.match_peaks(fa, fb, min_bp=ov["min_bp"],
                                 reciprocal=ov["reciprocal_shorter"],
                                 mode="either")
            n_pairs_pt = 0
            kappas: list[float] = []
            accs: list[float] = []
            for _, pr in res.pairs.iterrows():
                ida, idb = pr["a_peak_id"], pr["b_peak_id"]
                if ida not in calls_a.index or idb not in calls_b.index:
                    continue
                va = calls_a.loc[ida].to_numpy(dtype=float)
                vb = calls_b.loc[idb].to_numpy(dtype=float)
                ok = np.isfinite(va) & np.isfinite(vb)
                if ok.sum() == 0:
                    continue
                acc = float((va[ok] == vb[ok]).mean())
                kappa = _weighted_kappa(va, vb)
                n_pairs_pt += 1
                accs.append(acc)
                kappas.append(kappa)
                peak_rows.append(dict(
                    comparison=comp, ptype=ptype,
                    a_peak_id=ida, b_peak_id=idb,
                    a_cytoband=pr.get("a_cytoband", ""),
                    b_cytoband=pr.get("b_cytoband", ""),
                    n_samples=int(ok.sum()),
                    exact_agreement=acc, weighted_kappa=kappa,
                    a_freq=float(np.nanmean(va[ok] > 0)),
                    b_freq=float(np.nanmean(vb[ok] > 0)),
                ))

                a_bin = (va[ok] > 0).astype(int)
                b_bin = (vb[ok] > 0).astype(int)
                b01 = int(((a_bin == 0) & (b_bin == 1)).sum())
                b10 = int(((a_bin == 1) & (b_bin == 0)).sum())
                pval = _exact_mcnemar_p(b10, b01)
                mcnemar_rows.append(dict(
                    comparison=comp, ptype=ptype,
                    a_peak_id=ida, b_peak_id=idb,
                    a_cytoband=pr.get("a_cytoband", ""),
                    b_cytoband=pr.get("b_cytoband", ""),
                    only_a_altered=b10, only_b_altered=b01,
                    mcnemar_p=pval,
                ))

            # per-sample agreement across this ptype's matched pairs
            pairs_for_pt = [r for r in peak_rows
                            if r["comparison"] == comp and r["ptype"] == ptype]
            rows_a, rows_b = [], []
            for r in pairs_for_pt:
                if (r["a_peak_id"] in calls_a.index
                        and r["b_peak_id"] in calls_b.index):
                    rows_a.append(calls_a.loc[r["a_peak_id"]].to_numpy(float))
                    rows_b.append(calls_b.loc[r["b_peak_id"]].to_numpy(float))
            if rows_a:
                A = np.vstack(rows_a)
                B = np.vstack(rows_b)
                for si, s in enumerate(common):
                    va, vb = A[:, si], B[:, si]
                    m = np.isfinite(va) & np.isfinite(vb)
                    if m.sum() == 0:
                        continue
                    sample_rows.append(dict(
                        comparison=comp, ptype=ptype, sample=s,
                        n_peaks=int(m.sum()),
                        exact_agreement=float((va[m] == vb[m]).mean()),
                        weighted_kappa=_weighted_kappa(va, vb),
                    ))

            if n_pairs_pt:
                summary_rows.append(dict(
                    comparison=comp, kind=kind, ptype=ptype,
                    n_matched_peaks=n_pairs_pt,
                    n_paired_samples=len(common),
                    mean_exact_agreement=float(np.mean(accs)),
                    median_exact_agreement=float(np.median(accs)),
                    mean_weighted_kappa=float(np.mean(kappas)),
                    median_weighted_kappa=float(np.median(kappas)),
                ))

    mc = pd.DataFrame(mcnemar_rows)
    if not mc.empty:
        mc["mcnemar_fdr"] = np.nan
        for ptype in mc["ptype"].unique():
            m = (mc["ptype"] == ptype).to_numpy()
            p = mc.loc[m, "mcnemar_p"].to_numpy(float)
            fdr = _fdr_bh(p)
            idx = np.flatnonzero(m)
            mc.iloc[idx, mc.columns.get_loc("mcnemar_fdr")] = fdr

    peak_df = pd.DataFrame(peak_rows)
    sample_df = pd.DataFrame(sample_rows)
    summary_df = pd.DataFrame(summary_rows)

    for name, df in (("peak_agreement", peak_df),
                     ("sample_agreement", sample_df),
                     ("mcnemar", mc),
                     ("summary", summary_df)):
        if df.empty:
            df = pd.DataFrame(columns=_COLS[name])
        df.to_csv(out_table(cfg, f"P2_{name}.tsv"), sep="\t", index=False)

    if skipped:
        print(f"  [P2] skipped comparisons: {', '.join(skipped)}")

    return dict(peak_agreement=peak_df, sample_agreement=sample_df,
                mcnemar=mc, summary=summary_df, skipped=skipped)
