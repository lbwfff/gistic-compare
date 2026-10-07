"""Module C: per-sample 0/1/2 call concordance."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from . import intervals as iv
from .config import out_table

try:
    from statsmodels.stats.multitest import multipletests as _multipletests
    _HAS_STATSMODELS = True
except Exception:  # pragma: no cover
    _HAS_STATSMODELS = False


def _weighted_kappa(y_true: np.ndarray, y_pred: np.ndarray,
                    weights: str = "quadratic") -> float:
    """Quadratic-weighted kappa for ordinal 0/1/2 calls."""
    # build confusion matrix
    cats = [0, 1, 2]
    n = len(cats)
    cm = np.zeros((n, n), dtype=int)
    idx = {c: i for i, c in enumerate(cats)}
    ok = np.isfinite(y_true) & np.isfinite(y_pred)
    for t, p in zip(y_true[ok].astype(int), y_pred[ok].astype(int)):
        if t in idx and p in idx:
            cm[idx[t], idx[p]] += 1
    if cm.sum() == 0:
        return np.nan
    # weighted kappa
    O = cm / cm.sum()
    r = cm.sum(axis=1) / cm.sum()
    c = cm.sum(axis=0) / cm.sum()
    E = np.outer(r, c)
    idx_grid = np.arange(n)
    W = ((idx_grid[:, None] - idx_grid[None, :]) ** 2) / ((n - 1) ** 2)
    numer = np.sum(W * O)
    denom = np.sum(W * E)
    if denom == 0:
        return np.nan
    return float(1 - numer / denom)


def _exact_mcnemar_p(b01: int, b10: int) -> float:
    """Exact two-sided McNemar p-value from discordant counts (binomial)."""
    n = b01 + b10
    if n == 0:
        return 1.0
    k = min(b01, b10)
    # two-sided exact binomial test with p=0.5
    from scipy.stats import binomtest
    return float(binomtest(k, n, 0.5, alternative="two-sided").pvalue)


def _pivot_calls(calls: pd.DataFrame) -> pd.DataFrame:
    """peak_id × sample matrix of calls."""
    return calls.pivot_table(index="peak_id", columns="sample",
                             values="call", aggfunc="first")


# headers so empty tables still parse (star mode may have no paired rows)
_SUMMARY_COLS = {
    "peak_agreement": ["comparison", "ptype", "a_peak_id", "b_peak_id",
                       "a_cytoband", "b_cytoband", "n_samples",
                       "exact_agreement", "weighted_kappa", "a_freq", "b_freq"],
    "sample_agreement": ["comparison", "ptype", "sample", "n_peaks",
                         "exact_agreement", "weighted_kappa"],
    "mcnemar": ["comparison", "ptype", "a_peak_id", "b_peak_id",
                "a_cytoband", "b_cytoband", "only_a_altered",
                "only_b_altered", "mcnemar_p", "mcnemar_fdr"],
    "frequency": ["a_peak_id", "b_peak_id", "overlap_bp", "reciprocal_a",
                  "reciprocal_b", "a_cytoband", "b_cytoband", "comparison",
                  "kind", "ptype", "delta_freq", "a_freq", "b_freq",
                  "a_q", "b_q"],
    "summary": ["comparison", "ptype", "n_matched_peaks",
                "mean_exact_agreement", "median_exact_agreement",
                "mean_weighted_kappa", "median_weighted_kappa"],
}


def run_module_c(cfg: dict[str, Any], data: dict[str, dict]) -> dict[str, Any]:
    ana = cfg["analysis"]
    q_main = ana["q_thresholds"][0]
    ov = ana["overlap"]

    # ---------- paired: cue vs despiked ----------
    paired = [c for c in cfg["comparisons"] if c.get("kind") == "paired"]
    unpaired = [c for c in cfg["comparisons"] if c.get("kind") != "paired"]

    peak_level_rows = []
    sample_level_rows = []
    mcnemar_rows = []
    freq_rows = []

    for cmp_ in paired + unpaired:
        ca, cb = cmp_["a"], cmp_["b"]
        pa, pb = data[ca]["peaks"], data[cb]["peaks"]
        kind = cmp_.get("kind", "unpaired")

        # frequency comparison always possible
        for ptype in ("Amp", "Del"):
            fa = iv.filter_peaks(pa, q_main, ptype)
            fb = iv.filter_peaks(pb, q_main, ptype)
            res = iv.match_peaks(fa, fb, min_bp=ov["min_bp"],
                                 reciprocal=ov["reciprocal_shorter"],
                                 mode="either")
            if len(res.pairs) == 0:
                continue
            m = res.pairs.merge(
                fa[["peak_id", "freq_altered", "q"]]
                .rename(columns={"freq_altered": "a_freq", "q": "a_q"}),
                left_on="a_peak_id", right_on="peak_id").drop(columns="peak_id")
            m = m.merge(
                fb[["peak_id", "freq_altered", "q"]]
                .rename(columns={"freq_altered": "b_freq", "q": "b_q"}),
                left_on="b_peak_id", right_on="peak_id", how="left"
            ).drop(columns="peak_id")
            m["comparison"] = f"{ca}_vs_{cb}"
            m["kind"] = kind
            m["ptype"] = ptype
            m["delta_freq"] = m["a_freq"] - m["b_freq"]
            freq_rows.append(m)

        if kind != "paired":
            continue

        calls_a = _pivot_calls(data[ca]["calls"])
        calls_b = _pivot_calls(data[cb]["calls"])
        common_samples = sorted(set(calls_a.columns) & set(calls_b.columns))
        if not common_samples:
            continue
        calls_a = calls_a[common_samples]
        calls_b = calls_b[common_samples]

        for ptype in ("Amp", "Del"):
            fa = iv.filter_peaks(pa, q_main, ptype)
            fb = iv.filter_peaks(pb, q_main, ptype)
            res = iv.match_peaks(fa, fb, min_bp=ov["min_bp"],
                                 reciprocal=ov["reciprocal_shorter"],
                                 mode="either")
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
                peak_level_rows.append(dict(
                    comparison=f"{ca}_vs_{cb}", ptype=ptype,
                    a_peak_id=ida, b_peak_id=idb,
                    a_cytoband=pr["a_cytoband"], b_cytoband=pr["b_cytoband"],
                    n_samples=int(ok.sum()),
                    exact_agreement=acc,
                    weighted_kappa=kappa,
                    a_freq=float(np.nanmean(va > 0)),
                    b_freq=float(np.nanmean(vb > 0)),
                ))

                # McNemar on binary altered/not
                a_bin = (va[ok] > 0).astype(int)
                b_bin = (vb[ok] > 0).astype(int)
                b01 = int(((a_bin == 0) & (b_bin == 1)).sum())  # b altered, a not
                b10 = int(((a_bin == 1) & (b_bin == 0)).sum())  # a altered, b not
                table = np.array([[int(((a_bin == 1) & (b_bin == 1)).sum()), b10],
                                  [b01, int(((a_bin == 0) & (b_bin == 0)).sum())]])
                # 2x2 table of discordant counts for exact McNemar
                discord = np.array([[0, b10], [b01, 0]])
                pval = _exact_mcnemar_p(b10, b01)
                mcnemar_rows.append(dict(
                    comparison=f"{ca}_vs_{cb}", ptype=ptype,
                    a_peak_id=ida, b_peak_id=idb,
                    a_cytoband=pr["a_cytoband"], b_cytoband=pr["b_cytoband"],
                    only_a_altered=b10, only_b_altered=b01,
                    mcnemar_p=pval,
                ))

            # per-sample agreement across all matched peak pairs of this type
            if peak_level_rows:
                sub = [r for r in peak_level_rows
                       if r["comparison"] == f"{ca}_vs_{cb}" and r["ptype"] == ptype]
                # aggregate per sample
                peak_ids_a = [r["a_peak_id"] for r in sub]
                peak_ids_b = [r["b_peak_id"] for r in sub]
                if peak_ids_a:
                    mat_a = calls_a.loc[[p for p in peak_ids_a if p in calls_a.index]]
                    mat_b = calls_b.loc[[p for p in peak_ids_b if p in calls_b.index]]
                    # align by pair order
                    rows_a = []
                    rows_b = []
                    for r in sub:
                        if r["a_peak_id"] in calls_a.index and r["b_peak_id"] in calls_b.index:
                            rows_a.append(calls_a.loc[r["a_peak_id"]].to_numpy(float))
                            rows_b.append(calls_b.loc[r["b_peak_id"]].to_numpy(float))
                    if rows_a:
                        A = np.vstack(rows_a)
                        B = np.vstack(rows_b)
                        for si, s in enumerate(common_samples):
                            va, vb = A[:, si], B[:, si]
                            ok = np.isfinite(va) & np.isfinite(vb)
                            if ok.sum() == 0:
                                continue
                            sample_level_rows.append(dict(
                                comparison=f"{ca}_vs_{cb}", ptype=ptype,
                                sample=s,
                                n_peaks=int(ok.sum()),
                                exact_agreement=float((va[ok] == vb[ok]).mean()),
                                weighted_kappa=_weighted_kappa(va, vb),
                            ))

    # FDR on McNemar
    mc = pd.DataFrame(mcnemar_rows)
    if not mc.empty:
        mc["mcnemar_fdr"] = np.nan
        for ptype in mc["ptype"].unique():
            m = (mc["ptype"] == ptype).to_numpy()
            p = mc.loc[m, "mcnemar_p"].to_numpy(float)
            ok = np.isfinite(p)
            if ok.any():
                if _HAS_STATSMODELS:
                    _, fdr, _, _ = _multipletests(p[ok], method="fdr_bh")
                else:
                    order = np.argsort(p[ok])
                    ranked = p[ok][order]
                    fdr_sorted = ranked * len(ranked) / (np.arange(len(ranked)) + 1)
                    fdr_sorted = np.minimum.accumulate(fdr_sorted[::-1])[::-1]
                    fdr_sorted = np.clip(fdr_sorted, 0, 1)
                    fdr = np.empty_like(fdr_sorted)
                    fdr[order] = fdr_sorted
                idx = np.flatnonzero(m)[ok]
                mc.iloc[idx, mc.columns.get_loc("mcnemar_fdr")] = fdr

    peak_df = pd.DataFrame(peak_level_rows)
    sample_df = pd.DataFrame(sample_level_rows)
    freq_df = pd.concat(freq_rows, ignore_index=True) if freq_rows else pd.DataFrame()

    # summary
    summary_rows = []
    for (comp, ptype), g in peak_df.groupby(["comparison", "ptype"]) if not peak_df.empty else []:
        summary_rows.append(dict(
            comparison=comp, ptype=ptype, n_matched_peaks=len(g),
            mean_exact_agreement=g["exact_agreement"].mean(),
            median_exact_agreement=g["exact_agreement"].median(),
            mean_weighted_kappa=g["weighted_kappa"].mean(),
            median_weighted_kappa=g["weighted_kappa"].median(),
        ))
    summary_df = pd.DataFrame(summary_rows)

    for name, df in (("peak_agreement", peak_df), ("sample_agreement", sample_df),
                     ("mcnemar", mc), ("frequency", freq_df),
                     ("summary", summary_df)):
        # ensure header exists even when empty (e.g. star topology, no paired)
        if df.empty:
            df = pd.DataFrame(columns=_SUMMARY_COLS.get(name, df.columns.tolist()))
        df.to_csv(out_table(cfg, f"C_{name}.tsv"), sep="\t", index=False)

    return dict(peak_agreement=peak_df, sample_agreement=sample_df,
                mcnemar=mc, frequency=freq_df, summary=summary_df)
