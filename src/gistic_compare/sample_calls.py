"""Call-concordance analysis (shared by workflow C and pairwise module P2).

One analysis, two comparison modes:

  * benchmark ``C`` — every experiment vs the shared control cohort;
  * pairwise ``P2`` — the configured cohort pairs.

Per matched peak pair (GISTIC ``0/1/2`` calls on the same samples):

  * frequency layer   — matched-pair alteration frequencies and Δfreq,
    produced for **every** pair (paired and unpaired);
  * agreement layer   — per-pair exact agreement and quadratic-weighted κ,
    then per-sample agreement across that direction's matched pairs —
    only for ``kind=paired`` comparisons (the two cohorts must contain the
    same samples);
  * McNemar layer     — exact two-sided McNemar on the discordant pairs plus
    BH-FDR within each direction, again paired-only (and only when
    statsmodels is importable);
  * summary           — per comparison × direction.

Outputs: `<tag>_peak_agreement.tsv`, `<tag>_sample_agreement.tsv`,
         `<tag>_mcnemar.tsv`, `<tag>_frequency.tsv`, `<tag>_summary.tsv`
"""

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


def _fdr_bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    out = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    if _HAS_STATSMODELS:
        _, fdr, _, _ = _multipletests(p[ok], method="fdr_bh")
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
    "summary": ["comparison", "kind", "ptype", "n_matched_peaks",
                "n_paired_samples", "mean_exact_agreement",
                "median_exact_agreement", "mean_weighted_kappa",
                "median_weighted_kappa"],
}


def run_calls(cfg: dict[str, Any], data: dict[str, dict], tag: str
             ) -> dict[str, Any]:
    """The whole concordance analysis; ``tag`` is ``"C"`` or ``"P2"``."""
    ana = cfg["analysis"]
    q_main = ana["q_thresholds"][0]
    # analysis.overlap is part of the benchmark contract; a minimal pairwise
    # config may omit it, so fall back to the same defaults as the matcher
    ov = ana.get("overlap") or {}
    min_bp = int(ov.get("min_bp", 50))
    recip = float(ov.get("reciprocal_shorter", 0.5))

    peak_rows: list[dict[str, Any]] = []
    sample_rows: list[dict[str, Any]] = []
    mcnemar_rows: list[dict[str, Any]] = []
    freq_rows: list[pd.DataFrame] = []
    summary_rows: list[dict[str, Any]] = []

    for cmp_ in cfg["comparisons"]:
        ca, cb = cmp_["a"], cmp_["b"]
        comp = f"{ca}_vs_{cb}"
        kind = cmp_.get("kind", "unpaired")
        pa, pb = data[ca]["peaks"], data[cb]["peaks"]

        # ---- frequency layer: defined for every comparison ----
        for ptype in ("Amp", "Del"):
            fa = iv.filter_peaks(pa, q_main, ptype)
            fb = iv.filter_peaks(pb, q_main, ptype)
            res = iv.match_peaks(fa, fb, min_bp=min_bp,
                                 reciprocal=recip,
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
            m["comparison"] = comp
            m["kind"] = kind
            m["ptype"] = ptype
            m["delta_freq"] = m["a_freq"] - m["b_freq"]
            freq_rows.append(m)

        if kind != "paired":
            continue

        # ---- per-call layers: need the same samples in both cohorts ----
        calls_a = _pivot_calls(data[ca]["calls"])
        calls_b = _pivot_calls(data[cb]["calls"])
        common = sorted(set(calls_a.columns) & set(calls_b.columns))
        if not common:
            continue
        calls_a = calls_a[common]
        calls_b = calls_b[common]
        accs: list[float] = []
        kappas: list[float] = []

        for ptype in ("Amp", "Del"):
            fa = iv.filter_peaks(pa, q_main, ptype)
            fb = iv.filter_peaks(pb, q_main, ptype)
            res = iv.match_peaks(fa, fb, min_bp=min_bp,
                                 reciprocal=recip,
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

                # exact McNemar on the binary altered/not calls
                a_bin = (va[ok] > 0).astype(int)
                b_bin = (vb[ok] > 0).astype(int)
                b01 = int(((a_bin == 0) & (b_bin == 1)).sum())
                b10 = int(((a_bin == 1) & (b_bin == 0)).sum())
                mcnemar_rows.append(dict(
                    comparison=comp, ptype=ptype,
                    a_peak_id=ida, b_peak_id=idb,
                    a_cytoband=pr.get("a_cytoband", ""),
                    b_cytoband=pr.get("b_cytoband", ""),
                    only_a_altered=b10, only_b_altered=b01,
                    mcnemar_p=_exact_mcnemar_p(b10, b01),
                ))

            # per-sample agreement across this direction's matched pairs
            rows_a, rows_b = [], []
            for r in peak_rows:
                if r["comparison"] != comp or r["ptype"] != ptype:
                    continue
                if (r["a_peak_id"] in calls_a.index
                        and r["b_peak_id"] in calls_b.index):
                    rows_a.append(calls_a.loc[r["a_peak_id"]].to_numpy(float))
                    rows_b.append(calls_b.loc[r["b_peak_id"]].to_numpy(float))
            if rows_a:
                A = np.vstack(rows_a)
                B = np.vstack(rows_b)
                for si, s in enumerate(common):
                    va, vb = A[:, si], B[:, si]
                    ok = np.isfinite(va) & np.isfinite(vb)
                    if ok.sum() == 0:
                        continue
                    sample_rows.append(dict(
                        comparison=comp, ptype=ptype, sample=s,
                        n_peaks=int(ok.sum()),
                        exact_agreement=float((va[ok] == vb[ok]).mean()),
                        weighted_kappa=_weighted_kappa(va, vb),
                    ))

        for ptype in ("Amp", "Del"):
            sub = [r for r in peak_rows if r["comparison"] == comp
                   and r["ptype"] == ptype]
            if not sub:
                continue
            summary_rows.append(dict(
                comparison=comp, kind=kind, ptype=ptype,
                n_matched_peaks=len(sub),
                n_paired_samples=len(common),
                mean_exact_agreement=float(np.nanmean([r["exact_agreement"]
                                                       for r in sub])),
                median_exact_agreement=float(np.nanmedian(
                    [r["exact_agreement"] for r in sub])),
                mean_weighted_kappa=float(np.nanmean([r["weighted_kappa"]
                                                      for r in sub])),
                median_weighted_kappa=float(np.nanmedian(
                    [r["weighted_kappa"] for r in sub])),
            ))

    mc = pd.DataFrame(mcnemar_rows)
    if not mc.empty:
        mc["mcnemar_fdr"] = np.nan
        for ptype in mc["ptype"].unique():
            m = (mc["ptype"] == ptype).to_numpy()
            mc.loc[m, "mcnemar_fdr"] = _fdr_bh(
                mc.loc[m, "mcnemar_p"].to_numpy(float))

    peak_df = pd.DataFrame(peak_rows)
    sample_df = pd.DataFrame(sample_rows)
    freq_df = (pd.concat(freq_rows, ignore_index=True) if freq_rows
               else pd.DataFrame())

    for name, df in (("peak_agreement", peak_df),
                     ("sample_agreement", sample_df),
                     ("mcnemar", mc),
                     ("frequency", freq_df),
                     ("summary", pd.DataFrame(summary_rows))):
        if df.empty:
            df = pd.DataFrame(columns=_SUMMARY_COLS.get(name, df.columns))
        df.to_csv(out_table(cfg, f"{tag}_{name}.tsv"), sep="\t", index=False)

    return dict(peak_agreement=peak_df, sample_agreement=sample_df,
                mcnemar=mc, frequency=freq_df,
                summary=pd.DataFrame(summary_rows))


def run_module_c(cfg: dict[str, Any], data: dict[str, dict]) -> dict[str, Any]:
    """Benchmark form of the shared analysis (table prefix ``C_``)."""
    return run_calls(cfg, data, "C")


def run_module_p2(cfg: dict[str, Any], data: dict[str, dict]) -> dict[str, Any]:
    """Pairwise form of the shared analysis (table prefix ``P2_``)."""
    return run_calls(cfg, data, "P2")
