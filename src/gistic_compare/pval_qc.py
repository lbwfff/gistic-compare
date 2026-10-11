"""Module PQ (pairwise): p-value QC — Q-Q plots and genomic inflation.

Independent of every other module (like CF, it only *reads* saved tables).
It runs the quality check Yuan et al. 2018 apply to every genome-wide scan
(R05 §11): a Q-Q plot of each p-value set against the uniform null plus a
genomic inflation factor

    lambda = median(z^2) / 0.4549        (z = Phi^-1(1 - p/2))

where 0.4549 is the median of chi^2_1.  lambda ~ 1 is healthy; lambda > 1.1
signals systematic bias — batch effects, population stratification or an
unadjusted confounder — and everything downstream should be re-read before
being quoted.

P10's naive and permutation p-values are deliberately kept as **separate**
sets: the difference of their lambdas is a direct read-out of the
burden-induced inflation the permutation removed.

Set discovery: the module scans the run's table directory for the standard
outputs of P3/P4/P5/GS/P7/P10 and extracts every p-value column it knows.
Missing files are silently skipped (a module subset run still works); a
run with no readable p-value tables writes PQ_skipped.txt.

Outputs: PQ_inflation.tsv
Figures: PQ_qq_<table>__<column> (one Q-Q per p-value set)
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from .config import out_table

# (table stem, [p-value columns])
_SPECS: list[tuple[str, list[str]]] = [
    ("P4_regional_diff", ["p_amp", "p_del", "p_mw"]),
    ("P5_regional_diff", ["p_amp", "p_del", "p_mw"]),
    ("P4_cytoband_diff", ["p_amp", "p_del", "p_mw"]),
    ("P5_cytoband_diff", ["p_amp", "p_del", "p_mw"]),
    ("GS_window_gscore", ["p_amp", "p_del"]),
    ("GS_band_gscore", ["p_amp", "p_del"]),
    ("P7_window_amplitude", ["p_freq_amp", "p_freq_del", "p_amp_amp",
                             "p_amp_del", "p_high_amp", "p_high_del"]),
    ("P7_band_amplitude", ["p_freq_amp", "p_freq_del", "p_amp_amp",
                           "p_amp_del", "p_high_amp", "p_high_del"]),
    ("P7_peak_amplitude", ["p_freq_amp", "p_freq_del", "p_amp_amp",
                           "p_amp_del", "p_high_amp", "p_high_del"]),
    ("P3_peak_models", ["crude_p"]),
    ("P10_peak_perm", ["p_naive", "p_perm"]),
]


def _clean_p(p: np.ndarray) -> np.ndarray:
    """Finite p-values in [0, 1] — the set every PQ statistic uses.

    Exact zeros are **kept**: they are underflowed p-values (the most
    significant tests), not missing data, and dropping them biased lambda
    downward while the reported ``n`` still counted them.
    """
    p = np.asarray(p, dtype=float)
    return p[np.isfinite(p) & (p >= 0) & (p <= 1)]


def _lambda(p: np.ndarray) -> float:
    """Genomic inflation factor from a p-value set (genomic control)."""
    p = _clean_p(p)
    if len(p) < 10:
        return np.nan
    z = stats.norm.isf(np.clip(p, 1e-300, 1.0) / 2.0)
    return float(np.median(z * z) / 0.4549)


def _qq_arrays(p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(expected, observed) -log10 quantiles against the uniform null.

    Uses the same cleaned set as ``_lambda`` (so p > 1 corrupt values can
    never enter) minus the exact zeros, which have no finite position on a
    log scale — the Q-Q is a shape diagnostic; lambda and the reported n
    use the full set including zeros.
    """
    p = _clean_p(p)
    p = np.sort(p[p > 0])
    n = len(p)
    if n == 0:
        return np.array([]), np.array([])
    expected = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    observed = -np.log10(p)
    return expected, observed


def run_module_pq(cfg: dict[str, Any],
                  data: dict[str, dict]) -> dict[str, Any]:
    tdir = Path(cfg["output"]["_tables"])
    rows: list[dict] = []
    qq: list[dict] = []
    n_sets = 0

    for stem, pcols in _SPECS:
        path = tdir / f"{stem}.tsv"
        if not path.exists():
            continue
        try:
            df = pd.read_csv(path, sep="\t")
        except Exception:
            continue
        if df.empty:
            continue
        has_comp = "comparison" in df.columns
        for col in pcols:
            if col not in df.columns:
                continue
            n_sets += 1
            set_id = f"{stem}__{col}"
            pooled = df[col].to_numpy(float)
            exp, obs = _qq_arrays(pooled)
            if len(obs):
                qq.append(dict(set_id=set_id, expected=exp, observed=obs))
            comps = (list(dict.fromkeys(df["comparison"])) if has_comp
                     else ["all"])
            for comp in comps:
                sub = df[df["comparison"] == comp] if has_comp else df
                p = _clean_p(sub[col].to_numpy(float))
                if len(p) == 0:
                    continue
                lam = _lambda(p)          # same cleaned set as n/median
                rows.append(dict(
                    set_id=set_id, table=stem, column=col,
                    comparison=comp, n=len(p),
                    n_p_lt_005=int((p < 0.05).sum()),
                    median_p=round(float(np.median(p)), 5),
                    lambda_infl=(round(lam, 3) if np.isfinite(lam)
                                 else np.nan),
                ))

    inflation = pd.DataFrame(rows)
    if inflation.empty:
        skip_path = out_table(cfg, "PQ_skipped.txt")
        Path(skip_path).write_text(
            "no readable p-value tables found (run P3/P4/P5/GS/P7/P10 "
            "before PQ, or in the same run)\n", encoding="utf-8")
        print("  [PQ] skipped: no readable p-value tables found")
        empty = pd.DataFrame()
        return dict(inflation=empty, qq=[], skipped=True)

    Path(out_table(cfg, "PQ_skipped.txt")).unlink(missing_ok=True)
    inflation.to_csv(out_table(cfg, "PQ_inflation.tsv"), sep="\t",
                     index=False)
    lam = inflation.dropna(subset=["lambda_infl"])
    if len(lam):
        high = lam[lam["lambda_infl"] > 1.1]
        print(f"  [PQ] {n_sets} p-value set(s), {len(inflation)} "
              f"comparison rows; lambda "
              f"{lam['lambda_infl'].min():.2f}-"
              f"{lam['lambda_infl'].max():.2f}")
        if len(high):
            print("  [PQ] lambda > 1.1 (systematic bias — re-read the rest "
                  "of the run): "
                  + ", ".join(f"{r.set_id}({r.comparison})"
                              for r in high.itertuples()))
    return dict(inflation=inflation, qq=qq, skipped=False)
