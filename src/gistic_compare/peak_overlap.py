"""Module A: peak-level overlap between cohorts."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from . import intervals as iv
from .config import out_table

_STAT_COLS = ["q", "residual_q", "wide_len", "freq_altered", "type", "z"]


# ---------------------------------------------------------------------------
# correlation / agreement helpers
# ---------------------------------------------------------------------------

def rank_normalize_q(q: np.ndarray) -> np.ndarray:
    """Within-vector rank-normalize -log10(q) to normal scores (Blom).

    Ties get average ranks. Result ~ N(0,1) on a common scale, so
    cross-cohort Pearson / y=x become comparable despite q inflation.
    """
    q = np.asarray(q, dtype=float)
    eps = 1e-300
    nlq = -np.log10(np.clip(q, eps, None))
    r = stats.rankdata(nlq)
    n = len(nlq)
    if n < 2:
        return np.zeros(n)
    return stats.norm.ppf((r - 3.0 / 8.0) / (n + 1.0 / 4.0))


def pearson_r(x: np.ndarray, y: np.ndarray) -> float:
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 3 or np.std(x[m]) == 0 or np.std(y[m]) == 0:
        return np.nan
    return float(np.corrcoef(x[m], y[m])[0, 1])


def bootstrap_pearson(x: np.ndarray, y: np.ndarray,
                      n_boot: int = 2000, seed: int = 42,
                      ci: float = 0.95) -> dict[str, float]:
    """Pearson r with percentile bootstrap CI (pair resampling)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    n = int(len(x))
    r = pearson_r(x, y)
    if n < 3 or not np.isfinite(r):
        return dict(r=r, ci_lo=np.nan, ci_hi=np.nan, n=n)
    rng = np.random.default_rng(seed)
    rs: list[float] = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        rr = pearson_r(x[idx], y[idx])
        if np.isfinite(rr):
            rs.append(rr)
    if not rs:
        return dict(r=r, ci_lo=np.nan, ci_hi=np.nan, n=n)
    alpha = (1.0 - ci) / 2.0
    lo, hi = np.quantile(rs, [alpha, 1.0 - alpha])
    return dict(r=r, ci_lo=float(lo), ci_hi=float(hi), n=n)


def bootstrap_spearman(x: np.ndarray, y: np.ndarray,
                       n_boot: int = 2000, seed: int = 42,
                       ci: float = 0.95) -> dict[str, float]:
    """Spearman rho with percentile bootstrap CI."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    n = int(len(x))
    if n < 3:
        return dict(r=np.nan, p=np.nan, ci_lo=np.nan, ci_hi=np.nan, n=n)
    rho, p = stats.spearmanr(x, y)
    rho = float(rho)
    rng = np.random.default_rng(seed)
    rs: list[float] = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        if np.std(x[idx]) == 0 or np.std(y[idx]) == 0:
            continue
        rr = float(stats.spearmanr(x[idx], y[idx]).statistic)
        if np.isfinite(rr):
            rs.append(rr)
    if not rs:
        return dict(r=rho, p=float(p), ci_lo=np.nan, ci_hi=np.nan, n=n)
    alpha = (1.0 - ci) / 2.0
    lo, hi = np.quantile(rs, [alpha, 1.0 - alpha])
    return dict(r=rho, p=float(p), ci_lo=float(lo), ci_hi=float(hi), n=n)


def lin_ccc(x: np.ndarray, y: np.ndarray) -> float:
    """Lin's concordance correlation coefficient (absolute agreement)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 3:
        return np.nan
    vx = float(np.var(x, ddof=1))
    vy = float(np.var(y, ddof=1))
    mx = float(np.mean(x))
    my = float(np.mean(y))
    cov = float(np.cov(x, y, ddof=1)[0, 1])
    denom = vx + vy + (mx - my) ** 2
    return float(2.0 * cov / denom) if denom > 0 else np.nan


def _attach_peak_stats(pairs: pd.DataFrame, pa: pd.DataFrame,
                       pb: pd.DataFrame, side: str = "a") -> pd.DataFrame:
    """Attach a_/b_-prefixed peak stats onto match pairs.

    cytoband already exists on pairs as a_cytoband/b_cytoband — skip it here.
    """
    out = pairs.copy()
    for pref, peaks, idcol in (("a", pa, "a_peak_id"), ("b", pb, "b_peak_id")):
        cols = [c for c in _STAT_COLS if c in peaks.columns]
        sub = peaks[["peak_id"] + cols].rename(
            columns={c: f"{pref}_{c}" for c in cols})
        out = out.merge(sub, left_on=idcol, right_on="peak_id", how="left")
        out = out.drop(columns="peak_id")
    return out


def _spearman(x, y) -> tuple[float, float, int]:
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < 3:
        return np.nan, np.nan, n
    r, p = stats.spearmanr(x[m], y[m])
    return float(r), float(p), n


def run_module_a(cfg: dict[str, Any], data: dict[str, dict]) -> dict[str, Any]:
    ana = cfg["analysis"]
    q_thresholds = ana["q_thresholds"]
    ov = ana["overlap"]
    perm = ana["permutation"]
    comparisons = cfg["comparisons"]

    jaccard_rows = []
    corr_rows = []
    match_rows = []
    unique_rows = []
    cytoband_rows = []
    perm_rows = []
    perm_nulls: list[dict[str, Any]] = []
    boot = ana.get("bootstrap", {})
    n_boot = int(boot.get("n", 2000))
    boot_seed = int(boot.get("seed", 42))
    boot_ci = float(boot.get("ci", 0.95))

    for cmp_ in comparisons:
        ca, cb = cmp_["a"], cmp_["b"]
        pa, pb = data[ca]["peaks"], data[cb]["peaks"]
        kind = cmp_.get("kind", "unpaired")

        for q in q_thresholds:
            for ptype in ("Amp", "Del"):
                fa = iv.filter_peaks(pa, q, ptype).copy()
                fb = iv.filter_peaks(pb, q, ptype).copy()
                if len(fa) == 0 or len(fb) == 0:
                    continue
                # within-set rank-normalized -log10(q) → common N(0,1) scale
                fa["z"] = rank_normalize_q(fa["q"].to_numpy())
                fb["z"] = rank_normalize_q(fb["q"].to_numpy())

                # --- interval matching: three tiers ---
                # moderate = one-sided (shorter peak covered >=50%) — primary
                # strict   = both-sided reciprocal >=50% — sensitivity
                # loose    = any overlap >= min_bp — lower bound
                tiers = (
                    ("moderate", ov["reciprocal_shorter"], "either"),
                    ("strict", ov["reciprocal_both"], "both"),
                    ("loose", ov["min_bp"], "any"),
                )
                for label, rec, mode in tiers:
                    res = iv.match_peaks(fa, fb, min_bp=ov["min_bp"],
                                         reciprocal=rec, mode=mode)
                    j = iv.jaccard(res.n_a, res.n_b, len(res.pairs))
                    jaccard_rows.append(dict(
                        comparison=f"{ca}_vs_{cb}", kind=kind, q_thresh=q,
                        ptype=ptype, method=f"interval_{label}",
                        reciprocal_mode=mode,
                        n_a=res.n_a, n_b=res.n_b, n_pairs=len(res.pairs),
                        jaccard=j,
                        overlap_frac_a=len(res.pairs) / res.n_a if res.n_a else np.nan,
                        overlap_frac_b=len(res.pairs) / res.n_b if res.n_b else np.nan,
                    ))
                    if label == "moderate" and len(res.pairs) > 0:
                        merged = _attach_peak_stats(res.pairs, fa, fb, "a")
                        merged["comparison"] = f"{ca}_vs_{cb}"
                        merged["q_thresh"] = q
                        match_rows.append(merged)

                        # Spearman on raw metrics
                        for metric in ("q", "residual_q", "wide_len", "freq_altered"):
                            r, p, n = _spearman(merged[f"a_{metric}"].to_numpy(float),
                                                merged[f"b_{metric}"].to_numpy(float))
                            corr_rows.append(dict(
                                comparison=f"{ca}_vs_{cb}", kind=kind, q_thresh=q,
                                ptype=ptype, metric=metric,
                                corr_type="spearman",
                                corr=r, pvalue=p, n_pairs=n,
                                ci_lo=np.nan, ci_hi=np.nan,
                                ccc=np.nan, median_abs_delta=np.nan,
                            ))

                        # Pearson on rank-normalized q (z) + bootstrap CI
                        # + Lin CCC and median |Δz| on the same scale
                        az = merged["a_z"].to_numpy(float)
                        bz = merged["b_z"].to_numpy(float)
                        bp = bootstrap_pearson(az, bz, n_boot=n_boot,
                                               seed=boot_seed, ci=boot_ci)
                        corr_rows.append(dict(
                            comparison=f"{ca}_vs_{cb}", kind=kind, q_thresh=q,
                            ptype=ptype, metric="ranknorm_neg_log10_q",
                            corr_type="pearson_bootstrap",
                            corr=bp["r"], pvalue=np.nan,
                            n_pairs=bp["n"],
                            ci_lo=bp["ci_lo"], ci_hi=bp["ci_hi"],
                            ccc=lin_ccc(az, bz),
                            median_abs_delta=float(np.nanmedian(np.abs(az - bz)))
                            if len(az) else np.nan,
                        ))

                        # frequency–frequency (sample-size robust) + bootstrap Spearman
                        af = merged["a_freq_altered"].to_numpy(float)
                        bf = merged["b_freq_altered"].to_numpy(float)
                        fs = bootstrap_spearman(af, bf, n_boot=n_boot,
                                                seed=boot_seed, ci=boot_ci)
                        corr_rows.append(dict(
                            comparison=f"{ca}_vs_{cb}", kind=kind, q_thresh=q,
                            ptype=ptype, metric="freq_altered",
                            corr_type="spearman_bootstrap",
                            corr=fs["r"], pvalue=fs["p"], n_pairs=fs["n"],
                            ci_lo=fs["ci_lo"], ci_hi=fs["ci_hi"],
                            ccc=lin_ccc(af, bf),
                            median_abs_delta=float(np.nanmedian(np.abs(af - bf)))
                            if len(af) else np.nan,
                        ))

                    # unique peaks (moderate — primary tier)
                    if label == "moderate":
                        for side, only in (("a", res.only_a), ("b", res.only_b)):
                            if len(only):
                                tmp = only[["peak_id", "type", "cytoband", "chrom",
                                            "start", "end", "q", "freq_altered"]].copy()
                                tmp["comparison"] = f"{ca}_vs_{cb}"
                                tmp["q_thresh"] = q
                                tmp["side"] = side
                                tmp["present_in"] = ca if side == "a" else cb
                                tmp["absent_in"] = cb if side == "a" else ca
                                unique_rows.append(tmp)

                # --- cytoband matching (sensitivity) ---
                res_c = iv.match_by_cytoband(fa, fb)
                j = iv.jaccard(res_c.n_a, res_c.n_b, len(res_c.pairs))
                cytoband_rows.append(dict(
                    comparison=f"{ca}_vs_{cb}", kind=kind, q_thresh=q,
                    ptype=ptype, method="cytoband",
                    n_a=res_c.n_a, n_b=res_c.n_b, n_pairs=len(res_c.pairs),
                    jaccard=j,
                    overlap_frac_a=len(res_c.pairs) / res_c.n_a if res_c.n_a else np.nan,
                    overlap_frac_b=len(res_c.pairs) / res_c.n_b if res_c.n_b else np.nan,
                ))

                # --- permutation null (moderate — primary tier) ---
                pr = iv.permutation_jaccard(
                    fa, fb, n_perm=perm["n"], min_bp=ov["min_bp"],
                    reciprocal=ov["reciprocal_shorter"], mode="either",
                    seed=perm["seed"],
                )
                null_j = pr.pop("null_jaccards")
                perm_nulls.append(dict(
                    comparison=f"{ca}_vs_{cb}", q_thresh=q, ptype=ptype,
                    observed_jaccard=pr["observed_jaccard"],
                    null_jaccards=null_j,
                ))
                pr.update(comparison=f"{ca}_vs_{cb}", kind=kind,
                          q_thresh=q, ptype=ptype)
                perm_rows.append(pr)

    # --- three-way peak presence (moderate matching, q per threshold) ---
    three_way = _three_way(data, q_thresholds[0])

    results = dict(
        jaccard=pd.DataFrame(jaccard_rows),
        cytoband_jaccard=pd.DataFrame(cytoband_rows),
        matched_pairs=pd.concat(match_rows, ignore_index=True) if match_rows else pd.DataFrame(),
        correlations=pd.DataFrame(corr_rows),
        unique_peaks=pd.concat(unique_rows, ignore_index=True) if unique_rows else pd.DataFrame(),
        permutation=pd.DataFrame(perm_rows),
        permutation_nulls=perm_nulls,
        three_way=three_way,
    )
    _write(cfg, results)
    return results


def _three_way(data: dict[str, dict], q: float) -> pd.DataFrame:
    """Long table: each peak × presence in each cohort (q<thresh)."""
    rows = []
    for cid, d in data.items():
        for _, r in d["peaks"].iterrows():
            if r["q"] >= q:
                continue
            rows.append(dict(cohort=cid, peak_id=r["peak_id"], type=r["type"],
                             cytoband=r["cytoband"], chrom=r["chrom"],
                             start=r["start"], end=r["end"], q=r["q"],
                             freq_altered=r["freq_altered"]))
    long = pd.DataFrame(rows)
    if long.empty:
        return long
    # link peaks across cohorts by interval overlap per type
    links = []
    for ptype in ("Amp", "Del"):
        ids = [c for c in data]
        for i in range(len(ids)):
            for j in range(i + 1, len(ids)):
                a = long[(long.cohort == ids[i]) & (long.type == ptype)]
                b = long[(long.cohort == ids[j]) & (long.type == ptype)]
                if a.empty or b.empty:
                    continue
                res = iv.match_peaks(a, b, min_bp=50, reciprocal=0.5,
                                     mode="either")
                for _, pr in res.pairs.iterrows():
                    links.append(dict(a=pr["a_peak_id"], b=pr["b_peak_id"], type=ptype))
    link_df = pd.DataFrame(links, columns=["a", "b", "type"])
    # union-find style grouping of peak_ids
    parent: dict[str, str] = {}

    def find(x: str) -> str:
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x: str, y: str) -> None:
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    for pid in long["peak_id"]:
        find(pid)
    for _, r in link_df.iterrows():
        union(r["a"], r["b"])
    long["group"] = long["peak_id"].map(find)

    # group summary
    g = (long.groupby("group")
         .agg(n_cohorts=("cohort", "nunique"),
              cohorts=("cohort", lambda s: ",".join(sorted(set(s)))),
              type=("type", "first"),
              cytobands=("cytoband", lambda s: ",".join(sorted(set(s)))),
              min_q=("q", "min"),
              peak_ids=("peak_id", lambda s: ",".join(s)))
         .reset_index())
    return g


def _write(cfg: dict[str, Any], res: dict[str, Any]) -> None:
    for name, df in res.items():
        if name == "permutation_nulls":
            continue
        if isinstance(df, pd.DataFrame):
            df.to_csv(out_table(cfg, f"A_{name}.tsv"), sep="\t", index=False)
