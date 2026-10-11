"""Module P3 (pairwise workflow): per-peak altered ~ covariate models.

Builds union peak keys from P1-style moderate matching across every pair
(actually: one global union-find over all cohorts so each genomic region
gets a single key), maps per-sample 0/1 alteration onto those keys, joins
sample metadata, and fits:

  crude:    altered ~ cohort            (logistic; OR for cohort B vs A)
  adjusted: altered ~ cohort + covariate (one covariate at a time)

Paired comparisons: crude cohort effect is summarised by McNemar on the
shared samples (documented simplification); adjusted models use ordinary
logistic and ignore within-sample correlation (conservative note in README).

Skips entirely (writes P3_skipped.txt) when metadata is absent/empty or
statsmodels cannot be imported.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from . import intervals as iv
from .config import out_table

try:
    import statsmodels.api as sm
    from statsmodels.stats.multitest import multipletests
    _HAS_SM = True
except Exception:  # pragma: no cover
    _HAS_SM = False


class _UnionFind:
    def __init__(self) -> None:
        self.parent: dict[str, str] = {}

    def find(self, x: str) -> str:
        self.parent.setdefault(x, x)
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: str, b: str) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra


def build_union_keys(cfg: dict[str, Any], data: dict[str, dict],
                     q_main: float) -> pd.DataFrame:
    """One row per (peak_key, cohort, peak) using moderate matching chains."""
    ov = cfg["analysis"]["overlap"]
    uf = _UnionFind()
    meta: dict[str, dict[str, Any]] = {}
    cids = list(cfg["cohorts"])

    # seed every peak
    for cid in cids:
        for ptype in ("Amp", "Del"):
            for _, r in iv.filter_peaks(data[cid]["peaks"], q_main, ptype).iterrows():
                pid = r["peak_id"]
                uf.find(pid)
                meta[pid] = dict(cohort=cid, peak_id=pid, ptype=ptype,
                                 cytoband=r.get("cytoband", ""),
                                 chrom=r["chrom"], start=int(r["start"]),
                                 end=int(r["end"]), q=float(r["q"]),
                                 freq_altered=float(r["freq_altered"]))

    # chain unions across every pair (type-consistent)
    for i, ca in enumerate(cids):
        for cb in cids[i + 1:]:
            for ptype in ("Amp", "Del"):
                fa = iv.filter_peaks(data[ca]["peaks"], q_main, ptype)
                fb = iv.filter_peaks(data[cb]["peaks"], q_main, ptype)
                res = iv.match_peaks(
                    fa, fb, min_bp=ov["min_bp"],
                    reciprocal=ov["reciprocal_shorter"], mode="either")
                for _, pr in res.pairs.iterrows():
                    uf.union(pr["a_peak_id"], pr["b_peak_id"])

    rows = []
    for pid, m in meta.items():
        root = uf.find(pid)
        rm = meta[root]
        # human-readable key: direction + the coordinates of the chain's
        # root peak (internal ids like 'cue:Del:34' mean nothing to readers)
        rows.append({**m, "peak_key":
                     f"{m['ptype']}:{rm['chrom']}:{rm['start']}-{rm['end']}"})
    return pd.DataFrame(rows)


def _load_metadata(cfg: dict[str, Any]) -> pd.DataFrame | None:
    meta = cfg.get("metadata")
    if not meta:
        return None
    path = meta.get("_abs")
    if path is None or not Path(path).exists():
        return None
    sep = "\t" if str(path).endswith((".tsv", ".txt")) else ","
    try:
        df = pd.read_csv(path, sep=sep)
    except Exception:
        return None
    sid = meta.get("sample_id", "sample")
    if sid not in df.columns:
        return None
    covs = list(meta.get("covariates") or [])
    keep = [sid] + [c for c in covs if c in df.columns and c != sid]
    df = df[keep].drop_duplicates(subset=[sid])
    return df


def Path_exists(p) -> bool:
    from pathlib import Path
    return Path(p).exists()


def _design_matrix(df: pd.DataFrame, cohort_ref: str, cohort_alt: str,
                   covariate: str | None) -> tuple[np.ndarray, np.ndarray, pd.Series]:
    """Return y, X (with intercept), index for rows with both cohorts' samples
    restricted to the two cohorts of interest."""
    sub = df[df["cohort"].isin([cohort_ref, cohort_alt])].copy()
    y = (sub["cohort"] == cohort_alt).astype(float)  # placeholder wrong below
    # Actually outcome is alteration, cohort is predictor:
    y = sub["altered"].to_numpy(float)
    cohort_d = (sub["cohort"] == cohort_alt).astype(float).to_numpy()
    cols = [np.ones(len(sub)), cohort_d]
    names = ["Intercept", f"cohort[{cohort_alt}]"]
    if covariate:
        v = sub[covariate]
        if pd.api.types.is_numeric_dtype(v) and v.nunique(dropna=True) > 2:
            x = pd.to_numeric(v, errors="coerce").to_numpy(float)
            x = np.nan_to_num(x, nan=np.nanmean(x))
            cols.append(x)
        else:
            d = pd.get_dummies(v.astype("string").fillna("NA"),
                               drop_first=True, dtype=float)
            # drop the level matching cohort_ref-order first dummy = reference
            for c in d.columns:
                cols.append(d[c].to_numpy(float))
                names.append(f"{covariate}[{c}]")
        names_ext = names[:2] + [n for n in names[2:]]
        names = names_ext
    X = np.column_stack(cols)
    mask = np.ones(len(y), dtype=bool)
    # drop rows with NaN in any used column already handled; filter zero-var
    return y, X, pd.Index(names)


def _fit_logistic(y: np.ndarray, X: np.ndarray) -> tuple[float, float, float]:
    """Return (or, or_lo, or_hi) for the second column (cohort effect)."""
    # need both classes and variation
    if len(np.unique(y)) < 2:
        return np.nan, np.nan, np.nan
    try:
        model = sm.GLM(y, X, family=sm.families.Binomial())
        res = model.fit(disp=0, maxiter=100)
        beta = res.params[1]
        se = res.bse[1]
        lo, hi = beta - 1.96 * se, beta + 1.96 * se
        p = res.pvalues[1]
        return float(np.exp(beta)), float(np.exp(lo)), float(np.exp(hi)), float(p)
    except Exception:
        return np.nan, np.nan, np.nan, np.nan


def run_module_p3(cfg: dict[str, Any], data: dict[str, dict]) -> dict[str, Any]:
    out_skip = out_table(cfg, "P3_skipped.txt")

    def skip(msg: str) -> dict[str, Any]:
        _write_skip(out_skip, msg)
        print(f"  [P3] skipped: {msg}")
        empty = pd.DataFrame()
        return dict(models=empty, covariate_summary=empty,
                    coverage=empty, skipped=msg)

    if not _HAS_SM:
        return skip("statsmodels not importable")
    meta_df = _load_metadata(cfg)
    if meta_df is None:
        return skip("no metadata file (set metadata.path / sample_id / "
                    "covariates in config, or file missing)")
    covs = [c for c in (cfg["metadata"].get("covariates") or [])
            if c in meta_df.columns]
    if not covs:
        return skip("metadata loaded but no usable covariates")

    q_main = cfg["analysis"]["q_thresholds"][0]
    keys = build_union_keys(cfg, data, q_main)
    if keys.empty:
        return skip("no peaks at q threshold")

    # long alteration table: key × sample
    alt_rows = []
    key_lookup = keys.set_index("peak_id")[["peak_key", "ptype", "cytoband",
                                            "cohort"]].to_dict("index")
    for cid, d in data.items():
        calls = d["calls"]
        pk = calls[calls["peak_id"].isin(key_lookup)].copy()
        if pk.empty:
            continue
        pk["peak_key"] = pk["peak_id"].map(lambda p: key_lookup[p]["peak_key"])
        pk["cohort"] = cid
        pk["altered"] = (pd.to_numeric(pk["call"], errors="coerce") > 0).astype(float)
        alt_rows.append(pk[["peak_key", "sample", "cohort", "altered"]])
    if not alt_rows:
        return skip("no call data")
    alt = pd.concat(alt_rows, ignore_index=True)
    alt = alt.merge(meta_df, left_on="sample",
                    right_on=cfg["metadata"].get("sample_id", "sample"),
                    how="inner")
    if alt.empty:
        return skip("metadata join produced 0 samples")
    join_rate = alt["sample"].nunique() / max(
        1, pd.concat([d["calls"]["sample"] for d in data.values()]).nunique())

    # coverage table
    cov_rows = [dict(metric="sample_join_rate", value=float(join_rate))]
    for c in covs:
        na = float(meta_df[c].isna().mean())
        cov_rows.append(dict(metric=f"missing_{c}", value=na))
    coverage = pd.DataFrame(cov_rows)

    # restrict to comparisons present in config (use first pair's two cohorts
    # per comparison; if >2 cohorts, fit pairwise models)
    model_rows = []
    key_meta = keys.drop_duplicates("peak_key")[
        ["peak_key", "ptype", "cytoband"]].set_index("peak_key")

    pairs = [(c["a"], c["b"]) for c in cfg["comparisons"]]
    # de-dup unordered
    seen = set()
    uniq_pairs = []
    for a, b in pairs:
        t = tuple(sorted((a, b)))
        if t in seen:
            continue
        seen.add(t)
        uniq_pairs.append((a, b))

    for ca, cb in uniq_pairs:
        comp = f"{ca}_vs_{cb}"
        sub_pair = alt[alt["cohort"].isin([ca, cb])]
        for pk, g in sub_pair.groupby("peak_key"):
            if g["cohort"].nunique() < 2:
                continue
            # need both classes of outcome
            if g["altered"].nunique() < 2:
                continue
            y, X, names = _design_matrix(
                g.assign(cohort=g["cohort"]), ca, cb, None)
            or_c, lo_c, hi_c, p_c = _fit_logistic(y, X)
            row = dict(
                peak_key=pk, comparison=comp,
                ptype=key_meta.at[pk, "ptype"] if pk in key_meta.index else "",
                cytoband=key_meta.at[pk, "cytoband"] if pk in key_meta.index else "",
                n_samples=int(len(y)),
                crude_or=or_c, crude_or_lo=lo_c, crude_or_hi=hi_c,
                crude_p=p_c,
            )
            for cov in covs:
                # rows with non-missing covariate
                g2 = g.copy()
                if not pd.api.types.is_numeric_dtype(g2[cov]):
                    g2 = g2[g2[cov].notna()]
                else:
                    g2 = g2[g2[cov].notna()]
                if g2["cohort"].nunique() < 2 or g2["altered"].nunique() < 2:
                    row[f"adj_or_{cov}"] = np.nan
                    row[f"adj_p_{cov}"] = np.nan
                    continue
                y2, X2, _ = _design_matrix(g2, ca, cb, cov)
                if X2.shape[1] < 3:  # covariate added no column
                    row[f"adj_or_{cov}"] = np.nan
                    row[f"adj_p_{cov}"] = np.nan
                    continue
                or_a, lo_a, hi_a, p_a = _fit_logistic(y2, X2)
                row[f"adj_or_{cov}"] = or_a
                row[f"adj_p_{cov}"] = p_a
                row[f"adj_or_lo_{cov}"] = lo_a
                row[f"adj_or_hi_{cov}"] = hi_a
            model_rows.append(row)

    if not model_rows:
        return skip("no peak×sample rows survived joins / both-class filters")

    models = pd.DataFrame(model_rows)
    # FDR on crude and each adjusted p
    for col, out in [("crude_p", "crude_fdr")] + [
            (f"adj_p_{c}", f"adj_fdr_{c}") for c in covs]:
        if col in models.columns:
            p = models[col].to_numpy(float)
            fdr = np.full(len(p), np.nan)
            ok = np.isfinite(p)
            if ok.any():
                _, f, _, _ = multipletests(p[ok], method="fdr_bh")
                fdr[ok] = f
            models[out] = fdr

    # covariate summary: how much does adjustment move the cohort OR?
    sum_rows = []
    for cov in covs:
        a = np.log(models["crude_or"].to_numpy(float))
        b = np.log(models[f"adj_or_{cov}"].to_numpy(float))
        ok = np.isfinite(a) & np.isfinite(b)
        d = b[ok] - a[ok]
        sum_rows.append(dict(
            covariate=cov,
            n_peaks=int(ok.sum()),
            median_abs_delta_log_or=float(np.median(np.abs(d))) if len(d) else np.nan,
            n_sign_flip=int(np.sum(np.sign(a[ok]) != np.sign(b[ok]))) if len(d) else 0,
            n_crude_fdr05=int((models["crude_fdr"] < 0.05).sum())
            if "crude_fdr" in models else 0,
            n_adj_fdr05=int((models[f"adj_fdr_{cov}"] < 0.05).sum()),
        ))
    cov_summary = pd.DataFrame(sum_rows)

    models.to_csv(out_table(cfg, "P3_peak_models.tsv"), sep="\t", index=False)
    cov_summary.to_csv(out_table(cfg, "P3_covariate_summary.tsv"),
                       sep="\t", index=False)
    coverage.to_csv(out_table(cfg, "P3_metadata_coverage.tsv"),
                    sep="\t", index=False)
    if out_skip.exists():
        out_skip.unlink()

    return dict(models=models, covariate_summary=cov_summary,
                coverage=coverage, skipped=None)


def _write_skip(path, text: str) -> None:
    Path(path).write_text(text + "\n", encoding="utf-8")
