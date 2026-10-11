"""Module P10 (pairwise): burden-controlled permutation test on peaks.

Adopted from the controlled-permutation methodology of Yuan et al., Cancer
Cell 2018 (reproduction manual R05 §5).  It upgrades the per-peak
comparison of module P3 with the one confounder P3 cannot handle: **global
burden**.  A cohort whose tumours carry more SCNAs overall shows higher
alteration frequency at *every* peak, so a naive per-peak test is
systematically inflated.  The permutation randomises *which* sample-peak
pairs are altered while holding fixed both margins of the binary matrix:

  row sums    = each sample's total number of altered peaks
                (the "global disruption" proxy)
  column sums = each peak's overall alteration frequency

so that only peak-specific *preferential hitting* remains testable.

Matrix: union peaks (P1-style matching chains) x samples, from each
cohort's own all_lesions ``- CN values`` rows; X[i, j] = 1 when the
sample's copy value crosses the peak's direction at +/- ``cn_thr``
(default 0.25).  Direction-matched: an amplification peak only counts a
sample that is actually amplified there.

Statistics: per peak the Wald Z of ``altered ~ cohort`` (+ covariates when
metadata is provided and ``stat: wald``).  Without covariates the default
``stat: zprop`` uses the two-proportion z statistic — asymptotically the
same Wald Z and fully vectorised, so 10,000 margin-preserving permutations
stay fast.  Empirical two-tailed p = (#{|Z_perm| >= |Z_obs|} + 1)/(B + 1),
BH across peaks.

Every peak also gets its **naive** p (the same statistic without margin
control) so the burden-induced inflation is quantified: peaks that are
naive-significant but **not** permutation-robust (``burden_flag``) are
exactly the hits a global burden difference would have manufactured.

Peak QC (R05 §5.4): invariable peaks, peaks fully inside the centromere or
the terminal 1 Mb of a chromosome (probe-poor), and samples with more than
``max_segments`` segments are excluded.

Outputs: P10_peak_perm.tsv, P10_summary.tsv
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from .config import out_table
from .pairwise_model import _load_metadata, build_union_keys
from .pairwise_regional import _regional_cfg


def _permtest_cfg(cfg: dict[str, Any], rcfg: dict) -> dict:
    p = (cfg.get("analysis") or {}).get("permtest") or {}
    return dict(
        cn_thr=float(p.get("cn_thr", 0.25)),
        n_perm=int(p.get("n_perm", 10_000)),
        stat=str(p.get("stat", "zprop")),          # zprop | wald
        exclude_centromere_telomere=bool(
            p.get("exclude_centromere_telomere", True)),
        max_segments=int(p.get("max_segments", 2_000)),
        seed=int(p.get("seed", rcfg.get("seed", 42))),
        fdr=float(p.get("fdr", 0.05)),
    )


# ---------------------------------------------------------------------------
# binary matrix
# ---------------------------------------------------------------------------

def _binary_matrix(cfg: dict[str, Any], data: dict[str, dict], ca: str,
                   cb: str, keys: pd.DataFrame, cn_thr: float
                   ) -> tuple[np.ndarray, list[str], pd.DataFrame] | None:
    """Union-peak x sample binary matrix from each cohort's own CN rows.

    Returns (X, samples, peak_meta) or None when no peak is present on both
    sides; shared sample ids (an unpaired-run accident) are dropped with a
    notice."""
    cn = {cid: data[cid]["peak_cn"] for cid in (ca, cb)}
    shared = (set(cn[ca]["sample"]) & set(cn[cb]["sample"]))
    if shared:
        print(f"  [P10] dropping {len(shared)} shared sample id(s) "
              "(unpaired comparison)")
    ksub = keys[keys["cohort"].isin([ca, cb])]
    peak_meta = []
    cols: list[tuple[pd.Series, float]] = []
    for key, g in ksub.groupby("peak_key", sort=False):
        ptype = str(g["ptype"].iloc[0])
        first = g.sort_values(["chrom", "start"]).iloc[0]
        present = {}
        vals = {}
        for cid in (ca, cb):
            gc = g[g["cohort"] == cid]
            if gc.empty:
                present[cid] = False
                continue
            sub = cn[cid][cn[cid]["peak_id"].isin(set(gc["peak_id"]))]
            v = (sub.groupby("sample")["cn"].mean() if not sub.empty
                 else pd.Series(dtype=float))
            present[cid] = True
            vals[cid] = v
        if not (present.get(ca) and present.get(cb)):
            continue          # a single-sided peak has no between-cohort test
        col = pd.concat([vals[ca], vals[cb]])
        col = col[~col.index.duplicated(keep="first")]
        cols.append((col, cn_thr if ptype == "Amp" else -cn_thr))
        peak_meta.append(dict(
            peak_key=key, ptype=ptype, cytoband=first["cytoband"],
            chrom=first["chrom"], start=int(first["start"]),
            end=int(first["end"])))
    if not cols:
        return None
    samples = sorted(set().union(*[set(v.index) for v, _ in cols])
                     - shared)
    sidx = {s: i for i, s in enumerate(samples)}
    X = np.zeros((len(samples), len(cols)), dtype=np.int8)
    for j, (col, thr) in enumerate(cols):
        # values must be selected by the *filtered* sample labels: dropped
        # shared ids used to shift every later value onto the wrong sample
        keep = [s for s in col.index if s in sidx]
        idx = [sidx[s] for s in keep]
        vv = col.loc[keep].to_numpy(float)
        if thr > 0:
            X[idx, j] = vv >= thr
        else:
            X[idx, j] = vv <= thr
    return X, samples, pd.DataFrame(peak_meta)


def _qc_filter(cfg: dict[str, Any] | None, X: np.ndarray,
               meta: pd.DataFrame, samples: list[str], data: dict[str, dict],
               pcfg: dict, rcfg: dict
               ) -> tuple[np.ndarray, pd.DataFrame, list[str]]:
    """Drop invariable peaks, centromere/telomere peaks, and samples with
    too many segments."""
    from .pairwise_arm import _arm_intervals
    from .pairwise_regional import _load_band_genes
    notes = []
    keep_col = X.std(axis=0) > 0
    n_const = int((~keep_col).sum())
    if n_const:
        notes.append(f"{n_const} invariable peak(s) dropped")
    if pcfg["exclude_centromere_telomere"]:
        bands, _ = _load_band_genes(cfg, rcfg, data)
        arms = _arm_intervals(bands)
        keep_ct = np.ones(X.shape[1], dtype=bool)
        if not arms.empty:
            for j, r in meta.iterrows():
                if _is_centromere_telomere(r, arms):
                    keep_ct[j] = False
        n_ct = int((~keep_ct).sum())
        if n_ct:
            notes.append(f"{n_ct} centromere/telomere peak(s) dropped")
        keep_col &= keep_ct
    X, meta = X[:, keep_col], meta[keep_col].reset_index(drop=True)
    keep_s = np.ones(len(samples), dtype=bool)
    for cid in data:
        sc = data[cid]["seg_counts"]
        bad = set(sc.loc[sc["segment_count"] > pcfg["max_segments"],
                         "sample"])
        if bad:
            keep_s &= ~np.array([s in bad for s in samples])
    n_drop = int((~keep_s).sum())
    if n_drop:
        notes.append(f"{n_drop} sample(s) with >{pcfg['max_segments']} "
                     "segments dropped")
    X = X[keep_s]
    samples = [s for s, k in zip(samples, keep_s) if k]
    return X, meta, samples


def _is_centromere_telomere(peak: pd.Series, arms: pd.DataFrame) -> bool:
    """Peak fully inside the centromere gap or the terminal 1 Mb."""
    g = arms[arms["chrom"] == peak["chrom"]]
    if g.empty:
        return False
    tel_lo = int(g["start"].min()) + 1_000_000
    tel_hi = int(g["end"].max()) - 1_000_000
    if peak["end"] <= tel_lo or peak["start"] >= tel_hi:
        return True
    gp = g[g["arm"].str.endswith("p")]
    gq = g[g["arm"].str.endswith("q")]
    if gp.empty or gq.empty:
        return False
    return (peak["start"] >= int(gp["end"].max())
            and peak["end"] <= int(gq["start"].min()))


# ---------------------------------------------------------------------------
# margin-preserving permutation (checkerboard swaps, vectorised)
# ---------------------------------------------------------------------------

def _perturb(X: np.ndarray, rng: np.random.Generator, n_attempt: int,
             batch: int | None = None) -> int:
    """Attempt margin-preserving 2x2 swaps on binary X (in place).

    Both legal swap patterns toggle the same four cells, so a swap is
    applied by toggling them; batches whose four-cell sets collide are
    dropped cell-wise, which keeps row and column sums **exactly** fixed.
    Returns the number of swaps applied.

    ``batch`` defaults to a size *relative to the matrix*: within one huge
    batch almost every valid swap shares a cell with another valid swap and
    the collision filter drops both sides, so the chain barely moves (on a
    60x8 matrix a 2000-proposal batch applied 2 swaps out of 4800 attempts,
    vs ~10% of attempts at batch 25) — the permutation p-values from a
    frozen chain are garbage (they collapse to 1/(B+1) or 1).
    """
    n, m = X.shape
    if batch is None:
        batch = int(np.clip(n * m // 16, 32, 2_000))
    done = 0
    remaining = int(n_attempt)
    while remaining > 0:
        k = min(batch, remaining)
        remaining -= k
        i1 = rng.integers(0, n, k)
        i2 = rng.integers(0, n, k)
        j1 = rng.integers(0, m, k)
        j2 = rng.integers(0, m, k)
        ok = (i1 != i2) & (j1 != j2)
        a = X[i1, j1]
        b = X[i1, j2]
        c = X[i2, j1]
        d = X[i2, j2]
        swap = ok & (((a == 1) & (d == 1) & (b == 0) & (c == 0))
                     | ((a == 0) & (d == 0) & (b == 1) & (c == 1)))
        idx = np.flatnonzero(swap)
        if len(idx) == 0:
            continue
        cells = np.column_stack([i1[idx] * m + j1[idx],
                                 i1[idx] * m + j2[idx],
                                 i2[idx] * m + j1[idx],
                                 i2[idx] * m + j2[idx]])
        uniq, counts = np.unique(cells, return_counts=True)
        bad_cells = uniq[counts > 1]
        keep = ~np.isin(cells, bad_cells).any(axis=1)
        sel = idx[keep]
        if len(sel) == 0:
            continue
        X[i1[sel], j1[sel]] = 1 - X[i1[sel], j1[sel]]
        X[i1[sel], j2[sel]] = 1 - X[i1[sel], j2[sel]]
        X[i2[sel], j1[sel]] = 1 - X[i2[sel], j1[sel]]
        X[i2[sel], j2[sel]] = 1 - X[i2[sel], j2[sel]]
        done += len(sel)
    return done


# ---------------------------------------------------------------------------
# statistics
# ---------------------------------------------------------------------------

def _zprop(X: np.ndarray, is_b: np.ndarray) -> np.ndarray:
    """Two-proportion z per peak (cohort b - cohort a), vectorised."""
    ib = is_b.astype(bool)
    ia = ~ib
    ka = X[ia].sum(axis=0)
    kb = X[ib].sum(axis=0)
    na, nb = int(ia.sum()), int(ib.sum())
    with np.errstate(invalid="ignore", divide="ignore"):
        pooled = (ka + kb) / (na + nb)
        se = np.sqrt(pooled * (1 - pooled) * (1 / na + 1 / nb))
        return (kb / nb - ka / na) / se


def _z_wald(X: np.ndarray, is_b: np.ndarray, cov: np.ndarray | None
            ) -> np.ndarray:
    """Per-peak logistic Wald Z of the cohort coefficient (slow path)."""
    import statsmodels.api as sm
    Z = np.full(X.shape[1], np.nan)
    for j in range(X.shape[1]):
        y = X[:, j].astype(float)
        if len(np.unique(y)) < 2:
            continue
        colset = [np.ones(len(y)), is_b.astype(float)]
        if cov is not None and cov.size:
            colset.extend(c for c in cov.T)
        try:
            fit = sm.GLM(y, np.column_stack(colset),
                         family=sm.families.Binomial()).fit(disp=0)
            Z[j] = float(fit.params[1] / fit.bse[1])
        except Exception:
            continue
    return Z


def _fdr(p: np.ndarray) -> np.ndarray:
    out = np.full_like(p, np.nan, dtype=float)
    ok = np.isfinite(p)
    if ok.any():
        from statsmodels.stats.multitest import multipletests
        out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out


def _covariate_design(X: np.ndarray, samples: list[str], is_b: np.ndarray,
                      meta_df: pd.DataFrame
                      ) -> tuple[np.ndarray, list[str], np.ndarray,
                                 np.ndarray | None, str, int]:
    """Subset X to metadata-complete samples and build the covariate design.

    The matrix is subset **by sample label** — never by position: samples
    without metadata are dropped together with their rows, and the returned
    ``X`` / ``samples`` / ``is_b`` stay aligned.  (Indexing X by the merged
    frame's positions over the *full* sample list, as this used to do,
    raises KeyError for samples without metadata and is a no-op otherwise.)
    Returns ``(X, samples, is_b, cov, cov_note, n_dropped)`` with
    ``cov = None`` when fewer than 10 samples carry metadata.
    """
    cov_names = [c for c in meta_df.columns if c != "sample"]
    df = pd.DataFrame(dict(sample=samples, is_b=is_b.astype(float)))
    m = df.merge(meta_df, on="sample", how="inner")
    if len(m) < 10:
        return (X, samples, is_b, None,
                "none (metadata matched <10 samples)", 0)
    pos = {s: i for i, s in enumerate(samples)}
    X = X[[pos[s] for s in m["sample"]]]
    dropped = len(samples) - len(m)
    dm = pd.get_dummies(m[cov_names].astype("string"),
                        drop_first=True, dtype=float).fillna(0.0)
    return (X, list(m["sample"]), m["is_b"].to_numpy(float),
            dm.to_numpy(float), "+".join(cov_names), dropped)


def _perm_z_sample(Xw: np.ndarray, is_b: np.ndarray, B: int,
                  rng: np.random.Generator, use_wald: bool = False,
                  cov: np.ndarray | None = None
                  ) -> tuple[np.ndarray, int]:
    """Burn in, then draw Z for B margin-preserving permutations (in place).

    Returns ``(Z (B x n_peaks), n_swaps_burnin)``.  The Z matrix is
    **float64 on purpose**: with float32 storage, ties at the observed
    state round to a slightly smaller |Z| and the ``>=`` comparison in
    ``p_perm`` counts them as non-exceedances — on small/medium matrices
    (where the chain revisits the observed state often) this collapsed
    ``p_perm`` to the 1/(B+1) floor for perfectly null peaks.
    """
    n_att = 4 * Xw.shape[0] * Xw.shape[1]
    n_swaps = _perturb(Xw, rng, n_attempt=10 * n_att)     # burn-in
    Z = np.empty((B, Xw.shape[1]), dtype=float)
    for b in range(B):
        _perturb(Xw, rng, n_attempt=n_att)
        Z[b] = _z_wald(Xw, is_b, cov) if use_wald else _zprop(Xw, is_b)
    return Z, n_swaps


# ---------------------------------------------------------------------------
# module entry
# ---------------------------------------------------------------------------

def run_module_p10(cfg: dict[str, Any],
                   data: dict[str, dict]) -> dict[str, Any]:
    rcfg = _regional_cfg(cfg)
    pcfg = _permtest_cfg(cfg, rcfg)
    out_skip = out_table(cfg, "P10_skipped.txt")

    def skip(msg: str) -> dict[str, Any]:
        from pathlib import Path
        Path(out_skip).write_text(msg + "\n", encoding="utf-8")
        print(f"  [P10] skipped: {msg}")
        empty = pd.DataFrame()
        return dict(peaks=empty, summary=empty, skipped=msg, Z_perm={})

    if pcfg["stat"] == "wald" and pcfg["n_perm"] > 2_000:
        print(f"  [P10] stat=wald with n_perm={pcfg['n_perm']} is slow; "
              "consider stat=zprop (vectorised) or a smaller B")

    perm_rows: list[dict] = []
    sum_rows: list[dict] = []
    Z_perm_store: dict[str, np.ndarray] = {}
    q_main = cfg["analysis"]["q_thresholds"][0]
    keys = build_union_keys(cfg, data, q_main)
    if keys.empty:
        return skip("no peaks at q threshold")
    meta_df = _load_metadata(cfg)

    for cmp_ in cfg["comparisons"]:
        ca, cb = cmp_["a"], cmp_["b"]
        if cmp_.get("kind") == "paired":
            print(f"  [P10] {ca}_vs_{cb}: paired comparison skipped — P10 "
                  "implements unpaired tests only (paired inference is on "
                  "the roadmap)")
            continue
        comp = f"{ca}_vs_{cb}"
        built = _binary_matrix(cfg, data, ca, cb, keys, pcfg["cn_thr"])
        if built is None:
            print(f"  [P10] {comp}: no union peaks present in both cohorts "
                  "— skipped")
            continue
        X, samples, meta = built
        X, meta, samples = _qc_filter(cfg, X, meta, samples, data, pcfg, rcfg)
        if X.shape[1] < 2 or X.shape[0] < 6:
            print(f"  [P10] {comp}: too few variable peaks or samples after "
                  "QC — skipped")
            continue
        cohorts = [ca if s in set(data[ca]["peak_cn"]["sample"]) else cb
                   for s in samples]
        is_b = np.array([c == cb for c in cohorts])

        # optional covariates (metadata) — wald path only
        cov = None
        cov_note = "none"
        if meta_df is not None and len(meta_df.columns) > 1:
            X, samples, is_b, cov, cov_note, dropped = _covariate_design(
                X, samples, is_b, meta_df)
            if dropped:
                print(f"  [P10] {comp}: {dropped} sample(s) without "
                      "metadata dropped from the adjusted model")

        use_wald = pcfg["stat"] == "wald" or cov is not None
        if use_wald:
            Z_obs = _z_wald(X, is_b, cov)
        else:
            Z_obs = _zprop(X, is_b)
        with np.errstate(invalid="ignore"):
            p_naive = 2 * stats.norm.sf(np.abs(Z_obs))

        # margin-preserving permutations: burn-in + Markov chain
        rng = np.random.default_rng(pcfg["seed"])
        Xw = X.copy()
        Z_perm, n_swaps = _perm_z_sample(Xw, is_b, pcfg["n_perm"], rng,
                                         use_wald, cov)
        Z_perm_store[comp] = Z_perm

        B = pcfg["n_perm"]
        with np.errstate(invalid="ignore"):
            exceed = (np.abs(Z_perm) >= np.abs(Z_obs)[None, :]).sum(axis=0)
            p_perm = (exceed + 1) / (B + 1)
        q_perm = _fdr(p_perm)

        zo = Z_obs[np.isfinite(Z_obs)]
        lam_naive = float(np.median(zo ** 2) / 0.4549) if len(zo) else np.nan
        zp = Z_perm[np.isfinite(Z_perm)]
        lam_perm = float(np.median(zp.astype(float) ** 2) / 0.4549) \
            if len(zp) else np.nan

        freq_a = X[~is_b.astype(bool)].mean(axis=0)
        freq_b = X[is_b.astype(bool)].mean(axis=0)
        for j, r in meta.iterrows():
            perm_rows.append(dict(
                comparison=comp, peak_key=r["peak_key"], ptype=r["ptype"],
                cytoband=r["cytoband"], chrom=r["chrom"],
                start=int(r["start"]), end=int(r["end"]),
                freq_a=round(float(freq_a[j]), 4),
                freq_b=round(float(freq_b[j]), 4),
                delta_freq=round(float(freq_b[j] - freq_a[j]), 4),
                Z_obs=(round(float(Z_obs[j]), 3)
                       if np.isfinite(Z_obs[j]) else np.nan),
                p_naive=(round(float(p_naive[j]), 5)
                         if np.isfinite(p_naive[j]) else np.nan),
                p_perm=round(float(p_perm[j]), 5),
                q_perm=(round(float(q_perm[j]), 5)
                        if np.isfinite(q_perm[j]) else np.nan),
                burden_flag=bool(np.isfinite(p_naive[j])
                                 and p_naive[j] < 0.05
                                 and p_perm[j] >= pcfg["fdr"]),
            ))
        n_naive = int(np.nansum(p_naive < 0.05))
        n_perm_sig = int(np.nansum(q_perm < pcfg["fdr"]))
        n_burden = int(sum(1 for r in perm_rows
                           if r["comparison"] == comp and r["burden_flag"]))
        sum_rows.append(dict(
            comparison=comp, stat=("wald" if use_wald else "zprop"),
            covariates=cov_note, n_peaks=X.shape[1],
            n_samples=X.shape[0], n_perm=B, n_swaps_burnin=n_swaps,
            n_sig_naive_05=n_naive, n_sig_perm_q05=n_perm_sig,
            n_burden_inflated=n_burden,
            lambda_naive=(round(lam_naive, 3)
                          if np.isfinite(lam_naive) else np.nan),
            lambda_perm=(round(lam_perm, 3)
                         if np.isfinite(lam_perm) else np.nan),
        ))
        print(f"  [P10] {comp}: peaks={X.shape[1]} samples={X.shape[0]} "
              f"stat={'wald' if use_wald else 'zprop'} "
              f"covariates={cov_note}; sig naive(p<0.05)={n_naive} -> "
              f"perm(q<{pcfg['fdr']})={n_perm_sig} "
              f"(burden-inflated: {n_burden}; "
              f"lambda {lam_naive:.2f} -> {lam_perm:.2f})")

    peaks_df = pd.DataFrame(perm_rows)
    summary = pd.DataFrame(sum_rows)
    peaks_df.to_csv(out_table(cfg, "P10_peak_perm.tsv"), sep="\t",
                    index=False)
    summary.to_csv(out_table(cfg, "P10_summary.tsv"), sep="\t", index=False)
    if peaks_df.empty:
        return skip("no comparable peak set (all comparisons paired or "
                    "single-sided)")
    if out_skip.exists():
        out_skip.unlink()
    return dict(peaks=peaks_df, summary=summary, skipped=None,
                Z_perm=Z_perm_store, pcfg=pcfg)
