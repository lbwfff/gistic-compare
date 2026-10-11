"""Module P8 (pairwise): multivariate separation between cohorts.

P3-P7 test locus by locus.  P8 asks the complementary, global question:
**"as a whole, are the two cohorts' copy-number profiles separable, how
strongly, and which loci drive the separation?"**

Feature sets (``analysis.multivariate.features``):
  grid      — 1 Mb window copy values (GISTIC marker matrix when every
              cohort has one, else the seg-derived window values P4 uses).
              Peak-call independent — the default.
  cytoband  — windows aggregated to cytobands.
  peaks     — union-peak 0/1/2 calls.  Convenient but **circular**: the
              peaks were called from the same cohorts being compared, so
              separability is expected even under the null.  Use for
              description, not inference.

Analyses per comparison:
  * PERMANOVA (pseudo-F on euclidean distances, Anderson 2001) with label
    permutation — unpaired designs permute labels freely, paired designs
    permute **within pairs** (restricted permutation).  Reports pseudo-F,
    R^2 (variance explained by cohort) and a permutation p.
  * PCA coordinates (PC1/PC2, variance explained) for the cohort-coloured
    ordination figure; paired samples are connected.
  * cross-validated out-of-fold AUC of an L2 logistic classifier
    (``penalty: elasticnet`` switches to elastic net) as a headline
    separability statistic, with a permutation p from the same statistic.
  * stability selection: how often each feature lands in the top-k by
    |coefficient| across bootstrap refits.

Separability is NOT evidence of biology: batch, platform and processing
differences separate cohorts just as well.  Read P8 together with module
CF (confounder bounds) and the stats.md caveats.

Outputs: P8_separation.tsv, P8_pca_coords.tsv, P8_stability.tsv,
         P8_summary.tsv
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

from .config import out_table
from .pairwise_regional import _regional_cfg

try:
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import GroupKFold, StratifiedKFold, cross_val_predict
    _HAS_SK = True
except Exception:  # pragma: no cover
    _HAS_SK = False


def _mv_cfg(cfg: dict[str, Any], rcfg: dict) -> dict:
    m = (cfg.get("analysis") or {}).get("multivariate") or {}
    return dict(
        features=str(m.get("features", "grid")),
        n_perm=int(m.get("n_perm", 999)),
        n_auc_perm=int(m.get("n_auc_perm", 100)),
        n_boot=int(m.get("n_boot", 50)),
        top_k=int(m.get("top_k", 20)),
        min_n=int(m.get("min_n", 5)),
        seed=int(m.get("seed", rcfg.get("seed", 42))),
        penalty=str(m.get("penalty", "l2")),
        l1_ratio=float(m.get("l1_ratio", 0.5)),
        max_nan_frac=float(m.get("max_nan_frac", 0.2)),
    )


# ---------------------------------------------------------------------------
# feature matrices
# ---------------------------------------------------------------------------

def _feature_grid(cfg: dict[str, Any], data: dict[str, dict], cids: list[str],
                  rcfg: dict, src: str) -> tuple[np.ndarray, dict[str, list[str]],
                                                 list[str]]:
    """Windows x samples for the given cohorts (one shared source)."""
    from .pairwise_amplitude import _matrix_for_cohort
    from .score_track import make_grid
    grid = make_grid(rcfg["bin_size"])
    mats = {cid: _matrix_for_cohort(cfg, data, cid, grid, rcfg, src)
            for cid in cids}
    samples = {cid: _sample_ids(cfg, data, cid, src, rcfg) for cid in cids}
    for cid in cids:
        if len(samples[cid]) != mats[cid].shape[1]:
            raise RuntimeError(
                f"[P8] {cid}: sample count mismatch between matrix "
                f"({mats[cid].shape[1]}) and ids ({len(samples[cid])})")
    X = np.hstack([mats[cid] for cid in cids]).T          # (N, F)
    feature_names = [f"{r.chrom}:{int(r.start)}-{int(r.end)}"
                     for r in grid.itertuples()]
    return X, samples, feature_names


def _sample_ids(cfg: dict[str, Any], data: dict[str, dict], cid: str,
                src: str, rcfg: dict) -> list[str]:
    """Sample order of the cohort's window matrix (D matrix order when the
    D matrix is the source; else seg order)."""
    if src == "d_matrix":
        c = cfg["cohorts"][cid]
        cap = rcfg.get("marker_cap_file", "D.cap1.5.mat")
        from .marker_matrix import load_marker_matrix
        return list(load_marker_matrix(c["_abs"] / f"{c['prefix']}.{cap}")
                    ["samples"])
    return sorted(data[cid]["seg"]["sample"].unique())


def _feature_cytoband(cfg: dict[str, Any], data: dict[str, dict],
                      cids: list[str], rcfg: dict, src: str
                      ) -> tuple[np.ndarray, dict[str, list[str]], list[str]]:
    from .pairwise_amplitude import _matrix_for_cohort
    from .score_track import make_grid
    grid = make_grid(rcfg["bin_size"])
    mats = {cid: _matrix_for_cohort(cfg, data, cid, grid, rcfg, src)
            for cid in cids}
    bands, _ = _band_tables(cfg, rcfg, data)
    gstart = grid["start"].to_numpy()
    gend = grid["end"].to_numpy()
    glen = (gend - gstart).astype(float)
    gchrom = grid["chrom"].to_numpy()
    names, cols = [], []
    for b in bands.itertuples():
        band_len = float(b.end - b.start)
        ov = (np.minimum(gend, b.end)
              - np.maximum(gstart, b.start)).astype(float)
        idx = np.flatnonzero((gchrom == b.chrom)
                             & (ov > 0.5 * np.minimum(glen, band_len)))
        if not len(idx):
            continue
        names.append(str(b.band))
        cols.append(idx)
    X = np.hstack([np.stack([
        np.nanmean(mats[cid][cols_], axis=0) for cols_ in cols], axis=0).T
        for cid in cids])                                    # (N, F)
    samples = {cid: _sample_ids(cfg, data, cid, src, rcfg) for cid in cids}
    return X, samples, names


def _feature_peaks(cfg: dict[str, Any], data: dict[str, dict],
                   ca: str, cb: str) -> tuple[np.ndarray, dict[str, list[str]],
                                              list[str]]:
    """Union-peak 0/1/2 calls (rows = samples, NaN/grey -> 0)."""
    from .pairwise_peak import build_onco_data
    q_main = cfg["analysis"]["q_thresholds"][0]
    blocks = build_onco_data(cfg, data, q_main)
    comp = f"{ca}_vs_{cb}"
    blk = next((b for b in blocks if b["comparison"] == comp), None)
    if blk is None or blk["matrix"].empty:
        return np.empty((0, 0)), {ca: [], cb: []}, []
    mat = blk["matrix"]
    cols_a = list(blk["blocks"][0]["columns"])
    cols_b = list(blk["blocks"][1]["columns"])
    X = mat.loc[:, cols_a + cols_b].T.to_numpy(float)         # (N, F)
    X = np.nan_to_num(X, nan=0.0)
    samples = {ca: cols_a, cb: cols_b}
    return X, samples, list(mat.index)


def _band_tables(cfg: dict[str, Any], rcfg: dict, data: dict[str, dict]):
    from .pairwise_regional import _load_band_genes
    return _load_band_genes(cfg, rcfg, data)


def _clean_features(X: np.ndarray, max_nan_frac: float
                    ) -> tuple[np.ndarray, np.ndarray]:
    """Drop features with too much missingness; median-impute; z-score."""
    keep = np.isfinite(X).mean(axis=0) >= (1.0 - max_nan_frac)
    X = X[:, keep]
    if X.shape[1] == 0:
        return X, keep
    med = np.nanmedian(X, axis=0)
    idx = np.where(np.isnan(X))
    X[idx] = np.take(med, idx[1])
    sd = X.std(axis=0)
    sd[sd == 0] = 1.0
    return (X - X.mean(axis=0)) / sd, keep


# ---------------------------------------------------------------------------
# PERMANOVA
# ---------------------------------------------------------------------------

def _permanova(D2: np.ndarray, labels: np.ndarray,
               n_perm: int, rng: np.random.Generator,
               paired_pairs: np.ndarray | None) -> tuple[float, float, float]:
    """Pseudo-F, R^2 and permutation p.  ``paired_pairs`` (n_pairs x 2)
    switches to within-pair label swaps."""
    n = len(labels)
    g = 2
    ss_t = 0.5 * D2.sum() / n

    def _ss_w(lab: np.ndarray) -> float:
        mm = (lab[:, None] == np.array([0, 1])[None, :]).astype(float)
        w = mm.T @ D2 @ mm
        n0 = float((lab == 0).sum())
        n1 = float((lab == 1).sum())
        if n0 < 2 or n1 < 2:
            return np.nan
        return 0.5 * (w[0, 0] / n0 + w[1, 1] / n1)

    ss_w = _ss_w(labels)
    ss_b = ss_t - ss_w
    f_obs = (ss_b / (g - 1)) / (ss_w / (n - g))
    r2 = ss_b / ss_t if ss_t > 0 else np.nan
    ge = 0
    for _ in range(n_perm):
        if paired_pairs is not None:
            lab = labels.copy()
            flip = rng.random(len(paired_pairs)) < 0.5
            lab[paired_pairs[flip, 0]] = 1 - lab[paired_pairs[flip, 0]]
            lab[paired_pairs[flip, 1]] = 1 - lab[paired_pairs[flip, 1]]
        else:
            lab = rng.permutation(labels)
        w = _ss_w(lab)
        if not np.isfinite(w) or w <= 0:
            continue
        f = ((ss_t - w) / (g - 1)) / (w / (n - g))
        ge += int(f >= f_obs)
    p = (1.0 + ge) / (1.0 + n_perm)
    return float(f_obs), float(r2), float(p)


# ---------------------------------------------------------------------------
# classifier
# ---------------------------------------------------------------------------

def _make_clf(mcfg: dict):
    """sklearn >=1.8 deprecates ``penalty``: L2 is the solver default
    (l1_ratio=0); elastic net is selected via a fractional l1_ratio."""
    if mcfg["penalty"] == "elasticnet":
        return LogisticRegression(solver="saga", l1_ratio=mcfg["l1_ratio"],
                                  max_iter=2000, tol=1e-3)
    return LogisticRegression(solver="liblinear", l1_ratio=0)


def _cv_auc(X: np.ndarray, y: np.ndarray, groups: np.ndarray | None,
            mcfg: dict) -> float:
    n = len(y)
    if groups is not None:
        n_splits = min(5, len(np.unique(groups)))
        if n_splits < 2:
            return np.nan
        cv = GroupKFold(n_splits=n_splits)
        pred = cross_val_predict(_make_clf(mcfg), X, y, cv=cv,
                                 groups=groups, method="predict_proba")[:, 1]
    else:
        n_splits = min(5, int(y.sum()), int((1 - y).sum()))
        if n_splits < 2:
            return np.nan
        cv = StratifiedKFold(n_splits=n_splits, shuffle=True,
                             random_state=mcfg["seed"])
        pred = cross_val_predict(_make_clf(mcfg), X, y, cv=cv,
                                 method="predict_proba")[:, 1]
    from sklearn.metrics import roc_auc_score
    return float(roc_auc_score(y, pred))


# ---------------------------------------------------------------------------
# module entry
# ---------------------------------------------------------------------------

def _complete_pairs(sample_list: list[str]
                    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    """Rows belonging to complete cross-cohort pairs (paired designs).

    ``assign_comparison_kinds`` marks a comparison paired at >= ``paired_frac``
    sample overlap, so a "paired" sample list can contain singletons — and
    the restricted within-pair permutations / GroupKFold need exactly one
    sample per cohort per pair.  Samples without a complete pair are
    dropped (with their rows) instead of crashing the pair indexing.

    Returns ``(keep_mask over sample_list, pair_ids per kept row,
    (n_pairs x 2) kept-row index array, dropped sample names)``.
    """
    from collections import Counter
    cnt = Counter(sample_list)
    complete = {s for s, c in cnt.items() if c >= 2}
    keep = np.array([s in complete for s in sample_list], dtype=bool)
    kept = [s for s, k in zip(sample_list, keep) if k]
    ids = sorted(complete)
    pos = {s: i for i, s in enumerate(ids)}
    pair_ids = np.array([pos[s] for s in kept], dtype=int)
    pairs_idx = np.array(
        [[int(np.flatnonzero(pair_ids == i)[0]),
          int(np.flatnonzero(pair_ids == i)[1])]
         for i in range(len(ids))], dtype=int).reshape(-1, 2)
    dropped = sorted(set(sample_list) - complete)
    return keep, pair_ids, pairs_idx, dropped

def run_module_p8(cfg: dict[str, Any], data: dict[str, dict]) -> dict[str, Any]:
    rcfg = _regional_cfg(cfg)
    mcfg = _mv_cfg(cfg, rcfg)
    out_skip = out_table(cfg, "P8_skipped.txt")

    def skip(msg: str) -> dict[str, Any]:
        from pathlib import Path
        Path(out_skip).write_text(msg + "\n", encoding="utf-8")
        print(f"  [P8] skipped: {msg}")
        empty = pd.DataFrame()
        return dict(separation=empty, pca=empty, stability=empty,
                    summary=empty, skipped=msg, auc_perm={})

    if not _HAS_SK:
        return skip("scikit-learn not importable (pip install scikit-learn)")
    if mcfg["features"] not in ("grid", "cytoband", "peaks"):
        return skip(f"unknown feature set: {mcfg['features']!r}")

    sep_rows, pca_rows, stab_rows, sum_rows = [], [], [], []
    auc_perm_store: dict[str, np.ndarray] = {}

    # one shared matrix source for the whole run (D matrix when every
    # cohort has one, else seg-derived) — never mix within a run
    from .pairwise_amplitude import _matrix_source
    all_cids = list(cfg["cohorts"])
    src = _matrix_source(cfg, data, set(all_cids), rcfg)

    for cmp_ in cfg["comparisons"]:
        ca, cb = cmp_["a"], cmp_["b"]
        comp = f"{ca}_vs_{cb}"
        kind = cmp_.get("kind", "unpaired")
        if kind == "paired":
            print(f"  [P8] paired comparison {comp}: within-pair "
                  "permutation / GroupKFold")
        # shared sample ids in an unpaired comparison would create dependence
        shared = set()
        if kind != "paired":
            sa = set(data[ca]["seg"]["sample"].unique())
            sb = set(data[cb]["seg"]["sample"].unique())
            shared = sa & sb
            if shared:
                print(f"  [P8] {comp}: dropping {len(shared)} shared sample "
                      "ids (unpaired comparison)")

        # ---- features ----
        if mcfg["features"] == "peaks":
            X, samples, feat_names = _feature_peaks(cfg, data, ca, cb)
            X_all = X
            sample_list = list(samples.get(ca, [])) + list(samples.get(cb, []))
            cohorts = ([ca] * len(samples.get(ca, []))
                       + [cb] * len(samples.get(cb, [])))
        else:
            X_all, samples, feat_names = (
                _feature_grid(cfg, data, [ca, cb], rcfg, src)
                if mcfg["features"] == "grid"
                else _feature_cytoband(cfg, data, [ca, cb], rcfg, src))
            sample_list = list(samples[ca]) + list(samples[cb])
            cohorts = [ca] * len(samples[ca]) + [cb] * len(samples[cb])
        if not sample_list or X_all.size == 0:
            sep_rows.append(dict(comparison=comp, kind=kind,
                                 feature_set=mcfg["features"],
                                 status="no_features"))
            continue
        keep_samples = np.array([s not in shared for s in sample_list])
        X_all = X_all[keep_samples]
        sample_list = [s for s, k in zip(sample_list, keep_samples) if k]
        cohorts = [c for c, k in zip(cohorts, keep_samples) if k]
        y = np.array([1 if c == cb else 0 for c in cohorts])
        pair_ids: np.ndarray | None = None
        pairs_idx: np.ndarray | None = None
        if kind == "paired":
            keep, pair_ids, pairs_idx, dropped = _complete_pairs(sample_list)
            if dropped:
                print(f"  [P8] {comp}: {len(dropped)} sample(s) without a "
                      "complete cross-cohort pair dropped from the paired "
                      "analysis")
            X_all = X_all[keep]
            y = y[keep]
            sample_list = [s for s, k in zip(sample_list, keep) if k]
            cohorts = [c for c, k in zip(cohorts, keep) if k]
        X, kept_feat = _clean_features(X_all, mcfg["max_nan_frac"])
        feat_names = [f for f, k in zip(feat_names, kept_feat) if k]
        na, nb = int((y == 0).sum()), int((y == 1).sum())
        if na < mcfg["min_n"] or nb < mcfg["min_n"]:
            sep_rows.append(dict(comparison=comp, kind=kind,
                                 feature_set=mcfg["features"], n_a=na, n_b=nb,
                                 status=f"min_n={mcfg['min_n']} not met"))
            continue

        # ---- PERMANOVA ----
        rng = np.random.default_rng(mcfg["seed"])
        D2 = squareform(pdist(X, metric="euclidean")) ** 2
        f_obs, r2, p_perm = _permanova(D2, y, mcfg["n_perm"], rng, pairs_idx)

        # ---- PCA ----
        Xc = X - X.mean(axis=0)
        try:
            u, s, vt = np.linalg.svd(Xc, full_matrices=False)
            var = s ** 2
            var_pct = var / max(var.sum(), 1e-12)
            coords = Xc @ vt[:2].T
            for i, (samp, co, cid) in enumerate(zip(sample_list, coords,
                                                    cohorts)):
                pca_rows.append(dict(
                    comparison=comp, kind=kind, sample=samp, cohort=cid,
                    pair_id=(int(pair_ids[i]) if pair_ids is not None else ""),
                    pc1=float(co[0]), pc2=float(co[1]),
                    var_pc1=round(float(var_pct[0]), 4),
                    var_pc2=round(float(var_pct[1]), 4)))
        except np.linalg.LinAlgError:
            var_pct = (np.nan, np.nan)

        # ---- CV AUC + permutation ----
        auc = _cv_auc(X, y, pair_ids, mcfg)
        aucs_perm = np.full(mcfg["n_auc_perm"], np.nan)
        for b in range(mcfg["n_auc_perm"]):
            if kind == "paired" and pairs_idx is not None:
                yb = y.copy()
                flip = rng.random(len(pairs_idx)) < 0.5
                yb[pairs_idx[flip, 0]] = 1 - yb[pairs_idx[flip, 0]]
                yb[pairs_idx[flip, 1]] = 1 - yb[pairs_idx[flip, 1]]
            else:
                yb = rng.permutation(y)
            try:
                aucs_perm[b] = _cv_auc(X, yb, pair_ids, mcfg)
            except Exception:
                aucs_perm[b] = np.nan
        ok = np.isfinite(aucs_perm)
        p_auc = ((1.0 + int((aucs_perm[ok] >= auc).sum())) / (1.0 + int(ok.sum()))
                 if np.isfinite(auc) and ok.any() else np.nan)
        auc_perm_store[comp] = aucs_perm

        # ---- stability selection ----
        sel = np.zeros(X.shape[1])
        coefs = np.zeros(X.shape[1])
        k = min(mcfg["top_k"], X.shape[1])
        for _ in range(mcfg["n_boot"]):
            idx = []
            for cls in (0, 1):
                pool = np.flatnonzero(y == cls)
                idx.append(rng.choice(pool, size=len(pool), replace=True))
            boot = np.concatenate(idx)
            try:
                clf = _make_clf(mcfg).fit(X[boot], y[boot])
            except Exception:
                continue
            coef = np.abs(clf.coef_).ravel()
            coefs += clf.coef_.ravel()
            top = np.argsort(-coef)[:k]
            sel[top] += 1
        sel_freq = sel / max(mcfg["n_boot"], 1)
        order = np.argsort(-sel_freq)[:mcfg["top_k"]]
        annot = _annotate_features(cfg, rcfg, data, mcfg["features"],
                                   feat_names)
        for rank, fi in enumerate(order, start=1):
            stab_rows.append(dict(
                comparison=comp, rank=rank, feature=feat_names[fi],
                annotation=annot.get(feat_names[fi], ""),
                selection_freq=round(float(sel_freq[fi]), 3),
                mean_coef=round(float(coefs[fi] / max(mcfg["n_boot"], 1)), 5)))

        sep_rows.append(dict(
            comparison=comp, kind=kind, feature_set=mcfg["features"],
            status="ok", n_features=X.shape[1], n_a=na, n_b=nb,
            permanova_F=round(f_obs, 4), permanova_R2=round(r2, 4),
            permanova_p=p_perm, n_perm=mcfg["n_perm"],
            auc_oof=(round(auc, 4) if np.isfinite(auc) else np.nan),
            auc_perm_p=p_auc, n_auc_perm=int(ok.sum()),
            var_pc1=round(float(var_pct[0]), 4)
            if np.isfinite(var_pct[0]) else np.nan,
            var_pc2=round(float(var_pct[1]), 4)
            if np.isfinite(var_pct[1]) else np.nan,
        ))
        sum_rows.append(dict(comparison=comp, kind=kind,
                             feature_set=mcfg["features"], n_features=X.shape[1],
                             permanova_R2=round(r2, 4), permanova_p=p_perm,
                             auc_oof=(round(auc, 4)
                                      if np.isfinite(auc) else np.nan),
                             auc_perm_p=p_auc,
                             top_feature=(stab_rows[-1]["feature"]
                                          if stab_rows else "")))
        print(f"  [P8] {comp}: F={f_obs:.2f} R2={r2:.3f} "
              f"p={p_perm:.3f} | AUC={auc:.3f} (perm p={p_auc:.3f})")

    sep = pd.DataFrame(sep_rows)
    pca = pd.DataFrame(pca_rows)
    stab = pd.DataFrame(stab_rows)
    summary = pd.DataFrame(sum_rows)
    sep.to_csv(out_table(cfg, "P8_separation.tsv"), sep="\t", index=False)
    pca.to_csv(out_table(cfg, "P8_pca_coords.tsv"), sep="\t", index=False)
    stab.to_csv(out_table(cfg, "P8_stability.tsv"), sep="\t", index=False)
    summary.to_csv(out_table(cfg, "P8_summary.tsv"), sep="\t", index=False)
    if out_skip.exists():
        out_skip.unlink()
    return dict(separation=sep, pca=pca, stability=stab, summary=summary,
                skipped=None, auc_perm=auc_perm_store, mcfg=mcfg)


def _annotate_features(cfg: dict[str, Any], rcfg: dict, data: dict[str, dict],
                       feature_set: str, names: list[str]) -> dict[str, str]:
    """Cytoband annotation per feature (best effort)."""
    if feature_set == "cytoband":
        return {n: n for n in names}
    out: dict[str, str] = {}
    bands, _ = _band_tables(cfg, rcfg, data)
    if bands.empty:
        return out
    from .config import CHROMS
    bchrom = bands["chrom"].to_numpy()
    bstart = bands["start"].to_numpy()
    bend = bands["end"].to_numpy()
    bname = bands["band"].astype(str).to_numpy()
    for n in names:
        try:
            chrom, span = n.split(":")
            s = int(span.split("-")[0])
        except ValueError:
            continue
        if chrom not in CHROMS:
            continue
        i = np.flatnonzero((bchrom == chrom) & (bstart <= s) & (s < bend))
        if len(i):
            out[n] = str(bname[i[0]])
    return out
