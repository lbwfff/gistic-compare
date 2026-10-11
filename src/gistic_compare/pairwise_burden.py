"""Module P9 (pairwise): three-level SCNA burden score + wGII.

Adopted from the SCNA-disparity methodology of Yuan et al., Cancer Cell
2018 (reproduction manual R05 §2/§3): a per-sample burden score that counts
copy-number events at three genomic scales **separately** — focal, arm,
chromosome — weights each event by its amplitude, rank-normalizes each
level within the pooled two-cohort sample set, and compares the levels and
their sum between cohorts.  Complemented by wGII, the GISTIC-independent
weighted genome instability index (Burrell et al. 2013).

Why this module exists: a cohort whose genomes are globally more
rearranged shows higher alteration frequency at *every* locus, which
inflates every per-locus comparison (P3-P8).  The burden panel makes that
background explicit, splits it into the focal vs broad scales (the
"broad vs focal decomposition"), and gives every downstream table the
context it needs to be read against.

Per-sample quantities
---------------------
S_focal   sum of |C| over focal events.  Focal events come from GISTIC's
          own marker matrix (``D.cap1.5.mat``; seg fallback): markers are
          thresholded at +/- ``focal_thr`` (default 0.25), adjacent markers
          of equal threshold class are merged into intervals, and intervals
          overlapping a same-direction broad event by more than
          ``focal_broad_overlap_max`` are dropped (no double counting).
S_arm     sum of |C| over arm events: > ``arm_min_frac`` (0.70) of the
          arm's physical length covered by same-direction segments beyond
          +/- ``arm_thr``.
S_chrom   same, after merging both arms of a chromosome into one event
          when they share sign and their length-weighted log2 values
          differ by < ``arm_match_diff`` (0.1).
C         amplitude weight of an event with log2 ratio c:
          2 (|c| >= 1), 1 (0.25 <= |c| < 1), 0 (noise band), symmetric
          negative for deletions.  Levels are summed over |C| so
          amplifications and deletions count equally.
S         S_focal* + S_arm* + S_chrom* where S_* are ranks / n computed on
          the **pooled** two-cohort sample set (rank-normalizing within a
          cohort would destroy the between-cohort signal by construction).
          S in (0, 3]; its interpretation is *relative positioning*, so the
          unnormalized raw sums are reported alongside.
wGII      per chromosome, the length fraction with |CN - ploidy| > 0.3
          (copy-number units, NOT log2), averaged over the 22 autosomes;
          ploidy = length-weighted median copy number.

Tests: per comparison and metric, OLS ``metric ~ cohort`` (beta = change
per unit of being cohort b, with 95% CI) plus Mann-Whitney and the
rank-biserial effect size.  Optional covariates from the metadata table
are added one model with all covariates (documented as
``covariates=...``); without metadata the comparison is unadjusted and
clinical confounding can only be bounded (module CF).

Outputs: P9_burden_scores.tsv, P9_burden_tests.tsv, P9_crosscheck.tsv,
         P9_summary.tsv
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from .config import CHROMS, out_table
from .pairwise_amplitude import _matrix_source
from .pairwise_arm import _arm_intervals
from .pairwise_model import _load_metadata
from .pairwise_regional import _regional_cfg


def _burden_cfg(cfg: dict[str, Any]) -> dict:
    b = (cfg.get("analysis") or {}).get("burden") or {}
    return dict(
        arm_thr=float(b.get("arm_thr", 0.25)),
        arm_min_frac=float(b.get("arm_min_frac", 0.70)),
        arm_match_diff=float(b.get("arm_match_diff", 0.10)),
        focal_thr=float(b.get("focal_thr", 0.25)),
        high_thr=float(b.get("high_thr", 1.0)),
        focal_broad_overlap_max=float(b.get("focal_broad_overlap_max", 0.5)),
        sensitivity=bool(b.get("sensitivity", False)),
    )


# ---------------------------------------------------------------------------
# amplitude weighting
# ---------------------------------------------------------------------------

def _threshold_c(c: np.ndarray, thr: float, high: float) -> np.ndarray:
    """Amplitude weight C of log2 ratio c (Yuan et al. 2018):

    2 (c >= high), 1 (thr <= c < high), 0 (-thr < c < thr),
    -1 (-high <= c <= -thr), -2 (c < -high).  Sign is kept; the level sums
    use |C| so amplifications and deletions count equally.
    """
    c = np.asarray(c, dtype=float)
    return np.select(
        [c >= high, c >= thr, c > -thr, c >= -high],
        [2.0, 1.0, 0.0, -1.0], default=-2.0)


# ---------------------------------------------------------------------------
# arm / chromosome level events (from tumors.seg)
# ---------------------------------------------------------------------------

def _arm_events(seg: pd.DataFrame, arms: pd.DataFrame, thr: float,
                min_frac: float) -> pd.DataFrame:
    """Arm-level broad events per sample.

    An arm counts when > ``min_frac`` of its **physical length** is covered
    by same-direction segments beyond +/- ``thr``.  Returns one row per
    event: cohort, sample, arm, chrom, direction, mean_log2 (length-weighted
    mean of the qualifying segments), frac, and the qualifying spans
    (broad_spans) used later for focal-event exclusion.
    """
    samples = sorted(seg["sample"].unique())
    sidx = {s: i for i, s in enumerate(samples)}
    arm_chrom = arms["chrom"].to_numpy()
    arm_start = arms["start"].to_numpy()
    arm_end = arms["end"].to_numpy()
    arm_len = np.maximum(arm_end - arm_start, 1.0).astype(float)
    n_arm = len(arms)
    n_s = len(samples)
    # length covered and value sum per (arm, sample) and direction
    cov = {"Amp": np.zeros((n_arm, n_s)), "Del": np.zeros((n_arm, n_s))}
    vsum = {"Amp": np.zeros((n_arm, n_s)), "Del": np.zeros((n_arm, n_s))}
    spans: dict[tuple[int, str], list[tuple[int, int]]] = {}
    for c in CHROMS:
        sel = np.flatnonzero(arm_chrom == c)
        if len(sel) == 0:
            continue
        sub = seg[seg["chrom"] == c]
        if sub.empty:
            continue
        a_s = arm_start[sel].astype(float)
        a_e = arm_end[sel].astype(float)
        for r in sub.itertuples():
            lo = np.maximum(a_s, r.start)
            hi = np.minimum(a_e, r.end)
            ov = np.clip(hi - lo, 0.0, None)
            if not ov.any():
                continue
            j = sidx[r.sample]
            for dirn, qual in (("Amp", r.seg_mean >= thr),
                               ("Del", r.seg_mean <= -thr)):
                if qual:
                    cov[dirn][sel, j] += ov
                    vsum[dirn][sel, j] += ov * r.seg_mean
                    idx = sel[ov > 0]
                    for ai in idx:
                        spans.setdefault((ai, dirn), []).append(
                            (int(max(r.start, arm_start[ai])),
                             int(min(r.end, arm_end[ai]))))
    rows = []
    for ai, arm_row in enumerate(arms.itertuples()):
        for dirn in ("Amp", "Del"):
            for si, s in enumerate(samples):
                ln = cov[dirn][ai, si]
                frac = ln / arm_len[ai] if arm_len[ai] else np.nan
                if np.isfinite(frac) and frac > min_frac and ln > 0:
                    rows.append(dict(
                        cohort="", sample=s, arm=str(arm_row.arm),
                        chrom=arm_row.chrom, start=int(arm_row.start),
                        end=int(arm_row.end), direction=dirn,
                        mean_log2=float(vsum[dirn][ai, si] / ln),
                        frac=float(frac),
                        spans=spans.get((ai, dirn), [])))
    return pd.DataFrame(rows, columns=[
        "cohort", "sample", "arm", "chrom", "start", "end", "direction",
        "mean_log2", "frac", "spans"])


def _merge_chromosome_events(arm_ev: pd.DataFrame,
                             match_diff: float) -> pd.DataFrame:
    """Merge both arms of a chromosome into one chromosome-level event when
    they share direction and their mean log2 values differ by < match_diff.

    Returns the **remaining** arm events (not merged) and the chromosome
    events as a tuple (arm_ev, chrom_ev); the chromosome event's value is
    the coverage-weighted mean of the two arms."""
    arm_ev = arm_ev.copy()
    arm_ev["_c"] = arm_ev.index
    chrom_rows = []
    merged_ids: set[int] = set()
    for (sample, chrom, dirn), g in arm_ev.groupby(
            ["sample", "chrom", "direction"], sort=False):
        if len(g) < 2:
            continue
        lo = g.loc[g["arm"].str.endswith("p")]
        hi = g.loc[g["arm"].str.endswith("q")]
        if lo.empty or hi.empty:
            continue
        r_p, r_q = lo.iloc[0], hi.iloc[0]
        if abs(r_p["mean_log2"] - r_q["mean_log2"]) < match_diff:
            w_p, w_q = r_p["frac"], r_q["frac"]
            mean_log2 = ((r_p["mean_log2"] * w_p + r_q["mean_log2"] * w_q)
                         / max(w_p + w_q, 1e-12))
            chrom_rows.append(dict(
                sample=sample, chrom=chrom,
                start=int(min(r_p["start"], r_q["start"])),
                end=int(max(r_p["end"], r_q["end"])),
                direction=dirn, mean_log2=float(mean_log2),
                arms=f"{r_p['arm']}+{r_q['arm']}"))
            merged_ids.update([r_p["_c"], r_q["_c"]])
    arm_out = arm_ev[~arm_ev.index.isin(merged_ids)].drop(columns="_c")
    chrom_ev = pd.DataFrame(chrom_rows, columns=[
        "sample", "chrom", "start", "end", "direction", "mean_log2", "arms"])
    return arm_out, chrom_ev


# ---------------------------------------------------------------------------
# focal level events
# ---------------------------------------------------------------------------

def _merge_focal_runs(events: list[dict]) -> list[dict]:
    """Merge adjacent same-chromosome events with equal C (in-place list
    assumed sorted by (chrom, start); touching intervals merge)."""
    out: list[dict] = []
    for ev in sorted(events, key=lambda e: (e["chrom"], e["start"])):
        if out:
            last = out[-1]
            if (last["chrom"] == ev["chrom"] and last["C"] == ev["C"]
                    and ev["start"] <= last["end"] + 1):
                last["end"] = max(last["end"], ev["end"])
                continue
        out.append(dict(ev))
    return out


def _focal_from_segments(segs: pd.DataFrame, thr: float, high: float
                         ) -> list[dict]:
    """Focal events from one sample's (chrom, start, end, value) segments:
    threshold -> merge adjacent equal-C runs -> non-zero runs are events."""
    segs = segs.sort_values(["chrom", "start"])
    chroms = segs["chrom"].to_numpy()
    cvals = _threshold_c(segs["value"].to_numpy(float), thr, high)
    events: list[dict] = []
    run_chrom: str | None = None
    run_c = 0.0
    run: list[int] = []

    def flush() -> None:
        if run_c != 0.0 and run:
            sub = segs.iloc[run]
            events.append(dict(chrom=str(run_chrom),
                               start=int(sub["start"].min()),
                               end=int(sub["end"].max()), C=float(run_c)))

    for i in range(len(segs)):
        c = float(cvals[i])
        contiguous = bool(run) and chroms[i] == run_chrom and \
            segs.iloc[i]["start"] <= segs.iloc[run[-1]]["end"] + 1
        if c == run_c and contiguous:
            run.append(i)
        else:
            flush()
            run_chrom, run_c, run = chroms[i], c, [i]
    flush()
    return events


def _drop_broad_covered(events: list[dict], broad_spans: list[dict],
                        overlap_max: float) -> list[dict]:
    """Drop focal intervals overlapping a same-direction broad event by more
    than ``overlap_max`` of their own length (no double counting)."""
    if not broad_spans:
        return events
    out = []
    for ev in events:
        ev_dir = "Amp" if ev["C"] > 0 else "Del"
        covered = 0.0
        length = max(ev["end"] - ev["start"], 1)
        for sp in broad_spans:
            if sp["chrom"] != ev["chrom"] or sp["direction"] != ev_dir:
                continue
            ov = (min(ev["end"], sp["end"]) - max(ev["start"], sp["start"]))
            if ov > 0:
                covered += ov
        if covered / length <= overlap_max:
            out.append(ev)
    return out


def _d_matrix_segments(mat_path, rcfg: dict) -> pd.DataFrame:
    """Per-sample segments of GISTIC's marker matrix as
    (sample, chrom, start, end, value), split at chromosome boundaries."""
    from .marker_matrix import load_marker_matrix
    mm = load_marker_matrix(mat_path)
    pos = mm["pos"]
    chrom = mm["chrom"]
    M = len(pos)
    # chromosome block boundaries in marker space (markers are already in
    # genomic order; adjacent-equality gives the blocks without assuming
    # any particular chromosome naming order)
    new_block = np.r_[True, chrom[1:] != chrom[:-1]]
    block_start = np.flatnonzero(new_block)
    block_id = np.cumsum(new_block) - 1
    n_blocks = len(block_start)
    segs = []
    ir, jc, segval = mm["ir"], mm["jc"], mm["segval"]
    for j, s in enumerate(mm["samples"]):
        a, b = int(jc[j]), int(jc[j + 1])
        if b <= a:
            continue
        rows = ir[a:b]
        vals = segval[a:b]
        for k in range(len(rows)):
            r0 = int(rows[k])
            r1 = int(rows[k + 1]) if k + 1 < len(rows) else M
            # split at chromosome boundaries
            bi = block_id[r0]
            cuts = [r0]
            while bi < n_blocks - 1 and block_start[bi + 1] < r1:
                cuts.append(int(block_start[bi + 1]))
                bi += 1
            cuts.append(r1)
            for c0, c1 in zip(cuts[:-1], cuts[1:]):  # noqa: B905
                if c1 <= c0:
                    continue
                segs.append(dict(
                    sample=s, chrom=str(chrom[c0]),
                    start=int(pos[c0]),
                    end=int(pos[c1 - 1] + (pos[1] - pos[0] if M > 1 else 1)),
                    value=float(segval[a + k])))
    return pd.DataFrame(segs, columns=["sample", "chrom", "start", "end",
                                       "value"])


# ---------------------------------------------------------------------------
# wGII and FGA
# ---------------------------------------------------------------------------

def _wgii_fga(seg: pd.DataFrame, cn_diff: float = 0.3,
              fga_thr: float = 0.2) -> pd.DataFrame:
    """Per-sample ploidy, wGII (copy-number units!) and FGA (log2 units)."""
    rows = []
    for s, g in seg.groupby("sample", sort=True):
        g = g[g["chrom"].isin(CHROMS)].dropna(
            subset=["seg_mean", "start", "end"])
        if g.empty:
            rows.append(dict(sample=s, ploidy=np.nan, wGII=np.nan,
                             fga=np.nan))
            continue
        length = (g["end"] - g["start"]).clip(lower=0).to_numpy(float)
        cn = 2.0 * np.power(2.0, g["seg_mean"].to_numpy(float))
        # length-weighted median copy number = ploidy
        o = np.argsort(cn)
        cum = np.cumsum(length[o])
        total = cum[-1]
        ploidy = float(cn[o][np.searchsorted(cum, 0.5 * total)]) \
            if total > 0 else np.nan
        gxs = []
        for c in CHROMS:
            m = (g["chrom"] == c).to_numpy()
            ln = length[m].sum()
            if ln <= 0:
                continue
            gxs.append(length[m][np.abs(cn[m] - ploidy) > cn_diff].sum() / ln)
        rows.append(dict(
            sample=s, ploidy=ploidy,
            wGII=float(np.mean(gxs)) if gxs else np.nan,
            fga=float(length[np.abs(g["seg_mean"].to_numpy(float))
                             >= fga_thr].sum() / length.sum())))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# statistics
# ---------------------------------------------------------------------------

def _ols_cohort(y: np.ndarray, is_b: np.ndarray,
                cov: np.ndarray | None) -> dict:
    """OLS y ~ 1 + is_b (+ covariates); returns beta/CI/p for the cohort
    coefficient."""
    import statsmodels.api as sm
    cols = [np.ones(len(y)), is_b.astype(float)]
    if cov is not None and cov.size:
        cols.extend(c for c in cov.T)
    X = np.column_stack(cols)
    try:
        fit = sm.OLS(y, X).fit()
    except Exception:
        return dict(beta=np.nan, ci_lo=np.nan, ci_hi=np.nan, p_ols=np.nan)
    ci = fit.conf_int()           # numpy array rows = params
    return dict(beta=float(fit.params[1]), ci_lo=float(ci[1, 0]),
                ci_hi=float(ci[1, 1]), p_ols=float(fit.pvalues[1]))


def _rank_biserial(xa: np.ndarray, xb: np.ndarray) -> float:
    try:
        u1 = stats.mannwhitneyu(xa, xb, alternative="two-sided")[0]
    except ValueError:
        return np.nan
    return 2.0 * u1 / (len(xa) * len(xb)) - 1.0


# ---------------------------------------------------------------------------
# module entry
# ---------------------------------------------------------------------------

METRICS = [
    ("S_focal_star", "S focal (percentile)"),
    ("S_arm_star", "S arm (percentile)"),
    ("S_chrom_star", "S chromosome (percentile)"),
    ("S_total", "SCNA score S (0-3)"),
    ("wGII", "wGII"),
]


def _cohort_scores(cfg: dict[str, Any], data: dict[str, dict],
                   arms: pd.DataFrame, bcfg: dict, src: str, rcfg: dict,
                   segs_pre: dict | None = None, verbose: bool = True
                   ) -> pd.DataFrame:
    """Per-sample raw burden metrics for every cohort (no ranking)."""
    score_frames = []
    segs_cache = dict(segs_pre or {})
    for cid, d in data.items():
        seg = d["seg"]
        c = cfg["cohorts"][cid]
        arm_ev = _arm_events(seg, arms, bcfg["arm_thr"], bcfg["arm_min_frac"])
        arm_ev["cohort"] = cid
        arm_out, chrom_ev = _merge_chromosome_events(
            arm_ev, bcfg["arm_match_diff"])
        broad_spans = [dict(sample=r["sample"], chrom=r["chrom"],
                            direction=r["direction"], start=r["start"],
                            end=r["end"])
                       for r in arm_out.to_dict("records")]
        broad_spans += [dict(sample=r["sample"], chrom=r["chrom"],
                             direction=r["direction"], start=r["start"],
                             end=r["end"])
                        for r in chrom_ev.to_dict("records")]
        if cid in segs_cache:
            segs_all = segs_cache[cid]
        elif src == "d_matrix":
            cap = rcfg.get("marker_cap_file", "D.cap1.5.mat")
            mat = c["_abs"] / f"{c['prefix']}.{cap}"
            segs_all = _d_matrix_segments(mat, rcfg)
            segs_cache[cid] = segs_all
        else:
            segs_all = seg.rename(columns={"seg_mean": "value"})[
                ["sample", "chrom", "start", "end", "value"]]
        arm_out["C"] = _threshold_c(arm_out["mean_log2"].to_numpy(float),
                                    bcfg["focal_thr"], bcfg["high_thr"])
        chrom_ev["C"] = _threshold_c(chrom_ev["mean_log2"].to_numpy(float),
                                     bcfg["focal_thr"], bcfg["high_thr"])
        s_arm = arm_out.groupby("sample")["C"].apply(lambda x: np.abs(x).sum())
        s_chrom = chrom_ev.groupby("sample")["C"].apply(
            lambda x: np.abs(x).sum())
        samples = sorted(seg["sample"].unique())
        focal_rows = []
        for s in samples:
            segs = segs_all[segs_all["sample"] == s]
            spans = [sp for sp in broad_spans if sp["sample"] == s]
            ev = _focal_from_segments(segs, bcfg["focal_thr"],
                                      bcfg["high_thr"])
            ev = _drop_broad_covered(ev, spans,
                                     bcfg["focal_broad_overlap_max"])
            focal_rows.append(dict(
                sample=s, n_focal_events=len(ev),
                s_focal=float(np.abs([e["C"] for e in ev]).sum())
                if ev else 0.0))
        focal_df = pd.DataFrame(focal_rows)
        wg = _wgii_fga(seg)
        sc = pd.DataFrame({"sample": samples})
        sc["S_focal"] = sc["sample"].map(
            focal_df.set_index("sample")["s_focal"]).fillna(0.0)
        sc["n_focal_events"] = sc["sample"].map(
            focal_df.set_index("sample")["n_focal_events"]).fillna(0).astype(int)
        sc["S_arm"] = sc["sample"].map(s_arm).fillna(0.0)
        sc["n_arm_events"] = sc["sample"].map(
            arm_out.groupby("sample").size()).fillna(0).astype(int)
        sc["S_chrom"] = sc["sample"].map(s_chrom).fillna(0.0)
        sc["n_chrom_events"] = sc["sample"].map(
            chrom_ev.groupby("sample").size()).fillna(0).astype(int)
        sc = sc.merge(wg, on="sample", how="left")
        sc.insert(0, "cohort", cid)
        score_frames.append(sc)
        if verbose:
            print(f"  [P9] {cid}: focal median={sc['S_focal'].median():.0f}  "
                  f"arm median={sc['S_arm'].median():.0f}  "
                  f"chrom median={sc['S_chrom'].median():.0f}  "
                  f"median wGII={sc['wGII'].median():.3f}")
    return pd.concat(score_frames, ignore_index=True)


def _rank_and_test(cfg: dict[str, Any], scores: pd.DataFrame,
                   meta_df: pd.DataFrame | None, verbose: bool = True
                   ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per-comparison pooled rank normalization + cohort OLS / MW.

    Percentile ranks are computed **within each comparison's own pooled
    sample set** and returned as a per-comparison long table (one row per
    comparison x sample).  Writing ranks back into shared columns would
    overwrite them comparison by comparison: with overlapping comparisons
    (A_vs_B, A_vs_C) cohort A's percentile would come from the A+B pool
    while cohort B's came from the same pool — or worse, A's would be
    re-ranked within A+C while B's stayed ranked within A+B, and the two
    groups would be compared on different scales.
    """
    cov_names: list[str] = []
    if meta_df is not None and len(meta_df.columns) > 1:
        cov_names = [c for c in meta_df.columns if c != "sample"]

    pool_frames: list[pd.DataFrame] = []
    test_rows = []
    for cmp_ in cfg["comparisons"]:
        ca, cb = cmp_["a"], cmp_["b"]
        if cmp_.get("kind") == "paired":
            if verbose:
                print(f"  [P9] {ca}_vs_{cb}: paired comparison skipped — P9 "
                      "implements unpaired tests only (paired inference is "
                      "on the roadmap)")
            continue
        comp = f"{ca}_vs_{cb}"
        pool = scores[scores["cohort"].isin([ca, cb])].copy()
        n_pool = len(pool)
        for lvl in ("S_focal", "S_arm", "S_chrom"):
            pool[f"{lvl}_star"] = (pool[lvl].rank(method="average")
                                   / max(n_pool, 1))
        pool["S_total"] = pool[["S_focal_star", "S_arm_star",
                                "S_chrom_star"]].sum(axis=1, min_count=3)
        pool.insert(0, "comparison", comp)
        pool_frames.append(pool)

        model = pool.copy()
        model["is_b"] = (model["cohort"] == cb).astype(float)
        cov = None
        cov_note = "none"
        if meta_df is not None and cov_names:
            m = model.merge(meta_df, on="sample", how="inner")
            if len(m) >= 10 and all(c in m.columns for c in cov_names):
                dm = pd.get_dummies(m[cov_names].astype("string"),
                                    drop_first=True, dtype=float)
                dm = dm.fillna(0.0)
                cov = dm.to_numpy(float)
                model = m
                cov_note = "+".join(cov_names)
        for metric, label in METRICS:
            d = model.dropna(subset=[metric])
            ya = d.loc[d["cohort"] == ca, metric].to_numpy(float)
            yb = d.loc[d["cohort"] == cb, metric].to_numpy(float)
            if len(ya) < 3 or len(yb) < 3:
                continue
            ols = _ols_cohort(d[metric].to_numpy(float),
                              d["is_b"].to_numpy(float), cov)
            try:
                p_mw = float(stats.mannwhitneyu(ya, yb,
                                                alternative="two-sided")[1])
            except ValueError:
                p_mw = np.nan
            r_rb = _rank_biserial(ya, yb)
            higher = ca if np.nanmedian(ya) > np.nanmedian(yb) else cb
            test_rows.append(dict(
                comparison=comp, metric=metric, label=label, n_a=len(ya),
                n_b=len(yb), covariates=cov_note,
                beta=round(ols["beta"], 4),
                ci_lo=round(ols["ci_lo"], 4), ci_hi=round(ols["ci_hi"], 4),
                p_ols=ols["p_ols"], p_mw=p_mw, r_rb=round(r_rb, 4),
                median_a=round(float(np.median(ya)), 4),
                median_b=round(float(np.median(yb)), 4),
                cohort_higher=higher))
    pooled = (pd.concat(pool_frames, ignore_index=True) if pool_frames
              else pd.DataFrame())
    return pooled, pd.DataFrame(test_rows)


# reproduction-specific parameter grid (R05 §13 #4-#8); 1.0 for the broad
# overlap cap = no exclusion
_SENS_GRID = {
    "arm_min_frac": [0.5, 0.7, 0.9],
    "arm_match_diff": [0.05, 0.10, 0.20],
    "focal_thr": [0.20, 0.25, 0.30],
    "focal_broad_overlap_max": [0.30, 0.50, 1.00],
}
_SENS_METRICS = ("S_total", "S_focal_star", "S_arm_star", "S_chrom_star",
                 "wGII")


def run_sensitivity(cfg: dict[str, Any], data: dict[str, dict],
                    arms: pd.DataFrame, src: str, rcfg: dict,
                    base_bcfg: dict, meta_df: pd.DataFrame | None
                    ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Full parameter grid over the reproduction-specific choices: does the
    cohort-effect conclusion survive the parameterisation?  Returns
    (per-combo rows, per-metric verdict)."""
    import itertools

    names = list(_SENS_GRID)
    combos = list(itertools.product(*[_SENS_GRID[n] for n in names]))
    segs_pre: dict[str, pd.DataFrame] = {}
    if src == "d_matrix":
        cap = rcfg.get("marker_cap_file", "D.cap1.5.mat")
        for cid, c in cfg["cohorts"].items():
            mat = c["_abs"] / f"{c['prefix']}.{cap}"
            segs_pre[cid] = _d_matrix_segments(mat, rcfg)

    rows = []
    for k, combo in enumerate(combos, start=1):
        bcfg = dict(base_bcfg, **dict(zip(names, combo)))  # noqa: B905
        try:
            scores = _cohort_scores(cfg, data, arms, bcfg, src, rcfg,
                                    segs_pre=segs_pre, verbose=False)
            _, tests = _rank_and_test(cfg, scores, meta_df, verbose=False)
        except Exception as e:  # a combo may yield no events at all
            print(f"  [P9-sens] combo {combo} failed: {e}")
            continue
        for metric in _SENS_METRICS:
            r = tests[tests["metric"] == metric]
            if len(r):
                rr = r.iloc[0]
                rows.append({**dict(zip(names, combo)),  # noqa: B905
                             "metric": metric, "beta": rr["beta"],
                             "p_ols": rr["p_ols"],
                             "cohort_higher": rr["cohort_higher"]})
        if k % 10 == 0:
            print(f"  [P9-sens] {k}/{len(combos)} combos done")

    sens = pd.DataFrame(rows)
    verdict_rows = []
    for metric in _SENS_METRICS:
        st = sens[sens["metric"] == metric]
        if st.empty:
            continue
        signs = np.sign(st["beta"].to_numpy(float))
        n_pos, n_neg = int((signs > 0).sum()), int((signs < 0).sum())
        verdict_rows.append(dict(
            metric=metric, n_combos=len(st),
            n_beta_positive=n_pos, n_beta_negative=n_neg,
            sign_consistent=bool(n_pos == 0 or n_neg == 0),
            p_min=round(float(st["p_ols"].min()), 5),
            p_max=round(float(st["p_ols"].max()), 5),
            n_p_lt_05=int((st["p_ols"] < 0.05).sum()),
        ))
    return sens, pd.DataFrame(verdict_rows)


def run_module_p9(cfg: dict[str, Any], data: dict[str, dict]) -> dict[str, Any]:
    rcfg = _regional_cfg(cfg)
    bcfg = _burden_cfg(cfg)
    bands, _ = _load_bands(cfg, rcfg, data)
    if bands.empty:
        raise RuntimeError("[P9] no cytoband table available — cannot derive "
                           "arm intervals")
    arms = _arm_intervals(bands)

    src = _matrix_source(cfg, data, set(cfg["cohorts"]), rcfg)
    if src == "seg":
        print("  [P9] marker matrix (D.cap1.5.mat) missing for some cohort "
              "— focal events from tumors.seg (coarser resolution)")

    # ---- per-cohort per-sample scores ----
    scores = _cohort_scores(cfg, data, arms, bcfg, src, rcfg)
    source_by_cohort = {cid: src for cid in cfg["cohorts"]}

    # ---- rank normalization + tests ----
    meta_df = _load_metadata(cfg)
    scores, tests = _rank_and_test(cfg, scores, meta_df)

    # ---- cross-checks: S vs wGII vs FGA ----
    cross_rows = []
    for cmp_ in cfg["comparisons"]:
        ca, cb = cmp_["a"], cmp_["b"]
        comp = f"{ca}_vs_{cb}"
        pool = (scores[scores["comparison"] == comp]
                if "comparison" in scores.columns
                else scores[scores["cohort"].isin([ca, cb])])
        for x, y in (("wGII", "S_total"), ("wGII", "S_focal_star"),
                     ("fga", "S_total")):
            sub = pool[[x, y]].dropna()
            if len(sub) >= 5:
                rho = float(stats.spearmanr(sub[x], sub[y]).statistic)
            else:
                rho = np.nan
            cross_rows.append(dict(comparison=f"{ca}_vs_{cb}", x=x, y=y,
                                   n=len(sub), spearman_rho=round(rho, 3)))
    cross = pd.DataFrame(cross_rows)

    # ---- summary ----
    sum_rows = []
    for comp in dict.fromkeys(tests.get("comparison", [])):
        t = tests[tests["comparison"] == comp]
        if t.empty:
            continue
        row = dict(comparison=comp,
                   focal_source="d_matrix" if src == "d_matrix" else "seg",
                   covariates=t["covariates"].iloc[0])
        for metric, _ in METRICS:
            r = t[t["metric"] == metric]
            if len(r):
                row[f"{metric}_beta"] = r["beta"].iloc[0]
                row[f"{metric}_p"] = r["p_ols"].iloc[0]
        sum_rows.append(row)
    summary = pd.DataFrame(sum_rows)

    scores.to_csv(out_table(cfg, "P9_burden_scores.tsv"), sep="\t",
                  index=False)
    tests.to_csv(out_table(cfg, "P9_burden_tests.tsv"), sep="\t", index=False)
    cross.to_csv(out_table(cfg, "P9_crosscheck.tsv"), sep="\t", index=False)
    summary.to_csv(out_table(cfg, "P9_summary.tsv"), sep="\t", index=False)

    if not tests.empty:
        st = tests[tests["metric"] == "S_total"]
        for r in st.itertuples():
            print(f"  [P9] {r.comparison}: S beta={r.beta:+.3f} "
                  f"(95% CI {r.ci_lo:+.3f}..{r.ci_hi:+.3f}, "
                  f"p={r.p_ols:.2e}), higher in {r.cohort_higher}")

    res_extra: dict[str, Any] = {}
    if bcfg.get("sensitivity"):
        sens, verdict = run_sensitivity(cfg, data, arms, src, rcfg, bcfg,
                                        meta_df)
        sens.to_csv(out_table(cfg, "P9_sensitivity.tsv"), sep="\t",
                    index=False)
        verdict.to_csv(out_table(cfg, "P9_sensitivity_summary.tsv"),
                       sep="\t", index=False)
        for r in verdict.itertuples():
            print(f"  [P9-sens] {r.metric}: beta sign consistent in all "
                  f"{r.n_combos} combos = {r.sign_consistent} "
                  f"(p {r.p_min}-{r.p_max})")
        res_extra = dict(sensitivity=sens, sensitivity_verdict=verdict)

    return dict(scores=scores, tests=tests, crosscheck=cross, summary=summary,
                source=src, bcfg=bcfg, **res_extra)


def _load_bands(cfg: dict[str, Any], rcfg: dict, data: dict[str, dict]):
    from .pairwise_regional import _load_band_genes
    return _load_band_genes(cfg, rcfg, data)
