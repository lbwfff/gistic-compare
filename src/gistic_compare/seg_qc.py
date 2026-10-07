"""Module E: segment-level QC and fraction-genome-altered (FGA)."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from .config import HG19_CHROM_SIZES, out_table

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
GENOME_LEN = sum(HG19_CHROM_SIZES[c] for c in AUTOSOMES)


def compute_fga(seg: pd.DataFrame, log2_threshold: float = 0.2) -> pd.DataFrame:
    """Per-sample FGA on autosomes; |seg_mean| >= threshold counts as altered."""
    s = seg[seg["chrom"].isin(AUTOSOMES)].copy()
    s = s.dropna(subset=["seg_mean", "start", "end"])
    s["length"] = (s["end"] - s["start"]).clip(lower=0)
    s["is_amp"] = s["seg_mean"] >= log2_threshold
    s["is_del"] = s["seg_mean"] <= -log2_threshold
    s["is_alt"] = s["is_amp"] | s["is_del"]

    rows = []
    for sample, g in s.groupby("sample"):
        # overlapping segments possible? assume non-overlapping per sample per chrom
        total_bp = g["length"].sum()
        amp_bp = g.loc[g["is_amp"], "length"].sum()
        del_bp = g.loc[g["is_del"], "length"].sum()
        alt_bp = g.loc[g["is_alt"], "length"].sum()
        rows.append(dict(
            sample=sample,
            n_segments=len(g),
            total_seg_bp=int(total_bp),
            fga=alt_bp / GENOME_LEN,
            fga_amp=amp_bp / GENOME_LEN,
            fga_del=del_bp / GENOME_LEN,
            median_abs_log2=float(g["seg_mean"].abs().median()),
            max_abs_log2=float(g["seg_mean"].abs().max()),
            n_altered_segments=int(g["is_alt"].sum()),
            coverage_frac=total_bp / GENOME_LEN,   # should be ~1 if segs tile genome
        ))
    return pd.DataFrame(rows)


def _safe_kruskal(*groups):
    """Kruskal-Wallis that NaNs out degenerate input.

    All-identical groups raise ``ValueError`` in scipy < 1.16 (wheels
    available for Python 3.10); scipy >= 1.16 returns NaN instead.  Return
    NaN either way so a constant-input cohort (tiny/empty synthetic data,
    a single sample per group, ...) skips the group test instead of failing
    the module.
    """
    try:
        return stats.kruskal(*groups)
    except ValueError:
        return np.nan, np.nan


def run_module_e(cfg: dict[str, Any], data: dict[str, dict]) -> dict[str, Any]:
    ana = cfg["analysis"]
    thr = ana["fga"]["log2_threshold"]
    thr_sens = ana["fga"]["log2_threshold_sens"]

    qc_rows = []
    cohort_summary = []
    for cid, d in data.items():
        seg = d["seg"]
        fga_df = compute_fga(seg, thr)
        fga_df["cohort"] = cid
        fga_sens = compute_fga(seg, thr_sens)[["sample", "fga", "fga_amp", "fga_del"]]
        fga_sens = fga_sens.rename(columns={
            "fga": "fga_sens", "fga_amp": "fga_amp_sens",
            "fga_del": "fga_del_sens"})
        fga_df = fga_df.merge(fga_sens, on="sample", how="left")

        # merge sample_seg_counts
        sc = d["seg_counts"].rename(columns={"segment_count": "segment_count_reported"})
        fga_df = fga_df.merge(sc[["sample", "segment_count_reported", "included"]],
                              on="sample", how="left")

        # peaks altered count per sample (from all_lesions calls)
        calls = d["calls"]
        if not calls.empty:
            alt_per_sample = (calls.assign(a=(calls["call"] > 0))
                              .groupby("sample")["a"].sum())
            fga_df["n_altered_peaks"] = fga_df["sample"].map(alt_per_sample)
        else:
            fga_df["n_altered_peaks"] = np.nan

        qc_rows.append(fga_df)

        counts = d["seg_counts"]
        cohort_summary.append(dict(
            cohort=cid,
            n_samples_seg=int(seg["sample"].nunique()),
            n_samples_reported=int(len(counts)),
            n_included=int((counts["included"] == "yes").sum()),
            n_excluded=int((counts["included"] != "yes").sum()),
            total_segments=int(len(seg)),
            median_segments_per_sample=float(counts["segment_count"].median()),
            mean_segments_per_sample=float(counts["segment_count"].mean()),
            max_segments_per_sample=float(counts["segment_count"].max()),
            median_fga=float(fga_df["fga"].median()),
            mean_fga=float(fga_df["fga"].mean()),
            median_fga_amp=float(fga_df["fga_amp"].median()),
            median_fga_del=float(fga_df["fga_del"].median()),
            fga_threshold=thr,
        ))

    qc_all = pd.concat(qc_rows, ignore_index=True)
    summary = pd.DataFrame(cohort_summary)

    # ---------- group comparisons ----------
    stat_rows = []
    cohorts = list(data.keys())
    groups = [qc_all.loc[qc_all["cohort"] == c, "fga"].dropna().to_numpy()
              for c in cohorts]
    if all(len(g) > 0 for g in groups):
        h, p_kw = _safe_kruskal(*groups)
        stat_rows.append(dict(metric="fga", test="kruskal_wallis",
                              comparison="all", statistic=float(h),
                              pvalue=float(p_kw)))
    for cmp_ in cfg["comparisons"]:
        ca, cb = cmp_["a"], cmp_["b"]
        xa = qc_all.loc[qc_all["cohort"] == ca, "fga"].dropna().to_numpy()
        xb = qc_all.loc[qc_all["cohort"] == cb, "fga"].dropna().to_numpy()
        if len(xa) and len(xb):
            u, p = stats.mannwhitneyu(xa, xb, alternative="two-sided")
            stat_rows.append(dict(metric="fga", test="mannwhitneyu",
                                  comparison=f"{ca}_vs_{cb}",
                                  statistic=float(u), pvalue=float(p),
                                  n_a=len(xa), n_b=len(xb),
                                  median_a=float(np.median(xa)),
                                  median_b=float(np.median(xb)),
                                  kind=cmp_.get("kind", "unpaired")))
        # segment counts too
        xa2 = qc_all.loc[qc_all["cohort"] == ca, "n_segments"].dropna().to_numpy()
        xb2 = qc_all.loc[qc_all["cohort"] == cb, "n_segments"].dropna().to_numpy()
        if len(xa2) and len(xb2):
            u, p = stats.mannwhitneyu(xa2, xb2, alternative="two-sided")
            stat_rows.append(dict(metric="n_segments", test="mannwhitneyu",
                                  comparison=f"{ca}_vs_{cb}",
                                  statistic=float(u), pvalue=float(p),
                                  n_a=len(xa2), n_b=len(xb2),
                                  median_a=float(np.median(xa2)),
                                  median_b=float(np.median(xb2)),
                                  kind=cmp_.get("kind", "unpaired")))
    stats_df = pd.DataFrame(stat_rows)

    # TCGA seg samples vs reported difference
    extra_rows = []
    for cid, d in data.items():
        seg_samples = set(d["seg"]["sample"].unique())
        rep_samples = set(d["seg_counts"]["sample"].unique())
        only_seg = sorted(seg_samples - rep_samples)
        only_rep = sorted(rep_samples - seg_samples)
        for s in only_seg:
            extra_rows.append(dict(cohort=cid, sample=s, status="in_seg_not_reported"))
        for s in only_rep:
            extra_rows.append(dict(cohort=cid, sample=s, status="in_reported_not_seg"))
    sample_diff = pd.DataFrame(extra_rows,
                               columns=["cohort", "sample", "status"])

    # correlation FGA vs altered peaks
    corr_rows = []
    for cid in cohorts:
        sub = qc_all[qc_all["cohort"] == cid].dropna(subset=["fga", "n_altered_peaks"])
        if len(sub) >= 3:
            r, p = stats.spearmanr(sub["fga"], sub["n_altered_peaks"])
            corr_rows.append(dict(cohort=cid, n=len(sub),
                                  spearman_r=float(r), pvalue=float(p),
                                  metric="fga_vs_n_altered_peaks"))
    corr_df = pd.DataFrame(corr_rows)

    for name, df in (("qc_per_sample", qc_all), ("cohort_summary", summary),
                     ("group_stats", stats_df), ("sample_set_diff", sample_diff),
                     ("fga_peak_corr", corr_df)):
        df.to_csv(out_table(cfg, f"E_{name}.tsv"), sep="\t", index=False)

    return dict(qc=qc_all, summary=summary, stats=stats_df,
                sample_diff=sample_diff, fga_corr=corr_df)
