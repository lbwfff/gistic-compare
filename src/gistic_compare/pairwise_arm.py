"""Module P6 (pairwise): chromosome-arm (broad) landscape comparison.

P1-P3 compare GISTIC *focal* peaks and P4/P5 work at 1 Mb resolution; both
scales miss the third scale GISTIC reports: **broad, arm-level events**,
which carry the aneuploidy background of a cohort.  P6 compares the two
cohorts on that scale.

Broad calls per sample per arm come from GISTIC's own
``broad_values_by_arm.txt`` when every cohort provides it, otherwise they
are **derived from tumors.seg**: an arm is called amplified (deleted) for a
sample when at least ``min_frac`` of its length is covered by segments with
``seg_mean >= amp_thr`` (``<= del_thr``).  The source used is recorded in
every output table (``gistic`` vs ``derived_seg``) — the two definitions are
related but not identical and must not be mixed within one comparison.

Statistics:
  * per arm x direction: alteration frequency per cohort, delta frequency,
    two-sided Fisher exact + BH-FDR across arms (per direction)
  * per sample: broad-amp count, broad-del count and the fraction of arms
    altered (aneuploidy score); Mann-Whitney between cohorts with the
    rank-biserial effect size
  * arm-level q-values (amp_qv/del_qv.conf_XX.txt) are attached for display
    only when present in every cohort — like GS, the difference of two
    cohort-specific q-values has no null and is not tested.

Outputs: P6_arm_diff.tsv, P6_aneuploidy.tsv, P6_aneuploidy_tests.tsv,
         P6_summary.tsv
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from .config import CHROMS, out_table
from .pairwise_regional import (
    _fdr,
    _fisher_p,
    _load_band_genes,
    _regional_cfg,
)


def _broad_cfg(cfg: dict[str, Any], rcfg: dict) -> dict:
    b = (cfg.get("analysis") or {}).get("broad") or {}
    return dict(
        amp_thr=float(b.get("amp_thr", 0.30)),
        del_thr=float(b.get("del_thr", -0.30)),
        min_frac=float(b.get("min_frac", 0.50)),
    )


# ---------------------------------------------------------------------------
# arm intervals from the cytoband table
# ---------------------------------------------------------------------------

def _arm_label(band: str) -> str | None:
    """'1p36.33' -> '1p'; None when the band has no p/q prefix."""
    s = str(band).strip()
    for i, ch in enumerate(s):
        if ch in ("p", "q") and i > 0 and s[:i].isdigit():
            return f"{int(s[:i])}{ch}"
    return None


def _arm_intervals(bands: pd.DataFrame) -> pd.DataFrame:
    """Collapse cytobands into autosomal arm intervals (chrom, start, end)."""
    rows = []
    for chrom, g in bands.groupby("chrom"):
        if chrom not in CHROMS:
            continue
        g = g.assign(_arm=g["band"].map(_arm_label)).dropna(subset=["_arm"])
        for arm, gg in g.groupby("_arm"):
            rows.append(dict(arm=str(arm), chrom=chrom,
                             start=int(gg["start"].min()),
                             end=int(gg["end"].max())))
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["_num"] = out["arm"].str[:-1].astype(int)
    out["_s"] = out["arm"].str[-1]
    return (out.sort_values(["_num", "_s"])
               .drop(columns=["_num", "_s"])
               .reset_index(drop=True))


# ---------------------------------------------------------------------------
# broad calls derived from tumors.seg
# ---------------------------------------------------------------------------

def _broad_from_seg(seg: pd.DataFrame, arms: pd.DataFrame, cohort: str,
                    amp_thr: float, del_thr: float, min_frac: float
                    ) -> pd.DataFrame:
    """Derive per-sample arm calls: an arm counts as broadly altered when
    >= min_frac of its length is covered by qualifying segments."""
    samples = sorted(seg["sample"].unique())
    sidx = {s: i for i, s in enumerate(samples)}
    n_cov = np.zeros((len(arms), len(samples)))
    n_amp = np.zeros((len(arms), len(samples)))
    n_del = np.zeros((len(arms), len(samples)))
    arm_chrom = arms["chrom"].to_numpy()
    arm_start = arms["start"].to_numpy()
    arm_end = arms["end"].to_numpy()
    for c in CHROMS:
        sel_arms = np.flatnonzero(arm_chrom == c)
        if len(sel_arms) == 0:
            continue
        sub = seg[seg["chrom"] == c]
        if sub.empty:
            continue
        a_s = arm_start[sel_arms].astype(float)
        a_e = arm_end[sel_arms].astype(float)
        a_len = np.maximum(a_e - a_s, 1.0)
        for r in sub.itertuples():
            lo = np.maximum(a_s, r.start)
            hi = np.minimum(a_e, r.end)
            ov = np.clip(hi - lo, 0.0, None)
            if not ov.any():
                continue
            j = sidx[r.sample]
            n_cov[sel_arms, j] += ov
            if r.seg_mean >= amp_thr:
                n_amp[sel_arms, j] += ov
            elif r.seg_mean <= del_thr:
                n_del[sel_arms, j] += ov
    with np.errstate(invalid="ignore", divide="ignore"):
        frac_amp = np.where(n_cov > 0, n_amp / np.maximum(n_cov, 1.0), np.nan)
        frac_del = np.where(n_cov > 0, n_del / np.maximum(n_cov, 1.0), np.nan)
    rows = []
    for ai, arm in enumerate(arms["arm"]):
        for si, s in enumerate(samples):
            if np.isfinite(frac_amp[ai, si]) and frac_amp[ai, si] >= min_frac:
                rows.append(dict(cohort=cohort, sample=s, arm=arm, type="Amp"))
            if np.isfinite(frac_del[ai, si]) and frac_del[ai, si] >= min_frac:
                rows.append(dict(cohort=cohort, sample=s, arm=arm, type="Del"))
    return pd.DataFrame(rows, columns=["cohort", "sample", "arm", "type"])


def _resolve_broad_calls(cfg: dict[str, Any], data: dict[str, dict],
                         rcfg: dict, arms: pd.DataFrame,
                         bcfg: dict) -> tuple[dict[str, pd.DataFrame], str]:
    """Per-cohort long arm-call frames + the source label used for all."""
    native = {}
    for cid, d in data.items():
        if d.get("broad") is not None and len(d["broad"]):
            native[cid] = d["broad"]
    if len(native) == len(data):
        return native, "gistic"
    if native:
        print("  [P6] broad_values_by_arm.txt present in only "
              f"{len(native)}/{len(data)} cohorts — deriving arm calls "
              "from tumors.seg for all cohorts (consistent definition)")
    calls = {}
    for cid, d in data.items():
        calls[cid] = _broad_from_seg(
            d["seg"], arms, cid, bcfg["amp_thr"], bcfg["del_thr"],
            bcfg["min_frac"])
    return calls, "derived_seg"


def _rank_biserial(xa: np.ndarray, xb: np.ndarray) -> float:
    """rank-biserial r = 2*U1/(n1*n2) - 1 (positive when A > B)."""
    try:
        u1 = stats.mannwhitneyu(xa, xb, alternative="two-sided")[0]
    except ValueError:
        return np.nan
    return 2.0 * u1 / (len(xa) * len(xb)) - 1.0


def _aneuploidy_table(broad: dict[str, pd.DataFrame], arms: pd.DataFrame,
                      data: dict[str, dict]) -> pd.DataFrame:
    """Per-sample broad-event counts / aneuploidy score.

    Every cohort sample is listed (zeros for samples without broad events)
    so the between-cohort distributions are complete."""
    rows = []
    n_arms = len(arms)
    for cid, b in broad.items():
        counts: dict[str, tuple[int, int]] = {}
        if not b.empty:
            for s, g in b.groupby("sample"):
                counts[str(s)] = (int((g["type"] == "Amp").sum()),
                                  int((g["type"] == "Del").sum()))
        for s in sorted(data[cid]["seg"]["sample"].unique()):
            n_amp, n_del = counts.get(str(s), (0, 0))
            rows.append(dict(cohort=cid, sample=str(s),
                             n_broad_amp=n_amp, n_broad_del=n_del,
                             broad_frac=(n_amp + n_del) / n_arms
                             if n_arms else np.nan))
    return pd.DataFrame(rows, columns=[
        "cohort", "sample", "n_broad_amp", "n_broad_del", "broad_frac"])


def run_module_p6(cfg: dict[str, Any], data: dict[str, dict]) -> dict[str, Any]:
    rcfg = _regional_cfg(cfg)
    bcfg = _broad_cfg(cfg, rcfg)
    bands, _ = _load_band_genes(cfg, rcfg, data)
    if bands.empty:
        raise RuntimeError("[P6] no cytoband table available — cannot derive "
                           "arm intervals")
    arms = _arm_intervals(bands)
    if arms.empty:
        raise RuntimeError("[P6] cytoband table yielded no chromosome arms")

    broad, source = _resolve_broad_calls(cfg, data, rcfg, arms, bcfg)
    arm_qv: dict[str, dict[str, pd.DataFrame]] = {
        cid: (d.get("arm_qv") or {}) for cid, d in data.items()
    }
    qv_complete = all(
        {"Amp", "Del"} <= set(arm_qv[cid]) and len(arm_qv[cid]["Amp"])
        and len(arm_qv[cid]["Del"])
        for cid in data)

    arm_meta = arms.set_index("arm")[["chrom", "start", "end"]].to_dict("index")
    fdr_thr = rcfg["fdr"]

    diff_frames, aneu_frames, test_rows, summary_rows = [], [], [], []
    for cmp_ in cfg["comparisons"]:
        ca, cb = cmp_["a"], cmp_["b"]
        comp = f"{ca}_vs_{cb}"
        if cmp_.get("kind") == "paired":
            print(f"  [P6] skipped paired comparison {comp}")
            continue
        la = cfg.get("cohorts", {}).get(ca, {}).get("label", ca)
        lb = cfg.get("cohorts", {}).get(cb, {}).get("label", cb)

        # ---- per-sample aneuploidy metrics ----
        aneu = _aneuploidy_table({ca: broad[ca], cb: broad[cb]}, arms, data)
        aneu = aneu.sort_values(["cohort", "sample"]).reset_index(drop=True)
        aneu.insert(0, "comparison", comp)
        aneu.insert(1, "kind", cmp_.get("kind", "unpaired"))
        aneu_frames.append(aneu)
        ga = aneu[aneu["cohort"] == ca]
        gb = aneu[aneu["cohort"] == cb]
        tests = dict(comparison=comp)
        for metric in ("n_broad_amp", "n_broad_del", "broad_frac"):
            xa = ga[metric].to_numpy(float)
            xb = gb[metric].to_numpy(float)
            if len(xa) >= 3 and len(xb) >= 3 and (np.unique(np.concatenate(
                    [xa, xb])).size > 1):
                p = float(stats.mannwhitneyu(xa, xb,
                                             alternative="two-sided")[1])
                r = _rank_biserial(xa, xb)
            else:
                p, r = np.nan, np.nan
            tests.update({
                f"{metric}_median_a": float(np.median(xa)) if len(xa) else np.nan,
                f"{metric}_median_b": float(np.median(xb)) if len(xb) else np.nan,
                f"{metric}_p_mw": p, f"{metric}_r_rb": r,
            })
        test_rows.append(tests)

        # ---- arm-level frequency differences ----
        rows = []
        ba = broad[ca].set_index(["sample", "arm", "type"]) if len(broad[ca]) \
            else pd.DataFrame()
        bb = broad[cb].set_index(["sample", "arm", "type"]) if len(broad[cb]) \
            else pd.DataFrame()
        na, nb = (data[ca]["seg"]["sample"].nunique(),
                  data[cb]["seg"]["sample"].nunique())
        for _, arm_row in arms.iterrows():
            arm = str(arm_row["arm"])
            rec = dict(comparison=comp, kind=cmp_.get("kind", "unpaired"),
                       source=source, arm=arm,
                       chrom=arm_row["chrom"], start=int(arm_row["start"]),
                       end=int(arm_row["end"]))
            for dirn in ("Amp", "Del"):
                ka = int(((ba.index.get_level_values("arm") == arm)
                          & (ba.index.get_level_values("type") == dirn)).sum()
                         ) if len(ba) else 0
                kb = int(((bb.index.get_level_values("arm") == arm)
                          & (bb.index.get_level_values("type") == dirn)).sum()
                         ) if len(bb) else 0
                d = dirn.lower()
                rec.update({
                    f"freq_{d}_a": ka / na if na else np.nan,
                    f"freq_{d}_b": kb / nb if nb else np.nan,
                    f"n_{d}_a": ka, f"n_{d}_b": kb,
                    f"p_{d}": _fisher_p(ka, na, kb, nb),
                })
                if qv_complete:
                    qa = arm_qv[ca][dirn]
                    qb = arm_qv[cb][dirn]
                    rec[f"q_{d}_a"] = float(
                        qa.loc[qa["arm"] == arm, "q"].min()) if len(qa.loc[
                            qa["arm"] == arm]) else np.nan
                    rec[f"q_{d}_b"] = float(
                        qb.loc[qb["arm"] == arm, "q"].min()) if len(qb.loc[
                            qb["arm"] == arm]) else np.nan
            rows.append(rec)
        diff = pd.DataFrame(rows)
        diff["delta_freq_amp"] = diff["freq_amp_a"] - diff["freq_amp_b"]
        diff["delta_freq_del"] = diff["freq_del_a"] - diff["freq_del_b"]
        for d in ("amp", "del"):
            diff[f"fdr_{d}"] = _fdr(diff[f"p_{d}"].to_numpy(float))
            diff[f"sig_{d}"] = diff[f"fdr_{d}"] < fdr_thr
        diff_frames.append(diff)

        summary_rows.append(dict(
            comparison=comp, kind=cmp_.get("kind", "unpaired"), source=source,
            n_arms=len(diff), fdr_thr=fdr_thr,
            n_sig_amp=int(diff["sig_amp"].sum()),
            n_sig_del=int(diff["sig_del"].sum()),
            max_abs_delta_amp=round(float(np.nanmax(np.abs(
                diff["delta_freq_amp"]))), 4),
            max_abs_delta_del=round(float(np.nanmax(np.abs(
                diff["delta_freq_del"]))), 4),
            n_broad_amp_median_a=tests["n_broad_amp_median_a"],
            n_broad_amp_median_b=tests["n_broad_amp_median_b"],
            n_broad_amp_p_mw=tests["n_broad_amp_p_mw"],
            n_broad_del_median_a=tests["n_broad_del_median_a"],
            n_broad_del_median_b=tests["n_broad_del_median_b"],
            n_broad_del_p_mw=tests["n_broad_del_p_mw"],
            broad_frac_median_a=tests["broad_frac_median_a"],
            broad_frac_median_b=tests["broad_frac_median_b"],
            broad_frac_p_mw=tests["broad_frac_p_mw"],
        ))
        _ = (la, lb)

    diff_all = (pd.concat(diff_frames, ignore_index=True) if diff_frames
                else pd.DataFrame())
    aneu_all = (pd.concat(aneu_frames, ignore_index=True) if aneu_frames
                else pd.DataFrame())
    tests_df = pd.DataFrame(test_rows)
    summary = pd.DataFrame(summary_rows)

    diff_all.to_csv(out_table(cfg, "P6_arm_diff.tsv"), sep="\t", index=False)
    aneu_all.to_csv(out_table(cfg, "P6_aneuploidy.tsv"), sep="\t", index=False)
    tests_df.to_csv(out_table(cfg, "P6_aneuploidy_tests.tsv"), sep="\t",
                    index=False)
    summary.to_csv(out_table(cfg, "P6_summary.tsv"), sep="\t", index=False)

    if source == "gistic":
        print("  [P6] arm calls from GISTIC broad_values_by_arm.txt")
    else:
        print(f"  [P6] arm calls derived from tumors.seg "
              f"(amp>={bcfg['amp_thr']}, del<={bcfg['del_thr']}, "
              f"min_frac>={bcfg['min_frac']})")
    if not summary.empty:
        print(f"  [P6] sig arms (FDR<{fdr_thr}): "
              + "; ".join(f"{r.comparison}: amp={r.n_sig_amp} "
                          f"del={r.n_sig_del}"
                          for r in summary.itertuples()))
    return dict(arm_diff=diff_all, aneuploidy=aneu_all,
                aneuploidy_tests=tests_df, summary=summary, arms=arms,
                source=source, bcfg=bcfg)
