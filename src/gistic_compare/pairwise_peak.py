"""Module P1 (pairwise workflow): peak-level set differences between experiments.

For each experiment pair: moderate-matched pairs with delta q/freq, one-sided
unique peaks, Jaccard across match tiers, multi-experiment membership for
UpSet, and delta-frequency rank tests. Score-track correlations (B machinery,
Pearson + Spearman) are reused as P1_track_correlations.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from . import intervals as iv
from .config import out_table
from .pairwise_model import build_union_keys
from .score_track import make_grid, resample_track


def _fdr_bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    out = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    from statsmodels.stats.multitest import multipletests
    _, fdr, _, _ = multipletests(p[ok], method="fdr_bh")
    out[ok] = fdr
    return out


def _tier_kwargs(ov: dict[str, Any]) -> dict[str, tuple]:
    """Match tier presets from config overlap block."""
    return {
        "moderate": dict(min_bp=ov["min_bp"],
                         reciprocal=ov["reciprocal_shorter"], mode="either"),
        "strict": dict(min_bp=ov["min_bp"],
                       reciprocal=ov["reciprocal_both"], mode="both"),
        "loose": dict(min_bp=ov["min_bp"], reciprocal=0.0, mode="any"),
    }


def build_onco_data(cfg: dict[str, Any], data: dict[str, dict],
                    q_main: float) -> list[dict[str, Any]]:
    """Union-peak × sample call matrices for OncoPrint figures.

    One entry per comparison.  Rows are union peak keys (moderate-match
    chains from ``build_union_keys``), columns are the samples of the two
    cohorts side by side; values are the 0/1/2 amplitude-threshold calls
    (max over the peaks chained into the key).  A cell is NaN when the
    cohort has **no member peak** for that key — i.e. GISTIC did not call a
    peak there, which is shown as grey rather than as "no alteration".
    """
    keys = build_union_keys(cfg, data, q_main)
    if keys.empty:
        return []

    out = []
    for cmp_ in cfg["comparisons"]:
        ca, cb = cmp_["a"], cmp_["b"]
        comp = f"{ca}_vs_{cb}"
        ksub = keys[keys["cohort"].isin([ca, cb])].copy()
        if ksub.empty:
            continue

        # per-cohort per-sample max call over the peaks chained into a key
        call_frames = []
        for cid in (ca, cb):
            kc = ksub[ksub["cohort"] == cid]
            if kc.empty:
                continue
            calls = data[cid]["calls"]
            sub = calls[calls["peak_id"].isin(set(kc["peak_id"]))].copy()
            if sub.empty:
                continue
            pid2key = dict(zip(kc["peak_id"], kc["peak_key"]))
            sub["peak_key"] = sub["peak_id"].map(pid2key)
            sub = sub[sub["peak_key"].notna()]
            sub["call"] = pd.to_numeric(sub["call"], errors="coerce")
            agg = (sub.groupby(["peak_key", "sample"], as_index=False)["call"]
                      .max())
            agg["cohort"] = cid
            call_frames.append(agg)
        if not call_frames:
            continue
        long = pd.concat(call_frames, ignore_index=True)

        samples = (sorted(data[ca]["calls"]["sample"].unique())
                   + sorted(data[cb]["calls"]["sample"].unique()))
        matrix = long.pivot_table(index="peak_key", columns="sample",
                                  values="call", aggfunc="max")
        matrix = matrix.reindex(columns=samples)

        # row metadata: coordinates from the first member, per-cohort q/freq
        # plus the chain's *narrow* (Peak Limits) and *wide* (Wide Peak
        # Limits) spans, so a row can be labelled at either granularity
        pk_idx = {cid: data[cid]["peaks"].set_index("peak_id")
                  for cid in (ca, cb) if len(data[cid]["peaks"])}
        meta_rows = []
        for key, g in ksub.groupby("peak_key", sort=False):
            first = g.sort_values(["chrom", "start"]).iloc[0]
            row = dict(peak_key=key, ptype=first["ptype"],
                       cytoband=first["cytoband"], chrom=first["chrom"],
                       start=int(first["start"]), end=int(first["end"]))
            nar_s: list[int] = []
            nar_e: list[int] = []
            for _, m in g.iterrows():
                pk = pk_idx.get(m["cohort"])
                if pk is None or m["peak_id"] not in pk.index:
                    continue
                ref = pk.loc[m["peak_id"]]
                try:
                    if np.isfinite(float(ref["peak_start"])):
                        nar_s.append(int(ref["peak_start"]))
                        nar_e.append(int(ref["peak_end"]))
                except (TypeError, ValueError):
                    continue
            if nar_s:                      # the chain's narrow-peak core
                row["narrow_start"] = int(min(nar_s))
                row["narrow_end"] = int(max(nar_e))
            for cid, side in ((ca, "a"), (cb, "b")):
                m = g[g["cohort"] == cid]
                row[f"q_{side}"] = (float(m["q"].min()) if len(m) else np.nan)
                row[f"freq_{side}"] = (float(m["freq_altered"].max())
                                       if len(m) else np.nan)
                row[f"present_{side}"] = bool(len(m))
            meta_rows.append(row)
        meta = pd.DataFrame(meta_rows).set_index("peak_key")

        # row order: Amp above Del, then genomic
        from .config import CHROM_RANK
        meta["_r"] = meta["chrom"].map(CHROM_RANK)
        meta = meta.sort_values(["ptype", "_r", "start"],
                                ascending=[False, True, True])
        order = meta.index
        meta = meta.drop(columns="_r")
        matrix = matrix.reindex(order)
        long_out = long[long["peak_key"].isin(set(order))]

        blocks = []
        for cid in (ca, cb):
            cols = [s for s in samples
                    if s in set(data[cid]["calls"]["sample"].unique())]
            blocks.append(dict(cohort=cid,
                               label=cfg["cohorts"][cid].get("label", cid),
                               columns=cols))
        out.append(dict(comparison=comp, kind=cmp_.get("kind", "unpaired"),
                        a=ca, b=cb, matrix=matrix, meta=meta,
                        calls=long_out, blocks=blocks))
    return out


def run_module_p1(cfg: dict[str, Any], data: dict[str, dict]) -> dict[str, Any]:
    ana = cfg["analysis"]
    q_main = ana["q_thresholds"][0]
    ov = ana["overlap"]
    tiers = _tier_kwargs(ov)

    match_rows: list[dict[str, Any]] = []
    unique_rows: list[dict[str, Any]] = []
    jaccard_rows: list[dict[str, Any]] = []
    delta_rows: list[dict[str, Any]] = []

    for cmp_ in cfg["comparisons"]:
        ca, cb = cmp_["a"], cmp_["b"]
        kind = cmp_.get("kind", "unknown")
        comp = f"{ca}_vs_{cb}"

        for tier_name, kwargs in tiers.items():
            for ptype in ("Amp", "Del"):
                fa = iv.filter_peaks(data[ca]["peaks"], q_main, ptype)
                fb = iv.filter_peaks(data[cb]["peaks"], q_main, ptype)
                res = iv.match_peaks(fa, fb, **kwargs)
                jaccard_rows.append(dict(
                    comparison=comp, kind=kind, ptype=ptype,
                    tier=tier_name, n_a=res.n_a, n_b=res.n_b,
                    n_pairs=len(res.pairs),
                    jaccard=iv.jaccard(res.n_a, res.n_b, len(res.pairs)),
                ))
                if tier_name != ov.get("default_tier", "moderate"):
                    continue

                # matched pairs with deltas
                if len(res.pairs):
                    def _ann(df: pd.DataFrame, side: str) -> pd.DataFrame:
                        out = df[["peak_id", "q", "freq_altered"]].copy()
                        if {"wide_chrom", "wide_start", "wide_end"} <= set(df.columns):
                            out[side] = (
                                df["wide_chrom"].astype(str) + ":"
                                + df["wide_start"].astype(str) + "-"
                                + df["wide_end"].astype(str))
                        else:
                            out[side] = ""
                        return out.rename(columns={
                            "peak_id": "a_peak_id" if side == "a_wide"
                            else "b_peak_id",
                            "q": "a_q" if side == "a_wide" else "b_q",
                            "freq_altered": ("a_freq" if side == "a_wide"
                                             else "b_freq"),
                        })

                    m = res.pairs.merge(_ann(fa, "a_wide"), on="a_peak_id",
                                        how="left")
                    m = m.merge(_ann(fb, "b_wide"), on="b_peak_id", how="left")
                    m["comparison"] = comp
                    m["kind"] = kind
                    m["ptype"] = ptype
                    m["delta_q"] = m["a_q"] - m["b_q"]
                    m["delta_freq"] = m["a_freq"] - m["b_freq"]
                    m["delta_neglogq"] = (
                        -np.log10(np.clip(m["a_q"].to_numpy(float), 1e-300, None))
                        + np.log10(np.clip(m["b_q"].to_numpy(float), 1e-300, None)))
                    match_rows.append(m)

                    dsub = m[["comparison", "kind", "ptype", "a_peak_id",
                              "b_peak_id", "a_cytoband", "b_cytoband",
                              "a_freq", "b_freq", "delta_freq"]].copy()
                    delta_rows.append(dsub)

                        # unique peaks
                for side, df, prefix in (("a", res.only_a, ca),
                                         ("b", res.only_b, cb)):
                    for _, r in df.iterrows():
                        if {"wide_chrom", "wide_start", "wide_end"} <= set(df.columns):
                            wp = (f"{r['wide_chrom']}:{int(r['wide_start'])}-"
                                  f"{int(r['wide_end'])}")
                        else:
                            wp = ""
                        unique_rows.append(dict(
                            comparison=comp, kind=kind, ptype=ptype,
                            side=side, cohort=prefix,
                            peak_id=r["peak_id"], cytoband=r.get("cytoband", ""),
                            chrom=r["chrom"], start=int(r["start"]),
                            end=int(r["end"]),
                            wide_peak=wp,
                            q=float(r["q"]),
                            freq_altered=float(r["freq_altered"]),
                            neglogq=float(-np.log10(max(float(r["q"]), 1e-300))),
                        ))

    # global union keys → membership matrix for UpSet (>=3 cohorts)
    memberships: dict[str, pd.DataFrame] = {"Amp": pd.DataFrame(), "Del": pd.DataFrame()}
    if len(cfg["cohorts"]) >= 3:
        keys = build_union_keys(cfg, data, q_main)
        if not keys.empty:
            memberships = {
                ptype: keys.loc[keys["ptype"] == ptype,
                                ["peak_key", "cohort", "peak_id", "cytoband"]]
                       .rename(columns={"cytoband": "cyto"})
                for ptype in ("Amp", "Del")
            }

    # ---- delta stats (primary tier): rank tests on matched delta_freq ----
    delta_rows_all = delta_rows
    stats_rows: list[dict[str, Any]] = []
    if delta_rows_all:
        ddf = pd.concat(delta_rows_all, ignore_index=True)
        for (comp, ptype), g in ddf.groupby(["comparison", "ptype"]):
            kind = g["kind"].iloc[0]
            dfreq = g["delta_freq"].to_numpy(float)
            if len(dfreq) < 3:
                p = np.nan
                test = "insufficient"
            else:
                # matched peak pairs are the paired units; H0: median delta = 0
                if np.allclose(dfreq, 0):
                    p = 1.0
                else:
                    p = float(stats.wilcoxon(dfreq, alternative="two-sided").pvalue)
                test = "wilcoxon_delta"
            stats_rows.append(dict(
                comparison=comp, ptype=ptype, kind=kind,
                n_pairs=len(dfreq), test=test,
                median_delta_freq=float(np.median(dfreq)),
                mean_delta_freq=float(np.mean(dfreq)),
                pvalue=p,
            ))
        sdf = pd.DataFrame(stats_rows)
        if not sdf.empty and sdf["pvalue"].notna().any():
            sdf["fdr"] = _fdr_bh(sdf["pvalue"].to_numpy(float))
        else:
            sdf = pd.DataFrame(stats_rows)
            sdf["fdr"] = np.nan
    else:
        sdf = pd.DataFrame(columns=[
            "comparison", "ptype", "kind", "n_pairs", "test",
            "median_delta_freq", "mean_delta_freq", "pvalue", "fdr"])
        ddf = pd.DataFrame(columns=[
            "comparison", "kind", "ptype", "a_peak_id", "b_peak_id",
            "a_cytoband", "b_cytoband", "a_freq", "b_freq", "delta_freq"])

    match_df = (pd.concat(match_rows, ignore_index=True)
                if match_rows else pd.DataFrame())
    unique_df = (pd.DataFrame(unique_rows) if unique_rows
                 else pd.DataFrame(columns=[
                     "comparison", "kind", "ptype", "side", "cohort",
                     "peak_id", "cytoband", "chrom", "start", "end",
                     "wide_peak", "q", "freq_altered", "neglogq"]))
    jaccard_df = pd.DataFrame(jaccard_rows)
    if not isinstance(memberships, dict):
        memberships = {"Amp": pd.DataFrame(), "Del": pd.DataFrame()}
    member_df = {
        ptype: (df if df is not None and not df.empty else pd.DataFrame(
            columns=["peak_key", "cohort", "peak_id", "cyto"]))
        for ptype, df in memberships.items()
    }

    # ---- score-track correlation (reuse B machinery) ----
    bin_size = ana["score_track"]["bin_size"]
    grid = make_grid(bin_size)
    metrics = ["neg_log10_q", "g_score", "frequency"]
    tracks: dict[tuple[str, str, str], np.ndarray] = {}
    for cid, d in data.items():
        for ptype in ("Amp", "Del"):
            for m in metrics:
                tracks[(cid, ptype, m)] = resample_track(
                    d["scores"], grid, m, ptype).to_numpy()
    corr_rows = []
    for cmp_ in cfg["comparisons"]:
        ca, cb = cmp_["a"], cmp_["b"]
        for ptype in ("Amp", "Del"):
            for m in metrics:
                xa = tracks[(ca, ptype, m)]
                xb = tracks[(cb, ptype, m)]
                mask = np.isfinite(xa) & np.isfinite(xb)
                n = int(mask.sum())
                if n >= 3:
                    r_p, p_p = stats.pearsonr(xa[mask], xb[mask])
                    r_s, p_s = stats.spearmanr(xa[mask], xb[mask])
                else:
                    r_p = p_p = r_s = p_s = np.nan
                corr_rows.append(dict(
                    comparison=f"{ca}_vs_{cb}",
                    kind=cmp_.get("kind", "unknown"),
                    ptype=ptype, metric=m, n_bins=n,
                    pearson_r=r_p, pearson_p=p_p,
                    spearman_r=r_s, spearman_p=p_s,
                ))
    corr_df = pd.DataFrame(corr_rows)

    # ---- delta score grid (A − B on −log10 q), 1 Mb ----
    delta_grid = grid.copy()
    for cmp_ in cfg["comparisons"]:
        ca, cb = cmp_["a"], cmp_["b"]
        for ptype in ("Amp", "Del"):
            xa = tracks[(ca, ptype, "neg_log10_q")]
            xb = tracks[(cb, ptype, "neg_log10_q")]
            delta_grid[f"{ca}_vs_{cb}_{ptype}"] = xa - xb
    # long form for plotting: chrom/start/end + one col per comp×ptype already

    # ---- union-peak OncoPrint data ----
    onco = build_onco_data(cfg, data, q_main)
    if onco:
        all_calls = pd.concat([o["calls"] for o in onco], ignore_index=True)
        all_calls.to_csv(out_table(cfg, "P1_onco_calls.tsv"), sep="\t",
                         index=False)
        print(f"  [P1] oncoprint: {len(all_calls)} call rows, "
              f"{sum(len(o['matrix']) for o in onco)} union peaks")

    jaccard_df.to_csv(out_table(cfg, "P1_jaccard.tsv"), sep="\t", index=False)
    match_df.to_csv(out_table(cfg, "P1_matched_pairs.tsv"), sep="\t", index=False)
    unique_df.to_csv(out_table(cfg, "P1_unique_peaks.tsv"), sep="\t", index=False)
    ddf.to_csv(out_table(cfg, "P1_delta_stats_input.tsv"), sep="\t", index=False)
    sdf.to_csv(out_table(cfg, "P1_delta_stats.tsv"), sep="\t", index=False)
    for ptype, df in member_df.items():
        df.to_csv(out_table(cfg, f"P1_membership_{ptype.lower()}.tsv"),
                  sep="\t", index=False)
    corr_df.to_csv(out_table(cfg, "P1_track_correlations.tsv"),
                   sep="\t", index=False)
    delta_grid.to_csv(out_table(cfg, "P1_delta_track.tsv"),
                      sep="\t", index=False)

    return dict(jaccard=jaccard_df, matched=match_df, unique=unique_df,
                delta=ddf, delta_stats=sdf, membership=member_df,
                correlations=corr_df, delta_grid=delta_grid, onco=onco)
