"""Module F: novelty of experiment peaks relative to control.

For each experiment vs control (star topology), measure the fraction of the
experiment's peak "signal" that has no moderate match in the control.

Weights:
  count      — n_unique / n_exp
  freq       — sum(freq_altered unique) / sum(freq_altered all)
  neglogq    — sum(-log10 q unique) / sum(-log10 q all)
  gscore     — sum(G-score unique) / sum(G-score all), G-score aggregated
               from scores.gistic intervals overlapping each peak
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from . import intervals as iv
from .config import out_table

WEIGHTS = ("count", "freq", "neglogq", "gscore")


def _aggregate_gscore(scores: pd.DataFrame, peaks: pd.DataFrame) -> pd.Series:
    """Length-weighted mean G-score of scores.gistic intervals per peak.

    Returns Series **indexed by peak_id** (peak_id -> gscore).  Downstream
    frames (q-filtered, one-to-one-match remainders) carry different row
    orders than the full peaks table, so scores must be attached by id —
    never by position.
    """
    rows = []
    for _, pk in peaks.iterrows():
        m = (
            (scores["type"] == pk["type"])
            & (scores["chrom"] == pk["chrom"])
            & (scores["start"] < pk["end"])
            & (scores["end"] > pk["start"])
        )
        sub = scores.loc[m]
        if sub.empty:
            rows.append(np.nan)
            continue
        ov = (
            np.minimum(sub["end"].to_numpy(), pk["end"])
            - np.maximum(sub["start"].to_numpy(), pk["start"])
        ).clip(min=0).astype(float)
        g = sub["g_score"].to_numpy(float)
        good = (ov > 0) & ~np.isnan(g)
        if not good.any() or ov[good].sum() <= 0:
            rows.append(float(np.nanmean(g)) if np.isfinite(g).any() else np.nan)
            continue
        rows.append(float(np.average(g[good], weights=ov[good])))
    return pd.Series(rows, index=peaks["peak_id"].to_numpy())


def _attach_gscore(peaks: pd.DataFrame, gscore: pd.Series) -> np.ndarray:
    """Attach per-peak G-scores to a peaks frame **by peak_id**.

    The frame may be a q-filtered subset or a match remainder with a reset
    RangeIndex; positional lookup would hand it another peak's G-score.
    """
    return gscore.reindex(peaks["peak_id"]).to_numpy()


def _frac(num: float, den: float) -> float:
    if not np.isfinite(den) or den <= 0:
        return np.nan
    return float(num / den)


def run_module_f(cfg: dict[str, Any], data: dict[str, dict]) -> dict[str, Any]:
    ana = cfg["analysis"]
    q_thresholds = ana["q_thresholds"]
    ov = ana["overlap"]
    control = cfg.get("control")
    if not control or control not in data:
        raise ValueError(
            "module F requires comparison_mode=control_vs_all with a valid "
            f"'control' (got {control!r})")

    # precompute peak-level G-score aggregates per cohort
    gscores: dict[str, pd.Series] = {}
    for cid, d in data.items():
        gscores[cid] = _aggregate_gscore(d["scores"], d["peaks"])

    # experiments = all non-control cohorts (star topology)
    experiments = [c for c in data if c != control]
    if not experiments:
        raise ValueError("module F: no experiment cohorts besides control")

    summary_rows: list[dict[str, Any]] = []
    peak_rows: list[dict[str, Any]] = []

    for exp in experiments:
        for q in q_thresholds:
            for ptype in ("Amp", "Del"):
                fe = iv.filter_peaks(data[exp]["peaks"], q, ptype).copy()
                fc = iv.filter_peaks(data[control]["peaks"], q, ptype).copy()
                if fe.empty:
                    continue
                # moderate = one-sided reciprocal (same as module A primary)
                res = iv.match_peaks(fe, fc, min_bp=ov["min_bp"],
                                     reciprocal=ov["reciprocal_shorter"],
                                     mode="either")
                n_exp = res.n_a  # a = experiment
                n_unique = len(res.only_a)
                n_matched = len(res.pairs)

                # attach per-peak weights for unique peaks
                uniq = res.only_a.copy()
                uniq["gscore"] = _attach_gscore(uniq, gscores[exp])
                uniq["neglogq"] = -np.log10(np.clip(uniq["q"].to_numpy(float),
                                                    1e-300, None))
                allp = fe.copy()
                allp["gscore"] = _attach_gscore(allp, gscores[exp])
                allp["neglogq"] = -np.log10(np.clip(allp["q"].to_numpy(float),
                                                    1e-300, None))
                allp["_unique"] = allp["peak_id"].isin(uniq["peak_id"])

                def _sum(df: pd.DataFrame, col: str, mask_col: str | None = None) -> float:
                    sub = df if mask_col is None else df[df[mask_col]]
                    v = sub[col].to_numpy(float)
                    v = v[np.isfinite(v)]
                    return float(v.sum())

                tot = {
                    "count": float(n_exp),
                    "freq": _sum(allp, "freq_altered"),
                    "neglogq": _sum(allp, "neglogq"),
                    "gscore": _sum(allp, "gscore"),
                }
                uni = {
                    "count": float(n_unique),
                    "freq": _sum(uniq, "freq_altered"),
                    "neglogq": _sum(uniq, "neglogq"),
                    "gscore": _sum(uniq, "gscore"),
                }
                fracs = {w: _frac(uni[w], tot[w]) for w in WEIGHTS}

                summary_rows.append(dict(
                    experiment=exp, control=control,
                    comparison=f"{control}_vs_{exp}",
                    q_thresh=q, ptype=ptype,
                    n_exp=n_exp, n_unique=n_unique, n_matched=n_matched,
                    overlap_frac=n_matched / n_exp if n_exp else np.nan,
                    **{f"frac_{k}": v for k, v in fracs.items()},
                    **{f"total_{k}": tot[k] for k in WEIGHTS},
                    **{f"unique_{k}": uni[k] for k in WEIGHTS},
                ))

                for _, r in uniq.iterrows():
                    peak_rows.append(dict(
                        experiment=exp, control=control,
                        comparison=f"{control}_vs_{exp}",
                        q_thresh=q, ptype=ptype,
                        peak_id=r["peak_id"], cytoband=r["cytoband"],
                        chrom=r["chrom"], start=int(r["start"]),
                        end=int(r["end"]),
                        q=float(r["q"]), freq_altered=float(r["freq_altered"]),
                        neglogq=float(r["neglogq"]),
                        gscore=float(r["gscore"]) if np.isfinite(r["gscore"]) else np.nan,
                    ))

    summary = pd.DataFrame(summary_rows)
    peaks_uniq = pd.DataFrame(peak_rows)
    summary.to_csv(out_table(cfg, "F_novelty_summary.tsv"), sep="\t", index=False)
    peaks_uniq.to_csv(out_table(cfg, "F_novelty_peaks.tsv"), sep="\t", index=False)
    return dict(summary=summary, unique_peaks=peaks_uniq)
