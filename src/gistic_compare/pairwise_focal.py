"""Module P11 (pairwise): genome-wide focal q-value landscape.

The arm- and chromosome-level figures (P6 / P9) summarise *where* broad events
happen; P11 shows GISTIC's own **focal** resolution: each cohort's
``scores.gistic`` ``-log10(q)`` track drawn over the whole genome, Amplifications
above the significance threshold and Deletions mirrored below it.

Source: ``<prefix>.scores.gistic`` (already parsed by the loader into
``data[cid]["scores"]`` with a ``neg_log10_q`` column).  Peak genes for the
annotations come from ``<prefix>.amp_genes.conf_90.txt`` /
``.del_genes.conf_90.txt`` (``data[cid]["genes"]``): one row per GISTIC peak
(top gene first), reduced here to the strongest gene per peak.

This is a **descriptive overlay**, not a between-cohort test: q-values are
computed within each cohort, so the difference of two cohorts' q-values has no
null distribution (the same caveat as P6's arm q-values and module GS).  The
only reference line is the significance threshold ``q_threshold`` (default
0.25, or the first entry of ``analysis.q_thresholds``).

Outputs: P11_focal_peaks.tsv
"""

from __future__ import annotations

import re
from typing import Any

import numpy as np
import pandas as pd

from .config import CHROMS, out_table

# GISTIC lists every gene in a wide peak; the label is more useful as a
# protein-coding driver, so prefer the first non-miRNA / non-LOC / non-linc
# name and fall back to the top gene when the peak only contains those.
_SKIP_LABEL = re.compile(r"^(hsa-(mir|let)|MIR|LOC|LINC|CTD-|AC\d)", re.I)


def _focal_cfg(cfg: dict[str, Any]) -> dict:
    analysis = cfg.get("analysis") or {}
    f = analysis.get("focal") or {}
    q = f.get("q_threshold")
    if q is None:
        qs = analysis.get("q_thresholds") or [0.25]
        q = float(qs[0])
    return dict(q_threshold=float(q))


def _parse_wide_peak(text) -> tuple[str, float, float]:
    """'chr11:69334308-69344762' -> (chrom, start, end); NaN on odd input."""
    s = str(text).strip()
    if ":" not in s or "-" not in s:
        return "", np.nan, np.nan
    chrom, rng = s.split(":", 1)
    a, _, b = rng.partition("-")
    try:
        return chrom, float(a), float(b)
    except ValueError:
        return chrom, np.nan, np.nan


def _clean_gene(g: str) -> str:
    """GISTIC brackets its known cancer genes ('[CCND1]'); strip them."""
    return str(g).strip().strip("[]").strip()


def _pick_gene(genes: list[str]) -> str:
    """First protein-coding-looking gene name, else the top gene."""
    for g in genes:
        cg = _clean_gene(g)
        if cg and not _SKIP_LABEL.match(cg):
            return cg
    return _clean_gene(genes[0]) if genes else ""


def _peak_records(genes: pd.DataFrame | None, cohort: str) -> pd.DataFrame:
    """One row per GISTIC peak: cohort, type, cytoband, locus, q, label gene.

    The genes table lists every gene of every wide peak in file order, the
    top gene first; the label gene is the first protein-coding-looking name
    (see ``_SKIP_LABEL``), falling back to the top gene.
    """
    cols = ["cohort", "type", "cytoband", "chrom", "start", "end", "q",
            "neg_log10_q", "gene"]
    if genes is None or len(genes) == 0 or "col_index" not in genes.columns:
        return pd.DataFrame(columns=cols)
    rows = []
    grouped = genes.groupby(["type", "col_index"], sort=False)
    for (ptype, _ci), g in grouped:
        first = g.iloc[0]
        chrom, start, end = _parse_wide_peak(first.get("wide_peak", ""))
        if chrom not in CHROMS:
            continue
        try:
            q = float(first["q"])
        except (TypeError, ValueError):
            q = np.nan
        rows.append(dict(
            cohort=cohort, type=str(ptype),
            cytoband=str(first.get("cytoband", "")),
            chrom=chrom, start=start, end=end, q=q,
            neg_log10_q=float(-np.log10(q)) if np.isfinite(q) and q > 0
            else np.nan,
            gene=_pick_gene(list(g["gene"]))))
    return pd.DataFrame(rows, columns=cols)


def run_module_p11(cfg: dict[str, Any], data: dict[str, dict]) -> dict[str, Any]:
    fcfg = _focal_cfg(cfg)
    scores: dict[str, pd.DataFrame] = {}
    peak_frames: list[pd.DataFrame] = []
    for cid, d in data.items():
        s = d.get("scores")
        if s is None or len(s) == 0:
            raise RuntimeError(f"[P11] cohort {cid} has no scores.gistic")
        scores[cid] = s[["chrom", "start", "end", "type", "neg_log10_q",
                          "g_score", "frequency", "avg_amplitude"]].copy()
        peak_frames.append(_peak_records(d.get("genes"), cid))

    peaks = (pd.concat(peak_frames, ignore_index=True) if peak_frames
             else pd.DataFrame())
    peaks.to_csv(out_table(cfg, "P11_focal_peaks.tsv"), sep="\t", index=False)

    q_thr = fcfg["q_threshold"]
    for cid, s in scores.items():
        n_amp = int((s["type"] == "Amp").sum())
        n_del = int((s["type"] == "Del").sum())
        print(f"  [P11] {cid}: {n_amp} amp / {n_del} del score segments")
    if not peaks.empty:
        n_sig = int((peaks["q"] < q_thr).sum())
        print(f"  [P11] focal peaks with q<{q_thr:g}: {n_sig}")

    return dict(scores=scores, peaks=peaks, q_threshold=q_thr)
