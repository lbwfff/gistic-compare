"""Interval overlap, one-to-one peak matching, and permutation nulls."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .config import CHROMS, HG19_CHROM_SIZES

# ---------------------------------------------------------------------------
# basic interval algebra
# ---------------------------------------------------------------------------

def overlap_len(a_start, a_end, b_start, b_end) -> np.ndarray:
    return np.maximum(0, np.minimum(a_end, b_end) - np.maximum(a_start, b_start))


def filter_peaks(peaks: pd.DataFrame, q_thresh: float, ptype: str) -> pd.DataFrame:
    m = (peaks["type"] == ptype) & (peaks["q"] < q_thresh)
    return peaks.loc[m].reset_index(drop=True)


# ---------------------------------------------------------------------------
# greedy one-to-one matching by reciprocal overlap
# ---------------------------------------------------------------------------
#
# reciprocal modes:
#   both   — require rec_a >= t AND rec_b >= t          (same-scale peaks)
#   either — require max(rec_a, rec_b) >= t             (== shorter peak covered >= t;
#            solves broad-contains-focal size asymmetry; DEFAULT moderate tier)
#   any    — only min_bp (loose lower bound)

@dataclass
class MatchResult:
    pairs: pd.DataFrame      # a_peak_id, b_peak_id, overlap_bp, reciprocal_a, reciprocal_b
    only_a: pd.DataFrame
    only_b: pd.DataFrame
    n_a: int
    n_b: int


def _passes_rec(rec_a: float, rec_b: float, reciprocal: float, mode: str) -> bool:
    if mode == "both":
        return rec_a >= reciprocal and rec_b >= reciprocal
    if mode == "either":
        return max(rec_a, rec_b) >= reciprocal
    if mode == "any":
        return True
    raise ValueError(f"unknown reciprocal mode: {mode!r}")


def match_peaks(
    a: pd.DataFrame,
    b: pd.DataFrame,
    min_bp: int = 50,
    reciprocal: float = 0.5,
    mode: str = "either",
) -> MatchResult:
    """One-to-one greedy matching of peaks on genomic overlap.

    Pairs are ranked by overlap_bp (desc), then by max(reciprocal) (desc),
    and assigned greedily so each peak is used at most once.
    `a`/`b` need columns: peak_id, chrom, start, end.
    """
    a = a.reset_index(drop=True)
    b = b.reset_index(drop=True)
    cand = []
    b_by_chrom: dict[str, list[int]] = {}
    for j, rb in b.iterrows():
        b_by_chrom.setdefault(rb["chrom"], []).append(j)

    for i, ra in a.iterrows():
        for j in b_by_chrom.get(ra["chrom"], ()):
            rb = b.loc[j]
            ov = min(ra["end"], rb["end"]) - max(ra["start"], rb["start"])
            if ov < min_bp:
                continue
            la = max(1, ra["end"] - ra["start"])
            lb = max(1, rb["end"] - rb["start"])
            ra_rec = ov / la
            rb_rec = ov / lb
            if not _passes_rec(ra_rec, rb_rec, reciprocal, mode):
                continue
            cand.append((ov, max(ra_rec, rb_rec), min(ra_rec, rb_rec),
                         i, j, ra_rec, rb_rec))

    cand.sort(key=lambda t: (-t[0], -t[1], -t[2]))
    used_a: set[int] = set()
    used_b: set[int] = set()
    matched_a_ids: set[str] = set()
    matched_b_ids: set[str] = set()
    rows = []
    for ov, _, _, i, j, ra_rec, rb_rec in cand:
        if i in used_a or j in used_b:
            continue
        used_a.add(i)
        used_b.add(j)
        matched_a_ids.add(a.loc[i, "peak_id"])
        matched_b_ids.add(b.loc[j, "peak_id"])
        rows.append(dict(
            a_peak_id=a.loc[i, "peak_id"],
            b_peak_id=b.loc[j, "peak_id"],
            overlap_bp=int(ov),
            reciprocal_a=float(ra_rec),
            reciprocal_b=float(rb_rec),
            a_cytoband=a.loc[i, "cytoband"] if "cytoband" in a.columns else "",
            b_cytoband=b.loc[j, "cytoband"] if "cytoband" in b.columns else "",
        ))

    pairs = pd.DataFrame(rows, columns=[
        "a_peak_id", "b_peak_id", "overlap_bp", "reciprocal_a", "reciprocal_b",
        "a_cytoband", "b_cytoband",
    ])
    only_a = a.loc[~a["peak_id"].isin(matched_a_ids)].reset_index(drop=True)
    only_b = b.loc[~b["peak_id"].isin(matched_b_ids)].reset_index(drop=True)
    return MatchResult(pairs=pairs, only_a=only_a, only_b=only_b,
                       n_a=len(a), n_b=len(b))


def jaccard(a_n: int, b_n: int, n_pairs: int) -> float:
    union = a_n + b_n - n_pairs
    return n_pairs / union if union else float("nan")


# ---------------------------------------------------------------------------
# cytoband (string) matching — sensitivity analysis
# ---------------------------------------------------------------------------

def match_by_cytoband(a: pd.DataFrame, b: pd.DataFrame) -> MatchResult:
    """Match peaks whose cytoband strings are identical (per type assumed filtered)."""
    a = a.reset_index(drop=True)
    b = b.reset_index(drop=True)
    rows = []
    used_b: set[int] = set()
    for i, ra in a.iterrows():
        hits = b.index[(b["cytoband"] == ra["cytoband"]) & (~b.index.isin(used_b))]
        if len(hits) == 0:
            continue
        j = hits[0]
        used_b.add(int(j))
        rows.append(dict(
            a_peak_id=ra["peak_id"], b_peak_id=b.loc[j, "peak_id"],
            overlap_bp=-1, reciprocal_a=np.nan, reciprocal_b=np.nan,
            a_cytoband=ra["cytoband"], b_cytoband=b.loc[j, "cytoband"],
        ))
    pairs = pd.DataFrame(rows, columns=[
        "a_peak_id", "b_peak_id", "overlap_bp", "reciprocal_a", "reciprocal_b",
        "a_cytoband", "b_cytoband",
    ])
    matched_a = set(pairs["a_peak_id"]) if len(pairs) else set()
    matched_b = set(pairs["b_peak_id"]) if len(pairs) else set()
    only_a = a.loc[~a["peak_id"].isin(matched_a)].reset_index(drop=True)
    only_b = b.loc[~b["peak_id"].isin(matched_b)].reset_index(drop=True)
    return MatchResult(pairs=pairs, only_a=only_a, only_b=only_b,
                       n_a=len(a), n_b=len(b))


# ---------------------------------------------------------------------------
# permutation null for peak-set overlap
# ---------------------------------------------------------------------------
#
# Null: re-place both peak sets on the genome preserving (i) exact length
# multiset and (ii) per-chromosome peak counts, with uniform start position
# on each chromosome (peaks must fit). Match with the same rule as observed.
# p-value is Monte-Carlo corrected: p = (#{J_k >= J_obs} + 1) / (n_perm + 1).

def _chrom_counts(peaks: pd.DataFrame) -> dict[str, int]:
    return {str(c): int(n) for c, n in peaks["chrom"].value_counts().items()}


def _random_peak_set(lengths: np.ndarray, rng: np.random.Generator,
                     chrom_counts: dict[str, int] | None = None,
                     ) -> pd.DataFrame:
    """Draw a random peak set.

    If chrom_counts is given, preserve that per-chromosome count vector;
    otherwise choose chromosome ∝ hg19 length. Lengths are shuffled onto
    the chromosome slots (exact length multiset retained when
    len(lengths) == sum(chrom_counts.values())).
    """
    n_slots = sum(chrom_counts.values()) if chrom_counts else len(lengths)
    lens = np.resize(np.asarray(lengths, dtype=int), n_slots)
    rng.shuffle(lens)

    if chrom_counts:
        slots: list[str] = []
        for chrom, cnt in chrom_counts.items():
            slots.extend([chrom] * int(cnt))
        # slots order follows dict insertion (value_counts order); shuffle pairing
        # already via lens shuffle — also shuffle slot order so lengths aren't
        # systematically tied to one chrom's block order
        slot_arr = np.array(slots, dtype=object)
        rng.shuffle(slot_arr)
        chroms_iter = slot_arr.tolist()
    else:
        all_chroms = CHROMS
        sizes = np.array([HG19_CHROM_SIZES[c] for c in all_chroms], dtype=float)
        probs = sizes / sizes.sum()
        idx = rng.choice(len(all_chroms), size=n_slots, p=probs)
        chroms_iter = [all_chroms[i] for i in idx]

    rows = []
    for chrom, L in zip(chroms_iter, lens):
        chrom_size = HG19_CHROM_SIZES[chrom]
        L = int(min(int(L), chrom_size))
        max_start = chrom_size - L
        start = int(rng.integers(0, max_start + 1)) if max_start > 0 else 0
        rows.append((chrom, start, start + L))
    return pd.DataFrame(rows, columns=["chrom", "start", "end"])


def _count_matches(a: pd.DataFrame, b: pd.DataFrame, min_bp: int,
                   reciprocal: float, mode: str = "either") -> int:
    """Count one-to-one matches without building full result frames (speed)."""
    b_by_chrom: dict[str, list[int]] = {}
    for j, rb in b.iterrows():
        b_by_chrom.setdefault(rb["chrom"], []).append(j)
    cand = []
    for i, ra in a.iterrows():
        for j in b_by_chrom.get(ra["chrom"], ()):
            rb = b.loc[j]
            ov = min(ra["end"], rb["end"]) - max(ra["start"], rb["start"])
            if ov < min_bp:
                continue
            ra_rec = ov / max(1, ra["end"] - ra["start"])
            rb_rec = ov / max(1, rb["end"] - rb["start"])
            if not _passes_rec(ra_rec, rb_rec, reciprocal, mode):
                continue
            cand.append((-ov, -max(ra_rec, rb_rec), -min(ra_rec, rb_rec), i, j))
    cand.sort()
    used_a: set[int] = set()
    used_b: set[int] = set()
    n = 0
    for _, _, _, i, j in cand:
        if i in used_a or j in used_b:
            continue
        used_a.add(i)
        used_b.add(j)
        n += 1
    return n


def permutation_jaccard(
    a: pd.DataFrame,
    b: pd.DataFrame,
    n_perm: int = 10000,
    min_bp: int = 50,
    reciprocal: float = 0.5,
    mode: str = "either",
    seed: int = 42,
) -> dict[str, float]:
    """Randomly re-place peaks (length multiset + per-chrom counts preserved)
    and recompute one-to-one match Jaccard null.

    p-value: (#{J_k >= J_obs} + 1) / (n_perm + 1)  [Phipson & Smyth].

    Returns summary stats plus ``null_jaccards`` (np.ndarray of length n_perm).
    """
    observed = match_peaks(a, b, min_bp=min_bp, reciprocal=reciprocal, mode=mode)
    obs_n = len(observed.pairs)
    union = observed.n_a + observed.n_b - obs_n
    obs_j = obs_n / union if union else float("nan")

    rng = np.random.default_rng(seed)
    len_a = (a["end"] - a["start"]).to_numpy()
    len_b = (b["end"] - b["start"]).to_numpy()
    counts_a = _chrom_counts(a)
    counts_b = _chrom_counts(b)

    null_n = np.empty(n_perm, dtype=int)
    for k in range(n_perm):
        ra = _random_peak_set(len_a, rng, chrom_counts=counts_a)
        rb = _random_peak_set(len_b, rng, chrom_counts=counts_b)
        null_n[k] = _count_matches(ra, rb, min_bp, reciprocal, mode)

    null_union = len(a) + len(b) - null_n
    with np.errstate(divide="ignore", invalid="ignore"):
        null_j = np.where(null_union > 0, null_n / null_union, np.nan)
    null_j = null_j[~np.isnan(null_j)]

    if np.isnan(obs_j):
        pvalue = np.nan
    else:
        k_ge = int((null_j >= obs_j).sum())
        pvalue = float((k_ge + 1) / (len(null_j) + 1))
    return dict(
        n_a=int(observed.n_a), n_b=int(observed.n_b), n_pairs=int(obs_n),
        observed_jaccard=float(obs_j),
        null_mean=float(np.mean(null_j)) if len(null_j) else np.nan,
        null_std=float(np.std(null_j)) if len(null_j) else np.nan,
        null_median=float(np.median(null_j)) if len(null_j) else np.nan,
        null_q95=float(np.quantile(null_j, 0.95)) if len(null_j) else np.nan,
        pvalue=pvalue,
        n_perm=int(n_perm),
        min_bp=int(min_bp),
        reciprocal=float(reciprocal),
        reciprocal_mode=mode,
        null_jaccards=null_j,
    )
