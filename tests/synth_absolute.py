"""Generate synthetic absolute-copy-number tables (FACETS-output-like).

Two cohorts with known truth, used to test the ABS comparison suite
without R/FACETS:

  cohort A: WGD+ fraction ~0.60, purity ~ N(0.75, 0.08), ploidy ~ N(3.0, 0.4)
  cohort B: WGD+ fraction ~0.15, purity ~ N(0.65, 0.08), ploidy ~ N(2.2, 0.4)

Writes, per cohort directory:

  absolute_cn_<cid>.tsv    per-sample table (schema of plan_absolute_cn.md)
  absolute_seg_<cid>.tsv   allele-specific segments
                           (sample, chrom, start, end, cn_total, cn_major,
                            cn_minor)

Segments tile the real hg19 chromosome sizes so arm-level statistics can be
computed against the bundled cytoband table.  WGD+ samples carry mostly
minor CN >= 1 / total CN >= 3 segments (the duplication signature); WGD-
samples are mostly neutral with sporadic gains/losses.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from gistic_compare.genomes import HG19_CHROM_SIZES

_CHROMS = [f"chr{i}" for i in range(1, 23)]

# cohort -> (wgd_frac, purity_mean, ploidy_mean)
DEFAULT_TRUTH = {
    "a": dict(wgd_frac=0.60, purity_mean=0.75, ploidy_mean=3.0),
    "b": dict(wgd_frac=0.15, purity_mean=0.65, ploidy_mean=2.2),
}


def _sample_segments(rng: np.random.Generator, is_wgd: bool) -> pd.DataFrame:
    rows = []
    for c in _CHROMS:
        size = HG19_CHROM_SIZES[c]
        n_seg = int(rng.integers(3, 6))
        cuts = np.sort(rng.choice(np.arange(1, size), size=n_seg - 1,
                                  replace=False))
        bounds = np.r_[0, cuts, size]
        for k in range(n_seg):
            if is_wgd and rng.random() < 0.85:
                total = int(rng.choice([3, 3, 4, 4, 5]))
                minor = int(rng.choice([1, 1, 2])) if total >= 4 else 1
                minor = min(minor, total // 2)
            elif rng.random() < 0.70:
                total, minor = 2, 1
            elif rng.random() < 0.55:
                total, minor = 3, 1
            elif rng.random() < 0.6:
                total, minor = 1, 1
            else:
                total, minor = 0, 0
            rows.append(dict(
                chrom=c, start=int(bounds[k]), end=int(bounds[k + 1]),
                cn_total=int(total), cn_major=int(total - minor),
                cn_minor=int(minor)))
    return pd.DataFrame(rows)


def synth_absolute_tables(root, cohort: str, n_samples: int = 40,
                          truth: dict | None = None, seed: int = 42,
                          sample_names: list[str] | None = None
                          ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Write absolute_cn_<cohort>.tsv + absolute_seg_<cohort>.tsv under
    ``root``; returns (per-sample table, segments table)."""
    rng = np.random.default_rng(seed)
    tr = dict(DEFAULT_TRUTH.get(cohort, DEFAULT_TRUTH["b"]), **(truth or {}))
    names = (list(sample_names) if sample_names
             else [f"{cohort.upper()}{i + 1:02d}"
                   for i in range(n_samples)])
    cn_rows, seg_frames = [], []
    for i in range(len(names)):
        sample = names[i]
        is_wgd = bool(rng.random() < tr["wgd_frac"])
        purity = float(np.clip(rng.normal(tr["purity_mean"],
                                          tr.get("purity_sd", 0.08)), 0.15,
                               0.95))
        ploidy = float(max(1.2, rng.normal(
            tr["ploidy_mean"] + (0.8 if is_wgd else 0.0), 0.4)))
        loh = float(np.clip((0.45 if is_wgd else 0.10)
                            + rng.normal(0, 0.08), 0.0, 1.0))
        cnloh = float(np.clip(rng.normal(0.05, 0.03), 0.0, 1.0))
        wfrac = float(np.clip((0.45 if is_wgd else 0.05)
                              + rng.normal(0, 0.08), 0.0, 1.0))
        segs = _sample_segments(rng, is_wgd)
        segs.insert(0, "sample", sample)
        seg_frames.append(segs)
        dup_len = segs.loc[(segs["cn_minor"] >= 1)
                           & (segs["cn_total"] >= 3)]
        n_dup = int(len(dup_len))
        dup_frac = float((dup_len["end"] - dup_len["start"]).sum()
                         / (segs["end"] - segs["start"]).sum())
        cn_rows.append(dict(
            sample=sample,
            purity=round(purity, 4),
            ploidy_abs=round(ploidy, 3),
            wgd_status=int(is_wgd) + int(is_wgd and rng.random() < 0.1),
            wgd_frac=round(dup_frac, 4),
            loh_frac=round(loh, 4),
            cnloh_frac=round(cnloh, 4),
            n_segments_allelic=int(len(segs)),
            flag_low_purity=bool(purity < 0.30),
            n_dup_segments=int(n_dup),
        ))
    cn = pd.DataFrame(cn_rows)
    segs = pd.concat(seg_frames, ignore_index=True)
    root.mkdir(parents=True, exist_ok=True)
    cn.to_csv(root / f"absolute_cn_{cohort}.tsv", sep="\t", index=False)
    segs.to_csv(root / f"absolute_seg_{cohort}.tsv", sep="\t", index=False)
    return cn, segs
