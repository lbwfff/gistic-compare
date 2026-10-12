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


# ---------------------------------------------------------------------------
# synthetic ASCAT pipeline inputs (Phase 2 R-wrapper tests)
# ---------------------------------------------------------------------------

def _marker_grid(n_chr: int = 3, per_chr: int = 800, step: int = 100_000
                 ) -> pd.DataFrame:
    rows = []
    for c in range(1, n_chr + 1):
        start = 1_000_000
        for k in range(per_chr):
            rows.append((c, start + k * step))
    return pd.DataFrame(rows, columns=["chr", "pos"])


def _allele_fraction(m_minor: int, total: int, purity: float) -> float:
    """Expected BAF for a germline-heterozygous SNP under tumor purity
    ``purity`` and tumor minor/total copy numbers."""
    f = (purity * m_minor + (1 - purity) * 1.0) \
        / (purity * total + (1 - purity) * 2.0)
    return float(np.clip(max(f, 1.0 - f), 0.0, 1.0))


def _logr(total: int, purity: float, gamma: float = 0.55) -> float:
    """ASCAT-convention logR: the observed logR is gamma-compressed
    (ASCAT's copy-number model is 2^(logR/gamma)); gamma must match the
    value passed to ascat.runAscat."""
    return float(gamma * np.log2((purity * total + (1 - purity) * 2.0) / 2.0))


# per-chromosome-arm truth: (total, minor); arms = chr1p, chr1q, chr2p, ...
DEFAULT_ARM_TRUTH_WGD = [(4, 2), (4, 2), (5, 2), (4, 2), (4, 0), (4, 2)]
DEFAULT_ARM_TRUTH_NO_WGD = [(2, 1), (2, 1), (3, 1), (2, 1), (1, 0), (2, 1)]


def synth_ascat_inputs(root, cohort: str, n_samples: int = 4,
                       wgd: bool = True, purity: float = 0.75,
                       seed: int = 42, sample_names: list[str] | None = None,
                       n_markers_per_chr: int = 800,
                       gamma: float = 0.55,
                       het_calls: bool = False) -> pd.DataFrame:
    """Write per-sample marker-level logR/BAF TSVs (chr, pos, value) plus
    one pooled-germline pair under ``root``; returns the manifest rows
    (without cohort-ordering guarantees)."""
    rng = np.random.default_rng(seed)
    grid = _marker_grid(per_chr=n_markers_per_chr)
    n = len(grid)

    # pooled germline reference: mostly heterozygous, a few homozygous
    gl = rng.normal(0.0, 0.03, n)
    hom = rng.random(n) < 0.03
    gbaf = np.where(hom, rng.choice([0.01, 0.99], n),
                    np.clip(rng.normal(0.5, 0.02, n), 0, 1))
    root.mkdir(parents=True, exist_ok=True)
    glr_path = root / "germline_logr.tsv"
    gbaf_path = root / "germline_baf.tsv"
    pd.DataFrame({"chr": grid["chr"], "pos": grid["pos"],
                  "value": gl.round(4)}).to_csv(glr_path, sep="\t",
                                                index=False)
    pd.DataFrame({"chr": grid["chr"], "pos": grid["pos"],
                  "value": gbaf.round(4)}).to_csv(gbaf_path, sep="\t",
                                                  index=False)
    # probe-name grid, per cohort dir — the R wrapper resolves it as
    # <manifest_dir>/<cohort>/markers.tsv.gz (keyed by probe name for the
    # GC-content reference match)
    cdir = root / cohort
    cdir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"name": [f"SNP_A-{1000000 + i}" for i in range(n)],
                  "chr": grid["chr"], "pos": grid["pos"]}).to_csv(
        cdir / "markers.tsv.gz", sep="\t", index=False, compression="gzip")
    het_path = None
    if het_calls:
        # birdseed-style genotype calls derived from the germline BAF:
        # BAF ~ 0.5 -> AB (het), near 0/1 -> AA/BB (hom)
        call = np.where(np.abs(gbaf - 0.5) < 0.15, "AB",
                        np.where(gbaf < 0.5, "AA", "BB"))
        het_path = root / "germline_het.tsv"
        pd.DataFrame({"chr": grid["chr"], "pos": grid["pos"],
                      "call": call}).to_csv(het_path, sep="\t", index=False)

    truth = (DEFAULT_ARM_TRUTH_WGD if wgd else DEFAULT_ARM_TRUTH_NO_WGD)
    arms_per_chr = 2
    states = []
    for c in range(1, grid["chr"].max() + 1):
        for a in range(arms_per_chr):
            states.append(truth[(c - 1) * arms_per_chr + a])
    idx = np.arange(n)
    markers_per_arm = n_markers_per_chr // arms_per_chr
    arm_id = np.minimum((idx // n_markers_per_chr) * arms_per_chr
                        + (idx % n_markers_per_chr) // markers_per_arm,
                        len(states) - 1)

    names = (list(sample_names) if sample_names
             else [f"{cohort.upper()}{i + 1:02d}" for i in range(n_samples)])
    rows = []
    for i, s in enumerate(names):
        lr, bf = [], []
        for j in range(n):
            total, minor = states[arm_id[j]]
            lr.append(_logr(total, purity, gamma) + rng.normal(0, 0.02))
            bf.append(_allele_fraction(minor, total, purity)
                      + rng.normal(0, 0.01))
        lr_path = root / f"{s}_logr.tsv"
        baf_path = root / f"{s}_baf.tsv"
        pd.DataFrame({"chr": grid["chr"], "pos": grid["pos"],
                      "value": np.round(lr, 4)}).to_csv(lr_path, sep="\t",
                                                        index=False)
        pd.DataFrame({"chr": grid["chr"], "pos": grid["pos"],
                      "value": np.round(np.clip(bf, 0, 1), 4)}).to_csv(
            baf_path, sep="\t", index=False)
        row = dict(sample=s, cohort=cohort,
                   logr=lr_path.name, baf=baf_path.name,
                   germline_logr=glr_path.name,
                   germline_baf=gbaf_path.name)
        if het_path is not None:
            row["het_calls"] = het_path.name
        rows.append(row)
    return pd.DataFrame(rows)
