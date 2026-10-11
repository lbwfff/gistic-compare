"""Generate a tiny synthetic GISTIC 2.0 output bundle for tests.

The files follow the real GISTIC layout closely enough for every parser in
``gistic_compare.io_gistic`` and for modules P4/P5/GS (which need a marker
matrix ``D.cap1.5.mat`` in HDF5/MATLAB-v7.3 form).  Everything is deliberately
small: 2 cohorts, 8 samples each, a handful of peaks, and copy-number signal
on chr1 so the regional tests have something to find.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

LESIONS_COLS = [
    "Unique Name", "Descriptor", "Wide Peak Limits", "Peak Limits",
    "Region Limits", "q values",
    "Residual q values after removing segments shared with higher peaks",
    "Broad or Focal", "Amplitude Threshold",
]

# (type, descriptor, chrom, start, end, q, per-cohort call vector name)
def peaks_for(cohort: str) -> list[dict]:
    """A handful of peaks; some shared, some cohort-specific."""
    common = [
        dict(utype="amplification", desc="5p15.2", chrom="chr5",
             start=1_500_000, end=10_000_000, q=0.03,
             calls={"A": [2, 2, 2, 1, 0, 0, 0, 0],
                    "B": [2, 2, 0, 0, 0, 0, 0, 0]}),
        dict(utype="deletion", desc="3p14.2", chrom="chr3",
             start=50_300_000, end=61_000_000, q=0.08,
             calls={"A": [1, 0, 0, 0, 0, 0, 0, 0],
                    "B": [2, 1, 1, 0, 0, 0, 0, 0]}),
        dict(utype="amplification", desc="7p11.2", chrom="chr7",
             start=55_000_000, end=55_900_000, q=0.20,
             calls={"A": [0, 1, 0, 0, 0, 0, 0, 0],
                    "B": [0, 0, 2, 0, 0, 0, 0, 0]}),
    ]
    if cohort == "A":
        common += [
            dict(utype="amplification", desc="1q21.2", chrom="chr1",
                 start=147_000_000, end=150_000_000, q=0.01,
                 calls={"A": [2, 2, 2, 2, 1, 0, 0, 0],
                        "B": [0, 0, 0, 0, 0, 0, 0, 0]}),
            dict(utype="deletion", desc="17p13.1", chrom="chr17",
                 start=7_000_000, end=7_600_000, q=0.30,
                 calls={"A": [1, 1, 0, 0, 0, 0, 0, 0],
                        "B": [0, 0, 0, 0, 0, 0, 0, 0]}),
        ]
    else:
        common += [
            dict(utype="deletion", desc="9p21.3", chrom="chr9",
                 start=21_000_000, end=21_900_000, q=0.02,
                 calls={"A": [0, 0, 0, 0, 0, 0, 0, 0],
                        "B": [2, 2, 1, 1, 0, 0, 0, 0]}),
            dict(utype="amplification", desc="8q24.21", chrom="chr8",
                 start=128_000_000, end=128_800_000, q=0.12,
                 calls={"A": [0, 0, 0, 0, 0, 0, 0, 0],
                        "B": [1, 0, 0, 0, 0, 0, 0, 0]}),
        ]
    return common


def write_all_lesions(path: Path, cohort: str, samples: list[str]) -> None:
    rows = []
    # calls are keyed by "A"/"B" design slots; other cohorts reuse A's vector
    call_key = cohort if cohort in ("A", "B") else "A"
    for pk in peaks_for(cohort):
        coord = (f"{pk['chrom']}:{pk['start']}-{pk['end']}(probes 1:2)")
        summary = ([pk["utype"] + f"({pk['desc']})", pk["desc"], coord, coord,
                    coord, f"{pk['q']:.4g}", f"{pk['q'] / 2:.4g}",
                    "Broad or Focal", "Amplitude Threshold"])
        calls = pk["calls"][call_key]
        summary += [str(c) for c in calls]
        # CN values row: continuous per-sample copy change; altered samples
        # carry the direction-specific value so P7's peak layer has signal
        sign = 1.0 if pk["utype"] == "amplification" else -1.0
        cn_vals = [f"{sign * (0.6 if c > 0 else 0.0):.4g}" for c in calls]
        cn = ([pk["utype"] + f"({pk['desc']}) - CN values", pk["desc"], coord,
               coord, coord, "NA", "NA", "NA", "NA"] + cn_vals)
        rows.append(summary)
        rows.append(cn)
    with open(path, "w") as fh:
        fh.write("\t".join(LESIONS_COLS + samples) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")


def write_genes(path: Path, cohort: str, ptype: str) -> None:
    # transposed layout: row0 cytoband, row1 q, row2 residual q, row3 wide
    # peak, rows 4+ genes; first column is the row label
    pks = [p for p in peaks_for(cohort)
           if p["utype"] == ("amplification" if ptype == "amp" else "deletion")]
    if not pks:
        path.write_text("")
        return
    gene_pool = (["EGFR", "MYC", "CCND1", "NKX2-1"] if ptype == "amp"
                 else ["FHIT", "CDKN2A", "WWOX", "TP53"])
    rows = [["#"] + [p["desc"] for p in pks],
            ["q"] + [f"{p['q']:.4g}" for p in pks],
            ["residual_q"] + [f"{p['q'] / 2:.4g}" for p in pks],
            ["wide_peak"] + [f"{p['chrom']}:{p['start']}-{p['end']}"
                             for p in pks]]
    for gi, g in enumerate(gene_pool):
        rows.append([f"gene{gi}"] + [g if (gi + ci) % 2 == 0 else ""
                                     for ci, _p in enumerate(pks)])
    with open(path, "w") as fh:
        for r in rows:
            fh.write("\t".join(r) + "\n")


SCORES_COLS = ["Type", "Chromosome", "Start", "End", "-log10(q-value)",
               "G-score", "average amplitude", "frequency"]


def write_scores(path: Path, cohort: str) -> None:
    with open(path, "w") as fh:
        fh.write("\t".join(SCORES_COLS) + "\n")
        for pk in peaks_for(cohort):
            t = "Amp" if pk["utype"] == "amplification" else "Del"
            fh.write("\t".join([
                t, pk["chrom"].replace("chr", ""), str(pk["start"]),
                str(pk["end"]), "2.5", "0.75", "0.4", "0.3"]) + "\n")
        # a few tiled rows on chr1/chr2 for module B
        for i in range(10):
            for t, sign in (("Amp", 1), ("Del", -1)):
                fh.write("\t".join([
                    t, "1" if i % 2 == 0 else "2", str(i * 1_000_000),
                    str((i + 1) * 1_000_000 - 1), f"{0.5 + sign * 0.01 * i:.2f}",
                    f"{0.1 * i:.2f}", "0.2", f"{0.05 * i:.2f}"]) + "\n")


def write_seg(path: Path, cohort: str, samples: list[str]) -> None:
    """chr1 5-10 Mb carries the between-cohort focal signal for P4/P5/P7;
    chr2q is broadly amplified in every A sample (broad signal for P6)."""
    with open(path, "w") as fh:
        fh.write("Sample\tChromosome\tStart\tEnd\tNum_Probes\tSegment_Mean\n")
        for si, s in enumerate(samples):
            for c, size in [(1, 249250621), (2, 243199373), (3, 198022430)]:
                if c == 1:
                    # background + the differential region
                    fh.write(f"{s}\tchr1\t0\t5_000_000\t100\t0.0\n".replace("_", ""))
                    mid = 0.4 if cohort == "A" else -0.3
                    fh.write(f"{s}\tchr1\t5000000\t10000000\t200\t{mid}\n")
                    fh.write(f"{s}\tchr1\t10000000\t{size}\t800\t0.0\n")
                elif c == 2:
                    # whole chromosome at 0.5 in A (broad 2q amp), 0 in B
                    broad = 0.5 if cohort == "A" else 0.0
                    fh.write(f"{s}\tchr2\t0\t{size}\t300\t{broad}\n")
                else:
                    fh.write(f"{s}\tchr{c}\t0\t{size}\t300\t0.0\n")


def write_seg_counts(path: Path, samples: list[str]) -> None:
    with open(path, "w") as fh:
        fh.write("sample\tsegment_count\tincluded\n")
        for s in samples:
            fh.write(f"{s}\t6\tyes\n")


def write_regions_bed(path: Path, cohort: str) -> None:
    with open(path, "w") as fh:
        fh.write('track name="Any-AP Any-DP"\n')
        for pk in peaks_for(cohort):
            t = "AP" if pk["utype"] == "amplification" else "DP"
            fh.write(f"{pk['chrom']}\t{pk['start']}\t{pk['end']}\t"
                     f"Any-{t}\n")


def write_marker_mat(path: Path, cohort: str, samples: list[str]) -> None:
    """Write a D.cap1.5.mat (MATLAB v7.3 / HDF5 SegArray) with signal."""
    import h5py

    rng = np.random.default_rng(42)
    # markers every 1 Mb on chr1-3 (small subset of the real thing)
    chrom_sizes = {1: 249250621, 2: 243199373, 3: 198022430}
    pos, chrn = [], []
    for c, size in chrom_sizes.items():
        for p in range(500_000, size, 1_000_000):
            pos.append(p)
            chrn.append(c)
    M = len(pos)
    pos = np.asarray(pos, dtype=np.int64)
    chrn = np.asarray(chrn, dtype=np.float64)

    # segment boundaries: a coarse tiling plus exact breakpoints around the
    # chr1 5-10 Mb region so the differential segment covers its markers
    region_idx = np.flatnonzero((pos >= 5_000_000) & (pos <= 10_000_000)
                                & (chrn == 1))
    lo, hi = int(region_idx[0]), int(region_idx[-1]) + 1
    breaks = np.linspace(0, M, 13).astype(np.int64)
    breaks = np.array(sorted(set(breaks.tolist()) | {lo, hi}), dtype=np.int64)
    nseg_per_sample = len(breaks) - 1
    region_seg = int(np.flatnonzero(breaks == lo)[0])
    ir = breaks[:-1].astype(np.float64)

    sdesc = np.zeros((12, len(samples)))
    for j, s in enumerate(samples):
        for k, ch in enumerate(s[:12]):
            sdesc[k, j] = ord(ch)

    with h5py.File(path, "w") as h:
        d = h.create_group("D")
        d.create_dataset("pos", data=pos.reshape(1, -1))
        d.create_dataset("chrn", data=chrn.reshape(1, -1))
        d.create_dataset("sdesc", data=sdesc)
        refs = h.create_group("#refs#")
        pg = refs.create_group("p")
        # ir/jc are 1-D in real GISTIC v7.3 files; jc has nsamples+1 entries
        jc = np.arange(len(samples) + 1, dtype=np.float64) * nseg_per_sample
        pg.create_dataset("jc", data=jc)
        pg.create_dataset("ir", data=np.tile(ir, len(samples)))
        # per-sample segment values: A up / B down on the chr1 region markers
        vals = []
        for j in range(len(samples)):
            sv = rng.normal(0.0, 0.05, nseg_per_sample)
            sv[region_seg] = 0.6 if cohort == "A" else -0.6
            vals.append(sv)
        o = np.concatenate(vals).astype(np.float64)
        refs.create_dataset("o", data=o.reshape(1, -1))


def write_broad_values_by_arm(path: Path, cohort: str,
                              samples: list[str]) -> None:
    """Arm-level broad calls (0/1/2 convention: 1=amp, 2=del).

    Cohort A: 2q amp in every sample (matches the chr2 seg signal) and 1q
    amp in 5/8 (Fisher-differentiable); cohort B: 3q del in 3/8.
    """
    arms = [f"{n}{s}" for n in (1, 2, 3) for s in ("p", "q")]
    spec: dict[tuple[str, str], int] = {}
    if cohort == "A":
        for s in samples:
            spec[(s, "2q")] = 1
        for s in samples[:5]:
            spec[(s, "1q")] = 1
    else:
        for s in samples[:3]:
            spec[(s, "3q")] = 2
    with open(path, "w") as fh:
        fh.write("\t".join(["sample"] + arms) + "\n")
        for s in samples:
            vals = []
            for a in arms:
                v = spec.get((s, a), 0)
                vals.append({0: "0", 1: "1", 2: "2"}[v])
            fh.write("\t".join([s] + vals) + "\n")


def write_arm_qv(path: Path, cohort: str) -> None:
    """Minimal arm-level q-value table (arm + q), P6 display-only input."""
    arms = ["1p", "1q", "2p", "2q", "3p", "3q"]
    with open(path, "w") as fh:
        fh.write("arm\tq\n")
        for i, a in enumerate(arms):
            fh.write(f"{a}\t{0.01 * (i + 1):.4g}\n")


def make_cohort_dir(root: Path, cohort: str, prefix: str) -> Path:
    d = root / cohort
    d.mkdir(parents=True, exist_ok=True)
    samples = [f"{cohort}{i}" for i in range(1, 9)]
    write_all_lesions(d / f"{prefix}.all_lesions.conf_90.txt", cohort, samples)
    write_genes(d / f"{prefix}.amp_genes.conf_90.txt", cohort, "amp")
    write_genes(d / f"{prefix}.del_genes.conf_90.txt", cohort, "del")
    write_scores(d / f"{prefix}.scores.gistic", cohort)
    write_seg(d / f"{prefix}.tumors.seg", cohort, samples)
    write_seg_counts(d / f"{prefix}.sample_seg_counts.txt", samples)
    write_regions_bed(d / f"{prefix}.regions_track.conf_90.bed", cohort)
    write_marker_mat(d / f"{prefix}.D.cap1.5.mat", cohort, samples)
    write_broad_values_by_arm(d / f"{prefix}.broad_values_by_arm.txt",
                              cohort, samples)
    write_arm_qv(d / f"{prefix}.amp_qv.conf_90.txt", cohort)
    write_arm_qv(d / f"{prefix}.del_qv.conf_90.txt", cohort)
    return d
