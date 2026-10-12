#!/usr/bin/env python3
"""Split cohort-wide SNP6 Birdseed outputs into per-sample ABS inputs.

Inputs (per cohort, as produced by the apt-probeset-genotype / TCGA-style
chain from raw CEL files):

  <lrr_baf>   wide TSV: Name, Chr, Position, then per-sample pairs of
              ("<sample>.CEL.Log R Ratio", "<sample>.CEL.B Allele Freq")
  <calls>     birdseed.calls.txt: comment lines starting with '#', then
              probeset_id + one call column per sample (-1=NN, 0=AA,
              1=AB, 2=BB); NO chr/pos — joined onto the lrr_baf grid by
              probeset_id == Name

Outputs (under --out-dir, consumed by R/absolute_cn.R via the manifest):

  <cohort>/<sample>__logr.tsv.gz      chr, pos, value
  <cohort>/<sample>__baf.tsv.gz       chr, pos, value
  <cohort>/<sample>__het.tsv.gz       chr, pos, call (0/1/2; NN -> empty)
  manifest.tsv                        sample, cohort, logr, baf, het_calls

Sample columns containing 'hapmap' (case-insensitive) are reference
plates and are excluded — their role is already baked into the logR
(HAPMAP-reference normalisation).  Only autosomes (Chr 1-22) are kept:
ASCAT metrics and the WGD heuristics are autosomal, and integer chr
names keep the R-side chromosome filter exact.

Usage:
  python scripts/prepare_abs_inputs.py \
      --data-root GSEdata --out-dir GSEdata/raw_snp6 \
      --cohort cue:CUE/cue.lrr_baf.txt:CUE/birdseed.calls.txt \
      --cohort chb:CHB/CHB.lrr_baf.txt:CHB/birdseed.calls.txt
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import pandas as pd

AUTOSOMES = set(range(1, 23))


def split_spec(spec: str) -> tuple[str, Path, Path]:
    parts = spec.split(":")
    if len(parts) != 3:
        raise SystemExit(f"--cohort must be id:lrr_baf:calls, got {spec!r}")
    cid, lrr_baf, calls = parts
    return cid, Path(lrr_baf), Path(calls)


def sample_columns(columns: list[str]) -> dict[str, tuple[str, str]]:
    """Map tumour sample name -> (LRR column, BAF column)."""
    out: dict[str, tuple[str, str]] = {}
    lrr = {c.replace(".CEL.Log R Ratio", ""): c for c in columns
           if c.endswith(".CEL.Log R Ratio")}
    baf = {c.replace(".CEL.B Allele Freq", ""): c for c in columns
           if c.endswith(".CEL.B Allele Freq")}
    for name in lrr:
        if "hapmap" in name.lower():
            continue
        if name not in baf:
            raise SystemExit(f"LRR without BAF column for {name}")
        out[name] = (lrr[name], baf[name])
    if not out:
        raise SystemExit("no tumour sample columns found (all HapMap?)")
    return out


def process_cohort(cid: str, lrr_baf: Path, calls_path: Path,
                   out_dir: Path, chunksize: int) -> list[dict]:
    cohort_dir = out_dir / cid
    cohort_dir.mkdir(parents=True, exist_ok=True)

    # ---- birdseed calls: full load (small), indexed by probeset_id ----
    n_skip = 0
    with open(calls_path, errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                n_skip += 1
            else:
                break
    calls = pd.read_csv(calls_path, sep="\t", skiprows=n_skip)
    calls = calls.rename(columns={calls.columns[0]: "Name"})
    call_cols = {c.replace(".CEL", ""): c for c in calls.columns[1:]}
    calls = calls.set_index("Name")
    sys.stdout.write(f"[{cid}] calls: {len(calls)} markers x "
                     f"{len(call_cols)} samples\n")

    # ---- lrr_baf header -> tumour sample columns ----
    header = pd.read_csv(lrr_baf, sep="\t", nrows=0)
    samples = sample_columns(list(header.columns))
    missing = [s for s in samples if s not in call_cols]
    if missing:
        raise SystemExit(f"[{cid}] samples without birdseed calls: "
                         f"{missing[:5]}{'...' if len(missing) > 5 else ''}")
    sys.stdout.write(f"[{cid}] tumour samples: {len(samples)}\n")

    # ---- streaming split ----
    manifest_rows: list[dict] = []
    n_markers = 0
    n_dropped_xy = 0
    handles: dict[tuple[str, str], gzip.GzipFile] = {}
    first_done: set[tuple[str, str]] = set()

    def h(key: tuple[str, str]):
        if key not in handles:
            sample, kind = key
            handles[key] = gzip.open(
                cohort_dir / f"{sample}__{kind}.tsv.gz", "at",
                compresslevel=1)
        return handles[key]

    markers_h = gzip.open(cohort_dir / "markers.tsv.gz", "at",
                          compresslevel=1)
    markers_header_written = False

    reader = pd.read_csv(lrr_baf, sep="\t", chunksize=chunksize,
                         dtype={"Name": str, "Chr": str, "Position": "Int64"})
    for chunk in reader:
        n_all = len(chunk)
        chunk["ChrNum"] = pd.to_numeric(chunk["Chr"], errors="coerce")
        chunk = chunk[chunk["ChrNum"].isin(AUTOSOMES)]
        n_dropped_xy += n_all - len(chunk)
        if chunk.empty:
            continue
        grid = chunk[["ChrNum", "Position"]]
        joined = calls.reindex(chunk["Name"])
        # hygiene: birdseed calls exist for some CN/special probes whose
        # BAF column holds a placeholder (2.0), not an allele fraction —
        # their "calls" are meaningless and would poison the ASCAT het
        # set; a marker is only het-eligible when every sample's BAF is a
        # valid allele fraction (<= 1)
        baf_block = chunk[[baf_c for _, baf_c in samples.values()]]
        het_valid = (baf_block <= 1).all(axis=1).to_numpy()
        # probe-name grid (used by the R wrapper for GC-correction keys)
        mk = pd.DataFrame({"name": chunk["Name"].to_numpy(),
                           "chr": grid["ChrNum"].to_numpy(),
                           "pos": grid["Position"].to_numpy()})
        mk.to_csv(markers_h, sep="\t", index=False,
                  header=not markers_header_written)
        markers_header_written = True
        for sample, (lrr_col, baf_col) in samples.items():
            call_col = call_cols[sample]
            het_vals = joined[call_col].round().astype("Int64")
            het_vals = het_vals.where(het_vals.isin([0, 1, 2]) & het_valid)
            values = (
                ("logr", chunk[lrr_col].to_numpy(dtype="float64")),
                ("baf", chunk[baf_col].to_numpy(dtype="float64")),
                ("het", het_vals.to_numpy(dtype="object")),
            )
            for kind, series in values:
                df = pd.DataFrame({"chr": grid["ChrNum"].to_numpy(),
                                   "pos": grid["Position"].to_numpy(),
                                   "value": series})
                header_line = (sample, kind) not in first_done
                df.to_csv(h((sample, kind)), sep="\t", index=False,
                          header=header_line, na_rep="")
                if header_line:
                    first_done.add((sample, kind))
        n_markers += len(chunk)

    for handle in handles.values():
        handle.close()
    markers_h.close()

    for sample in samples:
        manifest_rows.append(dict(
            sample=sample, cohort=cid,
            logr=f"{cid}/{sample}__logr.tsv.gz",
            baf=f"{cid}/{sample}__baf.tsv.gz",
            het_calls=f"{cid}/{sample}__het.tsv.gz"))
    sys.stdout.write(f"[{cid}] autosomal markers written: {n_markers} "
                     f"(non-autosome rows dropped: {n_dropped_xy})\n")
    return manifest_rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-root", default="GSEdata", type=Path)
    ap.add_argument("--out-dir", default=None, type=Path,
                    help="default: <data-root>/raw_snp6")
    ap.add_argument("--cohort", action="append", required=True,
                    help="id:lrr_baf:calls (paths relative to data-root)")
    ap.add_argument("--chunksize", type=int, default=100_000)
    args = ap.parse_args()

    out_dir = args.out_dir or (args.data_root / "raw_snp6")
    out_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    for spec in args.cohort:
        cid, lrr_baf, calls = split_spec(spec)
        rows += process_cohort(cid, args.data_root / lrr_baf,
                               args.data_root / calls, out_dir,
                               args.chunksize)

    manifest = out_dir / "manifest.tsv"
    pd.DataFrame(rows).to_csv(manifest, sep="\t", index=False)
    sys.stdout.write(f"manifest: {manifest} ({len(rows)} samples)\n")


if __name__ == "__main__":
    main()
