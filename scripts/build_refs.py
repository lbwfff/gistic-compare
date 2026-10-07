#!/usr/bin/env python3
"""Export small reference tables from GISTIC's hg19.mat assets.

GISTIC 2.0 ships `hg19.mat` containing the cytoband table (`cyto`) and the
refgene annotation (`rg`).  gistic-compare bundles slim TSV exports of both so
users do not have to download the full GISTIC assets.  Re-run this script to
regenerate them from a newer asset bundle.

Usage::

    python scripts/build_refs.py --mat /path/to/hg19.mat [--out-dir <dir>]

The default --mat locations match the standard GISTIC binary distribution
layout; point it at whatever file you have.  Provenance: Broad Institute
GISTIC 2.0 reference bundle (hg19, affy U133plus2/gistic reference files).
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def load_mat(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    import numpy as np
    import scipy.io as sio
    m = sio.loadmat(str(path), squeeze_me=True, struct_as_record=False)

    def records(key: str) -> list:
        arr = m.get(key)
        if arr is None:
            return []
        return list(np.atleast_1d(arr))

    cyto_rows = []
    for c in records("cyto"):
        d = {f: getattr(c, f) for f in c._fieldnames}
        chrom = str(d["chr"])
        if not chrom.startswith("chr"):
            chrom = "chr" + chrom
        cyto_rows.append(dict(chrom=chrom, start=int(d["start"]),
                              end=int(d["end"]), band=str(d["name"])))
    cyto = pd.DataFrame(cyto_rows)

    gene_rows = []
    for e in records("rg"):
        d = {f: getattr(e, f) for f in e._fieldnames}
        sym = str(d.get("symb", "")).strip()
        if not sym:
            continue
        gene_rows.append(dict(chrom=str(d.get("chr", "")),
                              start=int(d["start"]), end=int(d["end"]),
                              gene=sym))
    refgene = pd.DataFrame(gene_rows)
    return cyto, refgene


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mat", required=True, help="GISTIC hg19.mat asset file")
    ap.add_argument("--out-dir", default=None,
                    help="output dir (default: src/gistic_compare/data/refs)")
    args = ap.parse_args()

    out_dir = (Path(args.out_dir) if args.out_dir else
               Path(__file__).resolve().parents[1] / "src" /
               "gistic_compare" / "data" / "refs")
    out_dir.mkdir(parents=True, exist_ok=True)

    cyto, refgene = load_mat(Path(args.mat))
    cyto.to_csv(out_dir / "hg19_cytobands.tsv", sep="\t", index=False)
    refgene.to_csv(out_dir / "hg19_refgene.tsv", sep="\t", index=False)
    print(f"cytobands: {len(cyto)} rows; refgene: {len(refgene)} rows")
    print(f"written to {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
