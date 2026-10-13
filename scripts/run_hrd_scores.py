"""Compute scarHRD genomic-scar scores (HRD-LOH / TAI / LST / HRD-sum)
for every sample of the ABS pipeline, via the vendored reference
implementation (R/scarHRD/), and write ABS_hrd_scores.tsv.

Usage: python scripts/run_hrd_scores.py [--config config.local.yaml]
"""
from __future__ import annotations

import argparse
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

from gistic_compare.config import load_config, out_table
from gistic_compare.pairwise_absolute import _abs_cfg, _regional_cfg, _wgd_class_columns

_R_SCRIPT = Path(__file__).resolve().parent.parent / "R" / "scarHRD" / "hrd_scores.R"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default="config.local.yaml")
    args = ap.parse_args()
    cfg = load_config(args.config)
    acfg = _abs_cfg(cfg, _regional_cfg(cfg))
    data_root = Path(cfg["data_root"])

    cn_tables, seg_tables = {}, {}
    for cid, c in (cfg.get("cohorts") or {}).items():
        cn_tables[cid] = pd.read_csv(
            data_root / c["path"] / f"absolute_cn_{cid}.tsv", sep="\t")
        seg_tables[cid] = pd.read_csv(
            data_root / c["path"] / f"absolute_seg_{cid}.tsv", sep="\t")
    cn_tables = _wgd_class_columns(seg_tables, cn_tables, acfg)
    cn = pd.concat(cn_tables.values(), ignore_index=True)
    coh = pd.concat([t.assign(cohort=cid)[["sample", "cohort"]]
                     for cid, t in cn_tables.items()], ignore_index=True)

    ploidy = cn.set_index("sample")["ploidy_abs"]
    inp = seg_all = pd.concat(seg_tables.values(), ignore_index=True)
    inp = pd.DataFrame({
        "sample": seg_all["sample"], "chrom": seg_all["chrom"],
        "start": seg_all["start"], "end": seg_all["end"],
        "totalCN": seg_all["cn_total"],
        "A": seg_all[["cn_major", "cn_minor"]].max(axis=1),
        "B": seg_all[["cn_major", "cn_minor"]].min(axis=1),
        "ploidy": seg_all["sample"].map(ploidy)})
    if inp["ploidy"].isna().any():
        raise SystemExit("ploidy missing for some samples — cn tables out of "
                         "date? re-run the ABS module first")

    with tempfile.TemporaryDirectory() as td:
        inp.to_csv(Path(td) / "input_seg.txt", sep="\t", index=False)
        subprocess.run(["Rscript", str(_R_SCRIPT),
                        str(Path(td) / "input_seg.txt"), td, "grch37"],
                       check=True)
        r = pd.read_csv(Path(td) / "all_HRDresults.tsv", sep="\t",
                        index_col=0)

    r.index.name = "sample"
    r.columns = ["HRD_LOH", "TAI", "LST", "HRD_sum"]
    out = cn.merge(r.reset_index(), on="sample", how="left")
    out = out.merge(coh, on="sample", how="left")
    out["HRD_adj"] = out["HRD_sum"] - 15.5 * out["ploidy_abs"]
    dst = out_table(cfg, "ABS_hrd_scores.tsv")
    out.to_csv(dst, sep="\t", index=False)
    print(f"[HRD] wrote {dst} ({len(out)} samples)")


if __name__ == "__main__":
    main()
