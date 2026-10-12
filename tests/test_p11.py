"""Tests for module P11 (genome-wide focal q-value landscape)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from gistic_compare import pairwise_focal


def test_parse_wide_peak():
    assert pairwise_focal._parse_wide_peak(
        "chr11:69334308-69344762") == ("chr11", 69334308.0, 69344762.0)
    chrom, s, e = pairwise_focal._parse_wide_peak("nonsense")
    assert chrom == "" and np.isnan(s) and np.isnan(e)


def test_pick_gene_skips_mirna_and_loc():
    assert pairwise_focal._pick_gene(["[CCND1]"]) == "CCND1"
    assert pairwise_focal._pick_gene(
        ["hsa-mir-4277", "LOC63930", "TPD52"]) == "TPD52"
    # all-skippable -> falls back to the top gene (cleaned)
    assert pairwise_focal._pick_gene(["[hsa-mir-21]"]) == "hsa-mir-21"


def test_peak_records_one_row_per_peak():
    genes = pd.DataFrame([
        dict(cohort="a", type="Amp", col_index=0, cytoband="8q24.21",
             q=1e-6, wide_peak="chr8:128734245-128800300", gene="MYC"),
        dict(cohort="a", type="Amp", col_index=0, cytoband="8q24.21",
             q=1e-6, wide_peak="chr8:128734245-128800300", gene="FAM84B"),
        # a peak whose only genes are non-labelable -> keeps the top gene
        dict(cohort="a", type="Del", col_index=1, cytoband="9p21.3",
             q=1e-9, wide_peak="chr9:21865498-22448737", gene="hsa-mir-491"),
        # off-autosome peaks are dropped
        dict(cohort="a", type="Amp", col_index=2, cytoband="Xq1",
             q=1e-3, wide_peak="chrX:1000-2000", gene="FOO"),
    ])
    out = pairwise_focal._peak_records(genes, "a")
    assert list(out.columns) == ["cohort", "type", "cytoband", "chrom", "start",
                                 "end", "q", "neg_log10_q", "gene"]
    assert len(out) == 2                       # one per peak, off-autosome gone
    myc = out[out["gene"] == "MYC"].iloc[0]
    assert myc["chrom"] == "chr8" and myc["neg_log10_q"] == pytest.approx(6.0)
    # bracketed/miRNA top gene is kept when nothing else is available
    assert out[out["type"] == "Del"].iloc[0]["gene"] == "hsa-mir-491"


def test_p11_end_to_end(pairwise_config, tmp_path):
    from gistic_compare.cli import main
    out = tmp_path / "out_p11"
    assert main(["--config", str(pairwise_config), "--modules", "P11",
                 "--outdir", str(out)]) == 0
    base = out / "test_pair"
    peaks = pd.read_csv(base / "tables" / "P11_focal_peaks.tsv", sep="\t")
    assert {"cohort", "type", "chrom", "start", "end", "q", "neg_log10_q",
            "gene"} <= set(peaks.columns)
    assert set(peaks["type"]) <= {"Amp", "Del"}
    assert peaks["neg_log10_q"].dropna().ge(0).all()
    for name in ("P11_focal_landscape.png", "P11_focal_landscape.pdf"):
        assert (base / "figures" / name).exists()
