"""Merged modules: P4's marker-matrix layer (ex-module P5) and the shared
call-concordance implementation behind benchmark C and pairwise P2."""

from __future__ import annotations

import pandas as pd


def _tables(out: object, run: str = "test_pair"):
    return out / run / "tables"


def _run(cfg, argv):
    from gistic_compare.cli import main
    assert main(["--config", str(cfg), *argv]) == 0


def test_p2_frequency_layer_on_unpaired(pairwise_config, tmp_path):
    """The synth pair is unpaired: P2 writes the frequency layer (defined for
    every comparison) and leaves the paired-only layers empty."""
    _run(pairwise_config, ["--modules", "P2", "--outdir", str(tmp_path)])
    tables = _tables(tmp_path)
    for name in ("frequency", "summary", "peak_agreement", "sample_agreement",
                 "mcnemar"):
        assert (tables / f"P2_{name}.tsv").exists(), name
    freq = pd.read_csv(tables / "P2_frequency.tsv", sep="\t")
    assert not freq.empty
    assert {"comparison", "kind", "ptype", "a_freq", "b_freq",
            "delta_freq", "a_q", "b_q"} <= set(freq.columns)
    # a paired-only layer is all header
    assert pd.read_csv(tables / "P2_peak_agreement.tsv", sep="\t").empty
    assert pd.read_csv(tables / "P2_summary.tsv", sep="\t").empty


def test_c_and_p2_share_the_summary_shape(benchmark_config, pairwise_config,
                                          tmp_path):
    _run(benchmark_config, ["--modules", "C", "--outdir", str(tmp_path / "bc")])
    _run(pairwise_config, ["--modules", "P2", "--outdir", str(tmp_path / "pw")])
    bench = _tables(tmp_path / "bc", "test_bench")
    pair = _tables(tmp_path / "pw")
    c_sum = pd.read_csv(bench / "C_summary.tsv", sep="\t")
    p2_sum = pd.read_csv(pair / "P2_summary.tsv", sep="\t")
    # the two forms share one implementation and one table shape
    assert list(c_sum.columns) == list(p2_sum.columns)
    assert {"comparison", "kind", "ptype", "n_matched_peaks",
            "n_paired_samples"} <= set(c_sum.columns)


def test_p4_marker_layer_tables_match_the_seg_layer(pairwise_config, tmp_path):
    _run(pairwise_config, ["--modules", "P4", "--outdir", str(tmp_path)])
    tables = _tables(tmp_path)
    for name in ("regional_diff", "cytoband_diff", "power"):
        p4 = pd.read_csv(tables / f"P4_{name}.tsv", sep="\t")
        p5 = pd.read_csv(tables / f"P5_{name}.tsv", sep="\t")
        assert list(p4.columns) == list(p5.columns), name


def test_p4_marker_layer_can_be_disabled(pairwise_config, tmp_path):
    import yaml
    cfg = yaml.safe_load(open(pairwise_config))
    cfg["analysis"]["regional"]["marker_matrix"] = False
    path = pairwise_config.parent / "config_no_marker.yaml"
    path.write_text(yaml.safe_dump(cfg))
    _run(path, ["--modules", "P4", "--outdir", str(tmp_path)])
    assert not (_tables(tmp_path) / "P5_regional_diff.tsv").exists()
    assert (_tables(tmp_path) / "P4_regional_diff.tsv").exists()


def test_p4_marker_layer_skipped_when_a_cohort_lacks_the_matrix(
        pairwise_config, tmp_path):
    _run(pairwise_config, ["--modules", "P4", "--outdir", str(tmp_path)])
    tables = _tables(tmp_path)
    # drop only 'a's matrix and re-run into a fresh directory
    import yaml
    cfg = yaml.safe_load(open(pairwise_config))
    (pairwise_config.parent / "A" / "ca.D.cap1.5.mat").unlink()
    path = pairwise_config.parent / "config_run2.yaml"
    path.write_text(yaml.safe_dump(cfg))
    out2 = tmp_path / "out2"
    _run(path, ["--modules", "P4", "--outdir", str(out2)])
    assert (_tables(out2) / "P4_regional_diff.tsv").exists()
    assert not (_tables(out2) / "P5_regional_diff.tsv").exists()
    del tables
