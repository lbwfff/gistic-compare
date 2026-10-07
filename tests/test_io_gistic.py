"""Parser tests for gistic_compare.io_gistic."""

from __future__ import annotations

import pytest

from gistic_compare import io_gistic


def test_parse_coord():
    chrom, s, e = io_gistic.parse_coord(
        "chr1:13884087-15253589(probes 5947:7043)")
    assert (chrom, s, e) == ("chr1", 13884087, 15253589)


def test_parse_coord_invalid():
    with pytest.raises(ValueError):
        io_gistic.parse_coord("not a coordinate")


def test_normalize_chrom():
    assert io_gistic.normalize_chrom("1") == "chr1"
    assert io_gistic.normalize_chrom("chr7") == "chr7"
    assert io_gistic.normalize_chrom("23") == "chr23"
    assert io_gistic.normalize_chrom("M") is None
    assert io_gistic.normalize_chrom("") is None


def test_load_all_lesions_skips_cn_rows(synth_root, benchmark_config):
    """Each peak is written twice (summary + CN values); only summary counts."""
    peaks, calls = io_gistic.load_all_lesions(
        synth_root / "A" / "ca.all_lesions.conf_90.txt", "a")
    assert len(peaks) == 5  # 3 shared + 2 cohort-specific
    assert set(peaks["type"]) == {"Amp", "Del"}
    assert (peaks["type"] == "Amp").sum() == 3
    # calls are long-format with 0/1/2 values
    assert set(calls["call"].unique()) <= {0, 1, 2}
    assert calls["sample"].nunique() == 8
    assert "q" in peaks.columns and "freq_altered" in peaks.columns


def test_load_cohort_bundle(synth_root, benchmark_config):
    from gistic_compare.config import load_config
    cfg = load_config(benchmark_config)
    d = io_gistic.load_cohort(cfg, "a")
    for key in ("peaks", "calls", "genes", "scores", "seg", "seg_counts",
                "bed"):
        assert key in d and d[key] is not None
    assert not d["scores"].empty
    assert {"cohort", "type", "chrom", "start", "end", "g_score"} <= \
        set(d["scores"].columns)


def test_load_all_and_sample_counts(synth_root, benchmark_config):
    from gistic_compare.config import load_config
    cfg = load_config(benchmark_config)
    data = io_gistic.load_all(cfg)
    assert set(data) == {"a", "b", "c"}
    assert cfg["_n_samples"] == {"a": 8, "b": 8, "c": 8}
