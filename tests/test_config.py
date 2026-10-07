"""Config / topology / presets tests."""

from __future__ import annotations

import pytest

from gistic_compare.config import build_comparisons, load_config
from gistic_compare.genomes import genome_build


def test_control_vs_all(synth_root, benchmark_config):
    cfg = load_config(benchmark_config)
    cmps = {(c["a"], c["b"]) for c in cfg["comparisons"]}
    assert cmps == {("c", "a"), ("c", "b")}
    assert all(c["kind"] == "unpaired" for c in cfg["comparisons"])


def test_experiment_pairs(synth_root, pairwise_config):
    cfg = load_config(pairwise_config)
    assert [(c["a"], c["b"]) for c in cfg["comparisons"]] == [("a", "b")]
    # kinds start unknown, filled from sample overlap after load
    assert cfg["comparisons"][0]["kind"] == "unknown"


def test_pairwise_with_control_rejected(pairwise_config):
    import yaml
    cfg = yaml.safe_load(open(pairwise_config))
    cfg["control"] = "a"
    with pytest.raises(ValueError, match="must not set 'control'"):
        build_comparisons({**cfg, "workflow": "pairwise"})


def test_benchmark_requires_control(synth_root, benchmark_config):
    import yaml
    cfg = yaml.safe_load(open(benchmark_config))
    cfg.pop("control")
    with pytest.raises(ValueError, match="requires 'control'"):
        build_comparisons({**cfg, "workflow": "benchmark"})


def test_unknown_genome_rejected():
    with pytest.raises(ValueError, match="hg19"):
        genome_build("hg38")


def test_output_dir_resolution(synth_root, benchmark_config, tmp_path):
    out = tmp_path / "custom_out"
    cfg = load_config(benchmark_config, outdir=out)
    assert cfg["output"]["_tables"] == (out / "test_bench" / "tables")
    assert cfg["output"]["_figures"] == (out / "test_bench" / "figures")
    assert cfg["output"]["_tables"].exists()


def test_output_never_written_inside_package(synth_root, benchmark_config,
                                             tmp_path):
    from pathlib import Path

    import gistic_compare
    pkg_dir = Path(gistic_compare.__file__).parent
    cfg = load_config(benchmark_config)
    assert pkg_dir not in cfg["output"]["_tables"].parents
    assert pkg_dir not in cfg["output"]["_figures"].parents


def test_loci_preset_loading(synth_root, pairwise_config, tmp_path):
    import yaml
    preset = tmp_path / "loci.yaml"
    preset.write_text(yaml.safe_dump(dict(
        id="test",
        loci=[dict(locus="test_locus", chrom="chr1", start=5000000,
                   end=10000000, type="Amp", expected_cohort="A")],
        bands=[dict(band="1q21.2", star=True)])))
    cfg = synth_root / "config_preset.yaml"
    base = yaml.safe_load(open(pairwise_config))
    base["analysis"]["regional"]["loci_set"] = str(preset)
    cfg.write_text(yaml.safe_dump(base))
    loaded = load_config(cfg)
    assert loaded["_preset_loci"]["loci"][0][0] == "test_locus"
    assert loaded["_preset_loci"]["bands"] == [("1q21.2", True)]


def test_missing_cohort_path(tmp_path):
    cfg = tmp_path / "bad.yaml"
    cfg.write_text(
        "workflow: pairwise\ndata_root: /nonexistent\ncohorts:\n"
        "  a: {label: A, path: nope, prefix: p}\n")
    with pytest.raises(FileNotFoundError, match="cohort path"):
        load_config(cfg)
