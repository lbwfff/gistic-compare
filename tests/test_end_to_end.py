"""End-to-end runs of both workflows on the synthetic project."""

from __future__ import annotations

from pathlib import Path

from gistic_compare.cli import main

REPO = Path(__file__).resolve().parents[1]


def _run(cfg, argv):
    rc = main(["--config", str(cfg), *argv])
    assert rc == 0
    return cfg


def test_benchmark_all_modules(benchmark_config, tmp_path):
    out = tmp_path / "out_bench"
    _run(benchmark_config, ["--all", "--outdir", str(out)])
    tables = out / "test_bench" / "tables"
    figures = out / "test_bench" / "figures"
    assert tables.exists() and figures.exists()
    names = {p.name for p in tables.iterdir()}
    for mid in "ABCDEF":
        assert any(n.startswith(f"{mid}_") for n in names), mid
    assert "00_cohort_overview.tsv" in names
    # G was skipped (no reference panel configured) and must have produced
    # no rows beyond the module-skip path
    g_summary = tables / "G_paper_summary.tsv"
    assert not g_summary.exists()
    assert figures.glob("A_jaccard_interval.png")
    # run manifest records the skip reason
    import json
    mf = json.loads((tables.parent / "run_manifest.json").read_text())
    skips = [m for m in mf["modules"] if m["id"] == "G"]
    assert skips and skips[0]["status"] == "skipped" and "panel" in skips[0]["reason"]


def test_benchmark_with_reference_panel(benchmark_config, tmp_path):
    import yaml

    cfg = yaml.safe_load(open(benchmark_config))
    cfg["paper_benchmark"] = {
        "enabled": True,
        "panel": str(REPO / "presets" / "reference_panels" / "we_lung.yaml"),
    }
    path = benchmark_config.parent / "config_panel.yaml"
    path.write_text(yaml.safe_dump(cfg))
    out = tmp_path / "out_panel"
    _run(path, ["--modules", "G", "--outdir", str(out)])
    tables = out / "test_bench" / "tables"
    assert (tables / "G_paper_summary.tsv").exists()
    assert (tables / "G_paper_genome_pct.tsv").exists()


def test_pairwise_all_modules(pairwise_config, tmp_path):
    out = tmp_path / "out_pair"
    _run(pairwise_config, ["--all", "--outdir", str(out)])
    tables = out / "test_pair" / "tables"
    names = {p.name for p in tables.iterdir()}
    for mid in ("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "CF", "GS"):
        assert any(n.startswith(f"{mid}_") for n in names), mid
    # P4 found the synthetic chr1 5-10 Mb differential region
    import pandas as pd
    diff = pd.read_csv(tables / "P4_regional_diff.tsv", sep="\t")
    hit = diff[(diff["chrom"] == "chr1") & (diff["start"] == 5_000_000)]
    assert not hit.empty
    assert (hit["freq_amp_a"] > hit["freq_amp_b"]).all()
    figures = out / "test_pair" / "figures"
    assert figures.glob("P1_jaccard.png")
    assert figures.glob("P4_genome_diff_fisher.png")
    assert figures.glob("GS_band_genome_amp.png")
    # MDE / power table: one row per comparison x direction x baseline
    power = pd.read_csv(tables / "P4_power.tsv", sep="\t")
    assert len(power) == 1 * 2 * 5          # comps x {Amp,Del} x baselines
    assert {"mde_up", "mde_down", "mde", "baseline_freq"} <= set(power)
    assert (power["mde"] > 0).all() and (power["mde"] < 1).all()
    summary = pd.read_csv(tables / "P4_summary.tsv", sep="\t")
    assert {"max_abs_diff_amp", "n_windows_diff_gt_mde_amp"} <= set(summary)
    # P5 shares the machinery
    assert (tables / "P5_power.tsv").exists()
    # OncoPrint: long call table + figure
    onco_calls = pd.read_csv(tables / "P1_onco_calls.tsv", sep="\t")
    assert {"peak_key", "sample", "cohort", "call"} <= set(onco_calls)
    assert figures.glob("P1_oncoprint*.png")
    # run manifest is written next to the tables
    import json
    mf = json.loads((tables.parent / "run_manifest.json").read_text())
    assert mf["run_id"] == "test_pair"
    assert all(m["status"] == "ok" for m in mf["modules"])


def test_pairwise_loci_preset_produces_paper_tables(pairwise_config, tmp_path):
    import yaml

    preset = pairwise_config.parent / "loci_test.yaml"
    preset.write_text(yaml.safe_dump(dict(
        id="test",
        loci=[dict(locus="test_region", chrom="chr1", start=5000000,
                   end=10000000, type="Amp", expected_cohort="A")],
        bands=[dict(band="1p36.33", star=False)])))
    cfg = yaml.safe_load(open(pairwise_config))
    cfg["analysis"]["regional"]["loci_set"] = str(preset)
    path = pairwise_config.parent / "config_loci.yaml"
    path.write_text(yaml.safe_dump(cfg))
    out = tmp_path / "out_loci"
    _run(path, ["--modules", "P4", "--outdir", str(out)])
    tables = out / "test_pair" / "tables"
    import pandas as pd
    loci = pd.read_csv(tables / "P4_paper_loci.tsv", sep="\t") \
        if (tables / "P4_paper_loci.tsv").exists() else None
    assert loci is not None and len(loci) == 1
    assert loci.iloc[0]["locus"] == "test_region"


def test_unknown_module_errors(pairwise_config):
    rc = main(["--config", str(pairwise_config), "--modules", "Z"])
    assert rc == 2


def test_skip_reason_g_without_panel(benchmark_config, capsys):
    _run(benchmark_config, ["--all", "--skip-figures"])
    captured = capsys.readouterr()
    assert "[skip] G:" in captured.out
