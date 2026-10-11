"""Tests for module ABS (absolute copy-number comparisons, table mode)."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from synth_absolute import synth_absolute_tables  # noqa: E402

# ---------------------------------------------------------------------------
# synth generator
# ---------------------------------------------------------------------------

def test_synth_absolute_tables(tmp_path):
    cn_a, seg_a = synth_absolute_tables(tmp_path, "a", n_samples=30,
                                        seed=1)
    assert {"sample", "purity", "ploidy_abs", "wgd_status", "wgd_frac",
            "loh_frac", "cnloh_frac", "flag_low_purity"} <= set(cn_a.columns)
    assert {"sample", "chrom", "start", "end", "cn_total", "cn_major",
            "cn_minor"} <= set(seg_a.columns)
    # duplication signature: WGD+ samples carry mostly minor>=1/total>=3
    dup_frac = cn_a["wgd_frac"]
    assert (dup_frac >= 0).all() and (dup_frac <= 1).all()
    assert (seg_a["cn_minor"] <= seg_a["cn_total"]).all()
    assert (seg_a["cn_total"] >= 0).all()


# ---------------------------------------------------------------------------
# comparison statistics
# ---------------------------------------------------------------------------

def test_abs_end_to_end(pairwise_config, tmp_path):
    """Write synth absolute tables into the cohort dirs, run ABS."""
    import yaml
    # cohort dirs of the synth project (path: A / B per conftest); sample
    # names must match the synth GISTIC bundle (A1..A8 / B1..B8)
    synth_root = pairwise_config.parent
    strong_a = dict(wgd_frac=0.60, purity_mean=0.80, ploidy_mean=3.2)
    strong_b = dict(wgd_frac=0.15, purity_mean=0.50, ploidy_mean=2.0)
    synth_absolute_tables(synth_root / "A", "a", seed=11, truth=strong_a,
                          sample_names=[f"A{i}" for i in range(1, 9)])
    synth_absolute_tables(synth_root / "B", "b", seed=12, truth=strong_b,
                          sample_names=[f"B{i}" for i in range(1, 9)])

    cfg = yaml.safe_load(open(pairwise_config))
    cfg.setdefault("analysis", {})["absolute_cn"] = {
        "tool": "facets", "min_purity": 0.30, "loh_frac_wgd": 0.30}
    path = pairwise_config.parent / "config_abs.yaml"
    path.write_text(yaml.safe_dump(cfg))

    from gistic_compare.cli import main
    out = tmp_path / "out_abs"
    assert main(["--config", str(path), "--modules", "ABS",
                 "--outdir", str(out), "--skip-figures"]) == 0
    tables = out / "test_pair" / "tables"
    assert not (tables / "ABS_skipped.txt").exists()

    wgd = pd.read_csv(tables / "ABS_wgd.tsv", sep="\t")
    r = wgd[wgd["comparison"] == "a_vs_b"].iloc[0]
    # synth truth: A WGD+ ~60%, B ~15% -> A higher, OR > 1, significant
    assert r["wgd_frac_a"] > r["wgd_frac_b"]
    assert r["odds_ratio"] > 3
    assert r["p_fisher"] < 0.05

    metrics = pd.read_csv(tables / "ABS_purity_ploidy.tsv", sep="\t")
    for metric, higher in (("purity", "a"), ("ploidy_abs", "a"),
                           ("loh_frac", "a")):
        rr = metrics[metrics["metric"] == metric].iloc[0]
        assert rr["cohort_higher"] == higher
        # n = 8 per cohort: loh_frac has limited power
        assert rr["p_mw"] < (0.01 if metric != "loh_frac" else 0.10)

    cn = pd.read_csv(tables / "ABS_cn_tests.tsv", sep="\t")
    st = cn[cn["metric"] == "S_chrom_cn_star"].iloc[0]
    # WGD+ cohort a: chromosome-level duplication events
    assert st["cohort_higher"] == "a"

    peaks = pd.read_csv(tables / "ABS_purity_adjusted_peaks.tsv", sep="\t")
    assert {"crude_or", "adj_or_purity", "adj_q_purity"} <= set(peaks.columns)
    assert (peaks["crude_or"] > 0).all()


def test_abs_skip_without_tables(pairwise_config, tmp_path):
    from gistic_compare.cli import main
    out = tmp_path / "out_abs_skip"
    assert main(["--config", str(pairwise_config), "--modules", "ABS",
                 "--outdir", str(out), "--skip-figures"]) == 0
    assert (out / "test_pair" / "tables" / "ABS_skipped.txt").exists()


def test_abs_low_purity_sensitivity(pairwise_config, tmp_path):
    import yaml
    synth_root = pairwise_config.parent
    cn_a, _ = synth_absolute_tables(
        synth_root / "A", "a", seed=21,
        truth=dict(wgd_frac=0.60, purity_mean=0.55, purity_sd=0.22,
                   ploidy_mean=3.2),
        sample_names=[f"A{i}" for i in range(1, 9)])
    synth_absolute_tables(synth_root / "B", "b", seed=22,
                          sample_names=[f"B{i}" for i in range(1, 9)])
    cfg = yaml.safe_load(open(pairwise_config))
    cfg.setdefault("analysis", {})["absolute_cn"] = {"tool": "facets"}
    path = pairwise_config.parent / "config_abs2.yaml"
    path.write_text(yaml.safe_dump(cfg))
    from gistic_compare.cli import main
    out = tmp_path / "out_abs2"
    assert main(["--config", str(path), "--modules", "ABS",
                 "--outdir", str(out), "--skip-figures"]) == 0
    tables = out / "test_pair" / "tables"
    wgd = pd.read_csv(tables / "ABS_wgd.tsv", sep="\t")
    assert (wgd["comparison"] == "a_vs_b [high-purity only]").any()
    summary = pd.read_csv(tables / "ABS_summary.tsv", sep="\t")
    assert summary["n_low_purity"].iloc[0] >= 0
    assert not cn_a.empty
