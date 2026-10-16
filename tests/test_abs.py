"""Tests for module ABS (absolute copy-number comparisons; table mode +
opt-in R/ASCAT pipeline mode)."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

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
        "tool": "ascat", "auto_run": True,
        "min_purity": 0.30, "loh_frac_wgd": 0.30}
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


def test_abs_loh_analyses(pairwise_config, tmp_path):
    """LOH x GISTIC-Del cross table + WGD-stratified LOH rows (and the two
    figures) on the synthetic tables."""
    import yaml
    synth_root = pairwise_config.parent
    synth_absolute_tables(synth_root / "A", "a", seed=81,
                          truth=dict(wgd_frac=0.5, purity_mean=0.80,
                                     ploidy_mean=3.2),
                          sample_names=[f"A{i}" for i in range(1, 9)])
    synth_absolute_tables(synth_root / "B", "b", seed=82,
                          truth=dict(wgd_frac=0.5, purity_mean=0.60,
                                     ploidy_mean=2.0),
                          sample_names=[f"B{i}" for i in range(1, 9)])
    cfg = yaml.safe_load(open(pairwise_config))
    cfg.setdefault("analysis", {})["absolute_cn"] = {"tool": "ascat"}
    path = pairwise_config.parent / "config_abs_loh.yaml"
    path.write_text(yaml.safe_dump(cfg))

    from gistic_compare.cli import main
    out = tmp_path / "out_abs_loh"
    assert main(["--config", str(path), "--modules", "ABS",
                 "--outdir", str(out)]) == 0
    tables = out / "test_pair" / "tables"
    figures = out / "test_pair" / "figures"
    assert not (tables / "ABS_skipped.txt").exists()

    # the cross: one row per (cohort, Del peak); the four classes partition
    # the covered samples and the LOH-only class is split by mechanism
    cross = pd.read_csv(tables / "ABS_loh_del_cross.tsv", sep="\t")
    assert len(cross), "no Del peak x sample pairs — synth geometry broken"
    for col in ("del_loh", "del_only", "loh_only", "neither", "n_loh",
                "n_cnloh", "loh_only_cnloh", "loh_only_cnloss", "core",
                "gistic_freq", "n_samples"):
        assert col in cross.columns, col
    total = (cross["del_loh"] + cross["del_only"] + cross["loh_only"]
             + cross["neither"])
    assert (total == cross["n_covered"]).all()
    assert (cross["n_covered"] <= cross["n_samples"]).all()
    assert (cross["n_loh"] == cross["del_loh"] + cross["loh_only"]).all()
    assert (cross["loh_only_cnloh"] + cross["loh_only_cnloss"]
            <= cross["loh_only"]).all()
    assert (cross["core"] == "narrow").all()
    assert {"a", "b"} <= set(cross["cohort"])

    # the stratified LOH rows sit next to the pooled ones, and exist
    # exactly for the strata that carry >= 3 samples per cohort (the
    # module's own threshold)
    metrics = pd.read_csv(tables / "ABS_purity_ploidy.tsv", sep="\t")
    strat = metrics[metrics["comparison"].str.contains("only]", regex=False)]
    sc = pd.read_csv(tables / "ABS_cn_scores.tsv", sep="\t")
    expected = set()
    for sval, slab in ((1, "WGD+ only"), (0, "WGD- only")):
        sub = sc[sc["wgd_status"] == sval]
        if (sub["cohort"] == "a").sum() >= 3 and (sub["cohort"] == "b").sum() >= 3:
            expected.add(f"a_vs_b [{slab}]")
    assert set(strat["comparison"]) == expected
    assert {"loh_frac", "cnloh_frac"} <= set(strat["metric"])
    pooled = metrics[~metrics["comparison"].str.contains("only]",
                                                         regex=False)]
    assert {"loh_frac", "cnloh_frac"} <= set(pooled["metric"])

    summary = pd.read_csv(tables / "ABS_summary.tsv", sep="\t")
    for col in ("loh_only_frac_a", "loh_only_frac_b",
                "loh_only_cnloh_frac_a", "loh_only_cnloh_frac_b"):
        assert col in summary.columns, col

    for name in ("ABS_loh_del_cross.png", "ABS_loh_wgd_strata.png",
                 "ABS_metric_dists.png"):
        assert (figures / name).exists(), name


def test_abs_skip_without_tables(pairwise_config, tmp_path):
    from gistic_compare.cli import main
    out = tmp_path / "out_abs_skip"
    assert main(["--config", str(pairwise_config), "--modules", "ABS",
                 "--outdir", str(out), "--skip-figures"]) == 0
    note = out / "test_pair" / "tables" / "ABS_skipped.txt"
    assert note.exists()
    # auto_run defaults to false: the skip must point at the server-side
    # handoff instead of shelling out to R (Rscript/ASCAT are installed in
    # this environment, so only the default gate keeps it from running)
    text = note.read_text()
    assert "ascat_handoff" in text and "auto_run" in text


def test_abs_schema_mismatch_skips(pairwise_config, tmp_path):
    """A table missing a frozen column skips ABS with a pointer, it does
    not compute on mislabelled data."""
    import yaml
    synth_root = pairwise_config.parent
    synth_absolute_tables(synth_root / "A", "a", seed=61, n_samples=8,
                          sample_names=[f"A{i}" for i in range(1, 9)])
    synth_absolute_tables(synth_root / "B", "b", seed=62, n_samples=8,
                          sample_names=[f"B{i}" for i in range(1, 9)])
    bad = pd.read_csv(synth_root / "A" / "absolute_cn_a.tsv", sep="\t")
    bad.drop(columns=["loh_frac"]).to_csv(
        synth_root / "A" / "absolute_cn_a.tsv", sep="\t", index=False)

    cfg = yaml.safe_load(open(pairwise_config))
    cfg.setdefault("analysis", {})["absolute_cn"] = {"tool": "ascat"}
    path = pairwise_config.parent / "config_abs_bad.yaml"
    path.write_text(yaml.safe_dump(cfg))

    from gistic_compare.cli import main
    out = tmp_path / "out_abs_bad"
    assert main(["--config", str(path), "--modules", "ABS",
                 "--outdir", str(out), "--skip-figures"]) == 0
    note = out / "test_pair" / "tables" / "ABS_skipped.txt"
    assert "loh_frac" in note.read_text()


def test_abs_numeric_chrom_labels(pairwise_config, tmp_path):
    """'1' (server-side tools often emit plain numbers) is normalised to
    'chr1' so the cytoband/arm lookups keep matching."""
    import yaml
    synth_root = pairwise_config.parent
    synth_absolute_tables(synth_root / "A", "a", seed=71, n_samples=8,
                          sample_names=[f"A{i}" for i in range(1, 9)])
    synth_absolute_tables(synth_root / "B", "b", seed=72, n_samples=8,
                          sample_names=[f"B{i}" for i in range(1, 9)])
    for cid in ("a", "b"):
        p = synth_root / cid.upper() / f"absolute_seg_{cid}.tsv"
        seg = pd.read_csv(p, sep="\t")
        seg["chrom"] = seg["chrom"].str.replace("chr", "", regex=False)
        seg.to_csv(p, sep="\t", index=False)

    cfg = yaml.safe_load(open(pairwise_config))
    cfg.setdefault("analysis", {})["absolute_cn"] = {"tool": "ascat"}
    path = pairwise_config.parent / "config_abs_numchrom.yaml"
    path.write_text(yaml.safe_dump(cfg))

    from gistic_compare.cli import main
    out = tmp_path / "out_abs_numchrom"
    assert main(["--config", str(path), "--modules", "ABS",
                 "--outdir", str(out), "--skip-figures"]) == 0
    tables = out / "test_pair" / "tables"
    assert not (tables / "ABS_skipped.txt").exists()
    wgd = pd.read_csv(tables / "ABS_wgd.tsv", sep="\t")
    assert len(wgd) >= 1


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
    cfg.setdefault("analysis", {})["absolute_cn"] = {"tool": "ascat"}
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


# ---------------------------------------------------------------------------
# Phase 2: R/ASCAT pipeline mode (skipif no Rscript/ASCAT)
# ---------------------------------------------------------------------------

def _has_ascat() -> bool:
    import shutil
    import subprocess
    rs = shutil.which("Rscript")
    if rs is None:
        return False
    try:
        r = subprocess.run(
            [rs, "-e", 'cat(requireNamespace("ASCAT", quietly=TRUE))'],
            capture_output=True, text=True, timeout=120)
    except (OSError, subprocess.TimeoutExpired):
        return False
    return "TRUE" in r.stdout


@pytest.mark.skipif(not _has_ascat(), reason="Rscript/ASCAT unavailable")
def test_abs_r_pipeline_end_to_end(pairwise_config, tmp_path):
    """ASCAT on synthetic marker-level logR/BAF -> harmonized tables ->
    comparison suite (full pipeline mode, Python subprocess -> R)."""
    import yaml

    from synth_absolute import synth_ascat_inputs
    synth_root = pairwise_config.parent
    raw = synth_root / "raw_snp6"
    ma = synth_ascat_inputs(raw, "a", n_samples=4, wgd=True, purity=0.80,
                            seed=31, sample_names=[f"A{i}"
                                                   for i in range(1, 9)])
    mb = synth_ascat_inputs(raw, "b", n_samples=4, wgd=False, purity=0.60,
                            seed=32, sample_names=[f"B{i}"
                                                   for i in range(1, 9)])
    manifest = pd.concat([ma, mb], ignore_index=True)
    manifest.to_csv(raw / "manifest.tsv", sep="\t", index=False)

    cfg = yaml.safe_load(open(pairwise_config))
    cfg.setdefault("analysis", {})["absolute_cn"] = {
        "tool": "ascat", "auto_run": True, "raw_dir": "raw_snp6",
        "manifest": "manifest.tsv"}
    path = pairwise_config.parent / "config_abs_r.yaml"
    path.write_text(yaml.safe_dump(cfg))

    from gistic_compare.cli import main
    out = tmp_path / "out_abs_r"
    assert main(["--config", str(path), "--modules", "ABS",
                 "--outdir", str(out), "--skip-figures"]) == 0
    tables = out / "test_pair" / "tables"
    assert not (tables / "ABS_skipped.txt").exists()

    # the R wrapper emitted the harmonized tables into the cohort dirs
    assert (synth_root / "A" / "absolute_cn_a.tsv").exists()
    assert (synth_root / "B" / "absolute_seg_b.tsv").exists()

    cn = pd.read_csv(synth_root / "A" / "absolute_cn_a.tsv", sep="\t")
    for col in ("sample", "purity", "ploidy_abs", "wgd_status", "wgd_frac",
                "loh_frac", "cnloh_frac", "n_segments_allelic",
                "flag_low_purity", "wgd_flag_b", "wgd_disagreement"):
        assert col in cn.columns, col
    # truth: WGD+ cohort with strong duplication signature
    assert (cn["fit_ok"] == 1).all()
    assert (cn["wgd_frac"] > 0.3).all()
    assert (cn["ploidy_abs"] > 2.5).all()
    seg = pd.read_csv(synth_root / "B" / "absolute_seg_b.tsv", sep="\t")
    assert (seg["cn_minor"] <= seg["cn_total"]).all()

    wgd = pd.read_csv(tables / "ABS_wgd.tsv", sep="\t")
    r = wgd[wgd["comparison"] == "a_vs_b"].iloc[0]
    assert r["wgd_frac_a"] > r["wgd_frac_b"]
    assert r["odds_ratio"] > 1

    metrics = pd.read_csv(tables / "ABS_purity_ploidy.tsv", sep="\t")
    for metric in ("purity", "ploidy_abs"):
        rr = metrics[metrics["metric"] == metric].iloc[0]
        assert rr["cohort_higher"] == "a"


@pytest.mark.skipif(not _has_ascat(), reason="Rscript/ASCAT unavailable")
def test_abs_r_pipeline_tumor_only_het_calls(pairwise_config, tmp_path):
    """Tumour-only mode: het SNPs from the samples' own Birdseed calls
    (het_calls manifest column) — truth recovery on par with matched
    mode; germline_mode/het_source recorded in output + summary."""
    import yaml

    from synth_absolute import synth_ascat_inputs
    synth_root = pairwise_config.parent
    raw = synth_root / "raw_snp6"
    ma = synth_ascat_inputs(raw, "a", n_samples=4, wgd=True, purity=0.80,
                            seed=41, het_calls=True,
                            sample_names=[f"A{i}" for i in range(1, 9)])
    mb = synth_ascat_inputs(raw, "b", n_samples=4, wgd=False, purity=0.60,
                            seed=42, het_calls=True,
                            sample_names=[f"B{i}" for i in range(1, 9)])
    manifest = pd.concat([ma, mb], ignore_index=True)
    manifest = manifest.drop(columns=["germline_logr", "germline_baf"])
    manifest.to_csv(raw / "manifest_tonly.tsv", sep="\t", index=False)

    cfg = yaml.safe_load(open(pairwise_config))
    cfg.setdefault("analysis", {})["absolute_cn"] = {
        "tool": "ascat", "auto_run": True, "raw_dir": "raw_snp6",
        "manifest": "manifest_tonly.tsv"}
    path = pairwise_config.parent / "config_abs_tonly.yaml"
    path.write_text(yaml.safe_dump(cfg))

    from gistic_compare.cli import main
    out = tmp_path / "out_abs_tonly"
    assert main(["--config", str(path), "--modules", "ABS",
                 "--outdir", str(out), "--skip-figures"]) == 0
    tables = out / "test_pair" / "tables"
    assert not (tables / "ABS_skipped.txt").exists()

    cn_a = pd.read_csv(synth_root / "A" / "absolute_cn_a.tsv", sep="\t")
    assert (cn_a["germline_mode"] == "tumor_only").all()
    assert (cn_a["het_source"] == "birdseed_calls").all()
    assert (cn_a["fit_ok"] == 1).all()
    # truth recovery on par with the matched-mode run
    assert (cn_a["wgd_frac"] > 0.3).all()
    assert (cn_a["ploidy_abs"] > 2.5).all()
    cn_b = pd.read_csv(synth_root / "B" / "absolute_cn_b.tsv", sep="\t")
    # the 3-chromosome synthetic genome occasionally lets one sample
    # converge to the correlated (low purity, high ploidy) optimum —
    # an ASCAT identifiability degeneracy (plan §7.2); assert the rate
    assert (cn_b["wgd_status"] == 0).mean() >= 0.75
    assert (cn_b["ploidy_abs"] < 2.5).mean() >= 0.75

    summary = pd.read_csv(tables / "ABS_summary.tsv", sep="\t")
    assert summary["germline_mode_a"].iloc[0] == "tumor_only"
    assert summary["het_source_a"].iloc[0] == "birdseed_calls"


@pytest.mark.skipif(not _has_ascat(), reason="Rscript/ASCAT unavailable")
def test_abs_r_pipeline_inferred_fallback(pairwise_config, tmp_path):
    """No germline, no het_calls -> ASCAT's built-in inference; on the
    synthetic data the AffySNP6 genotype priors do not match, so the
    per-sample failure tolerance kicks in (fit_ok=0 rows) instead of a
    crash."""
    import yaml

    from synth_absolute import synth_ascat_inputs
    synth_root = pairwise_config.parent
    raw = synth_root / "raw_snp6"
    ma = synth_ascat_inputs(raw, "a", n_samples=4, wgd=True, purity=0.80,
                            seed=51, sample_names=[f"A{i}"
                                                   for i in range(1, 9)])
    mb = synth_ascat_inputs(raw, "b", n_samples=4, wgd=False, purity=0.60,
                            seed=52, sample_names=[f"B{i}"
                                                   for i in range(1, 9)])
    manifest = pd.concat([ma, mb], ignore_index=True)
    manifest = manifest.drop(columns=["germline_logr", "germline_baf"])
    manifest.to_csv(raw / "manifest_inf.tsv", sep="\t", index=False)
    cfg = yaml.safe_load(open(pairwise_config))
    cfg.setdefault("analysis", {})["absolute_cn"] = {
        "tool": "ascat", "auto_run": True, "raw_dir": "raw_snp6",
        "manifest": "manifest_inf.tsv"}
    path = pairwise_config.parent / "config_abs_inf.yaml"
    path.write_text(yaml.safe_dump(cfg))
    from gistic_compare.cli import main
    out = tmp_path / "out_abs_inf"
    assert main(["--config", str(path), "--modules", "ABS",
                 "--outdir", str(out), "--skip-figures"]) == 0
    tables = out / "test_pair" / "tables"
    cn_a = pd.read_csv(synth_root / "A" / "absolute_cn_a.tsv", sep="\t")
    assert (cn_a["germline_mode"] == "tumor_only").all()
    assert (cn_a["het_source"] == "inferred").all()
    # the failure-tolerance rows are written either way
    assert "fit_ok" in cn_a.columns
    assert not (tables / "ABS_skipped.txt").exists()
