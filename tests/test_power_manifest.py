"""Tests for the MDE / power table (P4/P5) and run manifest."""

from __future__ import annotations

import json

import numpy as np
import pytest

from gistic_compare.pairwise_regional import _mde_two_props, _mde_window
from gistic_compare.run_manifest import _sanitize, build_manifest


# ---------------------------------------------------------------- MDE maths
def test_mde_decreases_with_sample_size():
    prev = np.inf
    for n in (10, 25, 50, 100, 500, 5000):
        _, _, m = _mde_two_props(n, n, 0.10)
        assert 0 < m < prev
        prev = m


def test_mde_symmetric_at_half():
    up, dn, m = _mde_two_props(100, 100, 0.5)
    assert up == pytest.approx(dn, abs=1e-9)
    assert m == pytest.approx(up)


def test_mde_positive_and_bounded():
    for base in (0.02, 0.1, 0.3, 0.7, 0.95):
        up, dn, m = _mde_two_props(80, 120, base)
        assert up >= 0 and dn >= 0 and m > 0
        assert base + up <= 1.0 + 1e-9
        assert base - dn >= -1e-9


def test_mde_matches_wald_power_for_moderate_n():
    """Analytic arcsine MDE should give ~80% power for a Wald z-test."""
    from math import asin, sin, sqrt

    from scipy.stats import norm
    na = nb = 500
    base = 0.1
    _, _, m = _mde_two_props(na, nb, base)
    za = norm.ppf(0.975)
    p_alt = sin(asin(sqrt(base)) + (za + norm.ppf(0.8)) * 0.5
                * (1 / na + 1 / nb) ** 0.5) ** 2
    # power of the two-proportion z test at p_alt vs base
    se0 = np.sqrt((p_alt + base) / 2 * (1 - (p_alt + base) / 2) * (2 / na))
    se1 = np.sqrt(p_alt * (1 - p_alt) / na + base * (1 - base) / nb)
    power = 1 - norm.cdf((za * se0 - (p_alt - base)) / se1)
    assert 0.75 < power < 0.85
    assert p_alt - base == pytest.approx(m, abs=1e-6)


def test_mde_window_vectorised():
    na = np.full(4, 50.0)
    nb = np.full(4, 50.0)
    p = np.array([0.05, 0.2, 0.5, 0.9])
    m = _mde_window(na, nb, p, 0.05, 0.8)
    assert m.shape == (4,)
    assert np.all(np.isfinite(m)) and np.all(m > 0)
    # mid-range baselines need the largest absolute shift; extremes need less
    assert m[2] == m.max()


# ------------------------------------------------------------- run manifest
def test_sanitize():
    from pathlib import Path as P
    out = _sanitize({"a": P("/tmp/x"), "b": {1, 2}, "c": (1, 2),
                     "d": np.float64(1.5), "e": [P("/y")]})
    assert out["a"] == "/tmp/x"
    assert out["b"] == ["1", "2"]
    assert out["c"] == [1, 2]
    assert out["d"] == 1.5
    assert out["e"] == ["/y"]


def test_manifest_content(synth_root, pairwise_config, tmp_path):
    from gistic_compare.config import load_config

    cfg = load_config(pairwise_config)
    mods = [dict(id="P1", label="peak set", status="ok", seconds=1.0),
            dict(id="G", label="panel", status="skipped", reason="no panel")]
    manifest = build_manifest(cfg, mods, total_seconds=9.5)
    assert manifest["run_id"] == "test_pair"
    assert manifest["workflow"] == "pairwise"
    assert manifest["tool_version"]
    assert manifest["dependencies"]["pandas"]
    assert {m["id"] for m in manifest["modules"]} == {"P1", "G"}
    # full config snapshot (minus output) is embedded and JSON-safe
    assert manifest["config"]["comparison_mode"] == "experiment_pairs"
    json.dumps(manifest)  # must not raise

    from gistic_compare.run_manifest import write_manifest
    path = write_manifest(cfg, mods, total_seconds=9.5)
    assert json.loads(path.read_text())["run_id"] == "test_pair"


def test_manifest_written_by_cli(pairwise_config, tmp_path):
    from gistic_compare.cli import main

    out = tmp_path / "out_m"
    assert main(["--config", str(pairwise_config), "--modules", "P4",
                 "--outdir", str(out), "--skip-figures"]) == 0
    mf = json.loads((out / "test_pair" / "run_manifest.json").read_text())
    assert mf["modules"][0]["status"] == "ok"
    assert mf["modules"][0]["seconds"] >= 0
    assert mf["cohorts"]["a"]["n_samples"] == 8
    assert mf["config"]["run_id"] == "test_pair"
