"""Unit tests for the P1 union-peak OncoPrint data builder."""

from __future__ import annotations

import numpy as np

from gistic_compare.pairwise_peak import build_onco_data


def test_build_onco_data_shape(synth_root, pairwise_config):
    from gistic_compare.config import load_config

    cfg = load_config(pairwise_config)
    from gistic_compare import io_gistic
    data = io_gistic.load_all(cfg)
    from gistic_compare.config import assign_comparison_kinds
    assign_comparison_kinds(cfg, data)

    onco = build_onco_data(cfg, data, cfg["analysis"]["q_thresholds"][0])
    assert len(onco) == 1  # one comparison (a vs b)
    o = onco[0]
    assert o["comparison"] == "a_vs_b"
    m = o["matrix"]
    # columns: the 8 samples of a then the 8 of b (order within a block is
    # not part of the contract — assert membership, not position)
    assert set(m.columns[:8]) == {f"A{i}" for i in range(1, 9)}
    assert set(m.columns[8:]) == {f"B{i}" for i in range(1, 9)}
    assert len(o["blocks"]) == 2
    assert o["blocks"][0]["cohort"] == "a"
    assert len(o["blocks"][0]["columns"]) == 8
    assert o["blocks"][1]["cohort"] == "b"
    assert len(o["blocks"][1]["columns"]) == 8
    # values are 0/1/2 or NaN only
    vals = set(np.unique(m.to_numpy(dtype=float)[~np.isnan(m.to_numpy(dtype=float))]))
    assert vals <= {0.0, 1.0, 2.0}
    # meta rows align with the matrix rows and carry per-cohort stats
    assert (o["meta"].index == m.index).all()
    assert {"ptype", "cytoband", "chrom", "start", "end",
            "q_a", "q_b", "freq_a", "freq_b",
            "present_a", "present_b"} <= set(o["meta"].columns)
    # Amp rows sort above Del rows
    ptypes = list(o["meta"]["ptype"])
    assert ptypes == sorted(ptypes, reverse=True)
    # every peak with members in exactly one cohort has NaN on the other
    solo = o["meta"][~(o["meta"]["present_a"] & o["meta"]["present_b"])]
    for k in solo.index:
        side = "a" if solo.loc[k, "present_a"] else "b"
        other = "b" if side == "a" else "a"
        row = m.loc[k, [c for c in m.columns
                        if c in (o["blocks"][0 if other == "a" else 1]["columns"])]]
        assert row.isna().all()
    # cohort-side cells with a member peak are never NaN
    pres = o["meta"][o["meta"]["present_a"]]
    cols_a = o["blocks"][0]["columns"]
    assert not m.loc[pres.index, cols_a].isna().all().all()


def test_build_onco_data_empty_q(synth_root, pairwise_config):
    from gistic_compare import io_gistic
    from gistic_compare.config import assign_comparison_kinds, load_config

    cfg = load_config(pairwise_config)
    cfg["analysis"]["q_thresholds"] = [1e-12, 1e-12]
    data = io_gistic.load_all(cfg)
    assign_comparison_kinds(cfg, data)
    onco = build_onco_data(cfg, data, 1e-12)
    assert onco == []


def test_oncoplot_loci_merge_and_ranking(pairwise_config, tmp_path):
    """OncoPlot rows: overlapping amp + del chains merge into one locus;
    rows rank by combined frequency; per-cohort stacked bars; figure."""
    from gistic_compare import io_gistic
    from gistic_compare.config import assign_comparison_kinds, load_config
    from gistic_compare.pairwise_peak import build_onco_data
    from gistic_compare.plotting.oncoplot import (
        _locus_table,
        _merge_loci,
        plot_oncoplot,
    )

    cfg = load_config(pairwise_config)
    data = io_gistic.load_all(cfg)
    assign_comparison_kinds(cfg, data)
    o = build_onco_data(cfg, data, cfg["analysis"]["q_thresholds"][0])[0]

    groups = _merge_loci(o["meta"])
    # union chains are type-consistent: a multi-key group is a merged
    # amp + del locus, and every single peak belongs to exactly one group
    assert sum(len(g) for g in groups) == len(o["meta"])
    blocks = o["blocks"]
    for b in blocks:
        b["_n"] = len(b["columns"])
    frame = _locus_table(o["meta"], o["matrix"], blocks)
    assert len(frame) == len(groups)
    assert frame["mixed"].sum() == sum(1 for g in groups if len(g) > 1)
    assert frame["freq_union"].is_monotonic_decreasing
    for col in ("amp_0", "del_0", "amp_1", "del_1", "n_0", "n_1", "cell"):
        assert col in frame.columns
    # per-cohort frequencies are counts / cohort size
    agg = frame[[c for c in frame.columns if c.startswith(("amp_", "del_"))]]
    assert agg.to_numpy().min() >= 0 and agg.to_numpy().max() <= 1
    assert (frame["amp_0"] - frame["n_amp_0"] / frame["n_0"]).abs().max() \
        < 1e-12
    assert (frame["del_1"] - frame["n_del_1"] / frame["n_1"]).abs().max() \
        < 1e-12

    plot_oncoplot(o, cfg)
    name = "P1_oncoplot_a_vs_b"
    for ext in ("png", "pdf"):
        assert (cfg["output"]["_figures"] / f"{name}.{ext}").exists()
