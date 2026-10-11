"""Command line interface for gistic-compare.

Run with::

    gistic-compare --config config.yaml --all
    python -m gistic_compare --config config.yaml --all
"""

from __future__ import annotations

import argparse
import sys
import time
import traceback

from . import (  # noqa: E402  # noqa: E402  # noqa: E402  # noqa: E402
    __version__,
    confounding,
    gene_sets,
    gs_score,
    io_gistic,  # noqa: E402
    novelty,
    pairwise_absolute,
    pairwise_amplitude,
    pairwise_arm,
    pairwise_burden,
    pairwise_calls,
    pairwise_diag,
    pairwise_marker,
    pairwise_model,
    pairwise_multivariate,
    pairwise_peak,
    pairwise_permutation,
    pairwise_regional,
    paper_reproduce,
    peak_overlap,
    plotting,  # noqa: E402
    pval_qc,
    sample_calls,
    score_track,
    seg_qc,
)
from .config import assign_comparison_kinds, load_config, out_table
from .run_manifest import write_manifest

BENCHMARK_MODULES = {
    "A": ("peak overlap", peak_overlap.run_module_a, plotting.figures_a),
    "B": ("score tracks", score_track.run_module_b, plotting.figures_b),
    "C": ("sample calls", sample_calls.run_module_c, plotting.figures_c),
    "D": ("gene sets", gene_sets.run_module_d, plotting.figures_d),
    "E": ("segment QC", seg_qc.run_module_e, plotting.figures_e),
    "F": ("novelty vs control", novelty.run_module_f, plotting.figures_f),
    "G": ("benchmark vs reference panel", paper_reproduce.run_module_g,
          plotting.figures_g),
}

PAIRWISE_MODULES = {
    "P1": ("peak set differences (detection)", pairwise_peak.run_module_p1,
           plotting.figures_p1),
    "P2": ("sample call concordance (detection)", pairwise_calls.run_module_p2,
           plotting.figures_p2),
    "P3": ("per-peak altered ~ covariate", pairwise_model.run_module_p3,
           plotting.figures_p3),
    "P4": ("regional between-cohort difference", pairwise_regional.run_module_p4,
           plotting.figures_p4),
    "P5": ("regional diff on GISTIC marker matrix", pairwise_marker.run_module_p5,
           plotting.figures_p5),
    "P6": ("arm-level broad landscape", pairwise_arm.run_module_p6,
           plotting.figures_p6),
    "P7": ("amplitude among altered samples", pairwise_amplitude.run_module_p7,
           plotting.figures_p7),
    "P8": ("multivariate cohort separation",
           pairwise_multivariate.run_module_p8, plotting.figures_p8),
    "P9": ("three-level SCNA burden + wGII",
           pairwise_burden.run_module_p9, plotting.figures_p9),
    "P10": ("burden-controlled permutation on peaks",
            pairwise_permutation.run_module_p10, plotting.figures_p10),
    "PQ": ("p-value QC (Q-Q + inflation)", pval_qc.run_module_pq,
           plotting.figures_pq),
    "PD": ("baseline/batch diagnostic (lambda inflation)",
           pairwise_diag.run_module_pd, plotting.figures_pd),
    "ABS": ("absolute copy-number comparisons (ASCAT tables)",
            pairwise_absolute.run_module_abs, plotting.figures_abs),
    "CF": ("confounder sensitivity (independent)", confounding.run_module_confounding,
           plotting.figures_confounding),
    "GS": ("G-score differential (permutation)", gs_score.run_module_gs,
           plotting.figures_gs),
}

WORKFLOWS = {"benchmark": BENCHMARK_MODULES, "pairwise": PAIRWISE_MODULES}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        prog="gistic-compare",
        description="Compare GISTIC 2.0 results across cohorts")
    ap.add_argument("--version", action="version",
                    version=f"gistic-compare {__version__}")
    ap.add_argument("--config", default=None, help="path to config yaml "
                    "(default: ./config.yaml)")
    ap.add_argument("--workflow", choices=["benchmark", "pairwise"],
                    default=None,
                    help="override the workflow from the config "
                         "(benchmark: A-G vs control; pairwise: P1-P5)")
    ap.add_argument("--modules", default=None,
                    help="comma-separated module IDs (default: all)")
    ap.add_argument("--all", action="store_true", help="run all modules "
                    "(default when --modules is not given)")
    ap.add_argument("--outdir", default=None,
                    help="output directory override "
                         "(default: 'results/<run_id>' next to the config)")
    ap.add_argument("--skip-figures", action="store_true")
    ap.add_argument("--fail-fast", action="store_true",
                    help="abort on first module error (default: continue)")
    ap.add_argument("--list-modules", action="store_true",
                    help="list available modules per workflow and exit")
    return ap.parse_args(argv)


def _marker_files_ok(cfg: dict) -> bool:
    """P5/GS need each cohort's marker matrix (D.cap1.5.mat)."""
    rcfg = (cfg.get("analysis") or {}).get("regional") or {}
    cap = rcfg.get("marker_cap_file", "D.cap1.5.mat")
    for cid, c in cfg["cohorts"].items():
        if not (c["_abs"] / f"{c['prefix']}.{cap}").exists():
            return False
    return True


def _skip_reasons(cfg: dict) -> dict[str, str]:
    """Modules that cannot run with this config, with a short reason."""
    reasons: dict[str, str] = {}
    if not cfg.get("_paper_panel"):
        reasons["G"] = ("no reference panel configured "
                        "(paper_benchmark.enabled + panel)")
    if cfg["workflow"] == "pairwise" and not _marker_files_ok(cfg):
        reasons["P5"] = "marker matrix (D.cap1.5.mat) missing in some cohort"
        reasons["GS"] = "marker matrix (D.cap1.5.mat) missing in some cohort"
    return reasons


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    if args.list_modules:
        for wf, mods in WORKFLOWS.items():
            print(f"workflow {wf}:")
            for mid, (label, _runner, _figs) in mods.items():
                print(f"  {mid:<3} {label}")
        return 0

    t0 = time.time()
    cfg = load_config(args.config, workflow_override=args.workflow,
                      outdir=args.outdir)
    workflow = cfg["workflow"]
    MODULES = WORKFLOWS[workflow]

    if args.modules:
        mods = [m.strip().upper() for m in args.modules.split(",") if m.strip()]
    else:
        mods = list(MODULES)
    bad = [m for m in mods if m not in MODULES]
    if bad:
        print(f"unknown module(s) for workflow={workflow}: {bad}; "
              f"choose from {list(MODULES)}", file=sys.stderr)
        other = [k for k, v in WORKFLOWS.items() if k != workflow
                 and any(m in v for m in bad)]
        if other:
            print(f"  (hint: {bad} may belong to workflow={other[0]})",
                  file=sys.stderr)
        return 2

    reasons = _skip_reasons(cfg)
    skipped_pre = [m for m in mods if m in reasons]
    mods = [m for m in mods if m not in reasons]
    modules_status: list[dict] = [
        dict(id=m, label=MODULES[m][0], status="skipped", reason=reasons[m])
        for m in skipped_pre
    ]
    for mid in skipped_pre:
        print(f"[skip] {mid}: {reasons[mid]}")
    if not mods:
        print("nothing to run", file=sys.stderr)
        write_manifest(cfg, modules_status, total_seconds=time.time() - t0)
        return 1

    print(f"[config] {cfg['_config_path']}")
    print(f"[config] genome={cfg['genome']}  workflow={workflow}  "
          f"run_id={cfg.get('_run_id', 'default')}  "
          f"mode={cfg.get('comparison_mode', '?')}"
          + (f"  control={cfg.get('control')}" if cfg.get("control") else ""))
    print(f"[output] {cfg['output']['_tables']}")
    print(f"[load] cohorts: {', '.join(cfg['cohorts'])}")
    t_load = time.time()
    data = io_gistic.load_all(cfg)

    # auto paired/unpaired for unknown kinds (experiment_pairs)
    paired_frac = float(cfg["analysis"].get("paired_frac", 0.8))
    assign_comparison_kinds(cfg, data, paired_frac=paired_frac)

    for cid, d in data.items():
        p = d["peaks"]
        role = " [CONTROL]" if cid == cfg.get("control") else ""
        print(f"  {cid}{role}: peaks={len(p)} amp={(p.type == 'Amp').sum()} "
              f"del={(p.type == 'Del').sum()} "
              f"samples={d['seg']['sample'].nunique()} "
              f"genes_rows={len(d['genes'])} scores={len(d['scores'])}")
    cmp_desc = ", ".join(
        f"{c['a']}_vs_{c['b']}({c.get('kind', 'unpaired')}"
        + (f",{c['_sample_overlap_frac']:.2f}"
           if "_sample_overlap_frac" in c else "")
        + ")"
        for c in cfg["comparisons"])
    print(f"[compare] {len(cfg['comparisons'])} pair(s): {cmp_desc}")
    print(f"[load] done in {time.time() - t_load:.1f}s")

    results: dict[str, dict] = {}
    failed: list[str] = []
    for mid in mods:
        label, runner, figs = MODULES[mid]
        print(f"\n=== Module {mid}: {label} ===")
        t = time.time()
        try:
            res = runner(cfg, data)
            results[mid] = res
            if not args.skip_figures:
                try:
                    figs(res, cfg)
                except Exception:
                    print(f"  [warn] figures failed for {mid}:")
                    traceback.print_exc()
            print(f"  done in {time.time() - t:.1f}s -> "
                  f"{cfg['output']['_tables']}/{mid}_*")
            modules_status.append(dict(id=mid, label=label, status="ok",
                                       seconds=round(time.time() - t, 1)))
        except Exception as e:
            failed.append(mid)
            print(f"  [ERROR] module {mid} failed: {e}")
            traceback.print_exc()
            modules_status.append(dict(id=mid, label=label, status="failed",
                                       error=str(e),
                                       seconds=round(time.time() - t, 1)))
            if args.fail_fast:
                break

    # overview table (always)
    try:
        _write_overview(cfg, data, results)
        print(f"\n[overview] {cfg['output']['_tables']}/00_cohort_overview.tsv")
    except Exception:
        traceback.print_exc()

    # run manifest (always) — config snapshot, versions, module status
    try:
        data_summary = {
            cid: dict(label=d["label"],
                      n_samples=int(d["seg"]["sample"].nunique()),
                      n_peaks=int(len(d["peaks"])),
                      n_amp=int((d["peaks"]["type"] == "Amp").sum()),
                      n_del=int((d["peaks"]["type"] == "Del").sum()))
            for cid, d in data.items()
        }
        manifest_path = write_manifest(cfg, modules_status,
                                       total_seconds=time.time() - t0,
                                       data_summary=data_summary)
        print(f"[manifest] {manifest_path}")
    except Exception:
        traceback.print_exc()

    print(f"\nTotal: {time.time() - t0:.1f}s; "
          f"modules ok={len(results)} "
          f"skipped={skipped_pre or 'none'} failed={failed or 'none'}")
    return 1 if failed else 0


def _write_overview(cfg: dict, data: dict, results: dict) -> None:
    import pandas as pd
    rows = []
    q_ts = cfg["analysis"]["q_thresholds"]
    for cid, d in data.items():
        p = d["peaks"]
        sc = d["seg_counts"]
        row = dict(
            cohort=cid,
            label=d["label"],
            n_samples=d["seg"]["sample"].nunique(),
            n_included=int((sc["included"] == "yes").sum()),
            n_peaks=len(p),
            n_amp=int((p["type"] == "Amp").sum()),
            n_del=int((p["type"] == "Del").sum()),
            n_segments=len(d["seg"]),
            median_segments=float(sc["segment_count"].median()),
        )
        for q in q_ts:
            row[f"n_peaks_q{q}"] = int((p["q"] < q).sum())
        # attach FGA if module E ran
        if "E" in results and not results["E"]["qc"].empty:
            qc = results["E"]["qc"]
            sub = qc[qc["cohort"] == cid]
            if len(sub):
                row["median_fga"] = float(sub["fga"].median())
                row["mean_fga"] = float(sub["fga"].mean())
        rows.append(row)
    pd.DataFrame(rows).to_csv(out_table(cfg, "00_cohort_overview.tsv"),
                              sep="\t", index=False)


if __name__ == "__main__":
    raise SystemExit(main())
