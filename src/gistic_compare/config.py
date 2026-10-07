"""Load and validate the YAML configuration.

All paths in the config are resolved relative to the config file's directory,
so a run directory is self-contained and relocatable (no machine-absolute
paths required).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

from .genomes import (  # noqa: F401  (re-exported for convenience)
    CHROM_RANK,
    CHROMS,
    HG19_CHROM_SIZES,
    genome_build,
)
from .presets import _resolve, load_loci_set, load_paper_panel

# tab10 colors for dynamic cohort palettes (fallback for unknown cohort ids)
TAB10 = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
]

# Okabe-Ito colour-blind-safe qualitative palette (publication default)
OKABE_ITO = {
    "black": "#000000",
    "orange": "#E69F00",
    "sky_blue": "#56B4E9",
    "bluish_green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "reddish_purple": "#CC79A7",
}

# fixed Amp/Del semantics used by every figure
PTYPE_COLORS = {"Amp": OKABE_ITO["vermillion"], "Del": OKABE_ITO["blue"]}


def default_config_path() -> Path | None:
    """Search for a config in the current working directory."""
    for name in ("config.yaml", "config.yml", "config.gistic.yaml"):
        p = Path.cwd() / name
        if p.exists():
            return p
    return None


def _hex_rgb(h: str) -> tuple[float, float, float]:
    h = h.lstrip("#")
    return tuple(int(h[i:i + 2], 16) / 255.0 for i in (0, 2, 4))  # type: ignore


def _colour_dist(a: str, b: str) -> float:
    ra, ga, ba = _hex_rgb(a)
    rb, gb, bb = _hex_rgb(b)
    return ((ra - rb) ** 2 + (ga - gb) ** 2 + (ba - bb) ** 2) ** 0.5


def _pick_fallback_colour(taken: list[str]) -> str:
    """tab10 entry farthest from every colour already in use.

    Taking tab10 in config order can hand a cohort a colour nearly identical
    to an earlier one (tab10 orange #ff7f0e vs Okabe orange #E69F00), which
    makes two cohorts indistinguishable in every figure.
    """
    if not taken:
        return TAB10[0]
    return max(TAB10, key=lambda c: min(_colour_dist(c, t) for t in taken))


def build_comparisons(cfg: dict[str, Any]) -> list[dict[str, Any]]:
    """Create cfg['comparisons'] from comparison_mode / control.

    experiment_pairs: all unordered pairs of listed cohorts (no control);
                      kind starts as 'unknown' and is filled in after load
                      from sample-ID overlap (assign_comparison_kinds).
    control_vs_all:   star topology, a = control (fixed samples) for every row;
                      all kind=unpaired (control shares no samples with experiments).
    all_pairs:        every unordered pair; kind from cohort paired_with.
    explicit:         use cfg['comparisons'] as written.
    """
    mode = cfg.get("comparison_mode")
    cids = list(cfg["cohorts"])
    raw = cfg.get("comparisons")
    workflow = cfg.get("workflow", "benchmark")

    if mode is None:
        if workflow == "pairwise":
            mode = "experiment_pairs"
        else:
            mode = "explicit" if raw else "control_vs_all"
    cfg["comparison_mode"] = mode

    if mode == "explicit":
        if not raw:
            raise ValueError(
                "comparison_mode=explicit requires a non-empty 'comparisons' list")
        return list(raw)

    if mode == "all_pairs":
        out = []
        for i, a in enumerate(cids):
            for b in cids[i + 1:]:
                kind = "unpaired"
                if cfg["cohorts"][a].get("paired_with") == b:
                    kind = "paired"
                out.append({"a": a, "b": b, "kind": kind})
        return out

    if mode == "experiment_pairs":
        if len(cids) < 2:
            raise ValueError(
                "comparison_mode=experiment_pairs needs at least 2 cohorts")
        if cfg.get("control"):
            raise ValueError(
                "comparison_mode=experiment_pairs must not set 'control' "
                "(pairwise workflow has no control cohort)")
        return [
            {"a": a, "b": b, "kind": "unknown"}
            for i, a in enumerate(cids)
            for b in cids[i + 1:]
        ]

    if mode == "control_vs_all":
        control = cfg.get("control")
        if not control:
            raise ValueError("comparison_mode=control_vs_all requires 'control'")
        if control not in cfg["cohorts"]:
            raise ValueError(
                f"control '{control}' not in cohorts: {cids}")
        experiments = [c for c in cids if c != control]
        if not experiments:
            raise ValueError("control_vs_all needs at least one non-control cohort")
        # a = control always; control has fixed samples disjoint from experiments
        return [
            {"a": control, "b": exp, "kind": "unpaired"}
            for exp in experiments
        ]

    raise ValueError(
        f"unknown comparison_mode: {mode!r} "
        "(expected experiment_pairs | control_vs_all | all_pairs | explicit)")


def _setup_output(cfg: dict[str, Any], cfg_dir: Path,
                  outdir_override: str | Path | None) -> None:
    """Resolve where tables/figures go.

    Priority:  --outdir  >  output.results  >  ``results/`` next to the config.
    Inside the results root each run gets ``<run_id>/tables|figures`` unless
    ``output.tables`` / ``output.figures`` are set explicitly (relative to the
    config file or absolute).
    """
    out = cfg.setdefault("output", {})
    run_id = str(cfg.get("run_id", "default"))
    cfg["_run_id"] = run_id

    if outdir_override is not None:
        root = Path(outdir_override).expanduser().resolve()
    else:
        root = _resolve(out.get("results", "results"), cfg_dir)
    root.mkdir(parents=True, exist_ok=True)

    tables = out.get("tables")
    figures = out.get("figures")
    tdir = (Path(tables).expanduser() if tables
            else root / run_id / "tables")
    fdir = (Path(figures).expanduser() if figures
            else root / run_id / "figures")
    if not tdir.is_absolute():
        tdir = (cfg_dir / tdir)
    if not fdir.is_absolute():
        fdir = (cfg_dir / fdir)
    tdir.mkdir(parents=True, exist_ok=True)
    fdir.mkdir(parents=True, exist_ok=True)

    out["_tables"] = tdir.resolve()
    out["_figures"] = fdir.resolve()
    out["_run_id"] = run_id


def _setup_presets(cfg: dict[str, Any], cfg_dir: Path) -> None:
    """Attach preset dicts (loci set / paper panel); all optional."""
    regional = ((cfg.get("analysis") or {}).get("regional") or {})
    loci_set = regional.get("loci_set")
    cfg["_preset_loci"] = (load_loci_set(loci_set, cfg_dir)
                           if loci_set
                           else dict(id="", loci=[], bands=[], path=None))

    pb = cfg.get("paper_benchmark") or {}
    panel_path = pb.get("panel") if pb.get("enabled") else None
    cfg["_paper_panel"] = (load_paper_panel(panel_path, cfg_dir)
                           if panel_path else None)


def load_config(path: str | Path | None = None,
                workflow_override: str | None = None,
                outdir: str | Path | None = None) -> dict[str, Any]:
    cfg_path = (Path(path) if path else default_config_path())
    if cfg_path is None:
        raise FileNotFoundError(
            "no config given and no config.yaml in the current directory — "
            "pass --config path/to/config.yaml")
    cfg_path = cfg_path.expanduser().resolve()
    if not cfg_path.exists():
        raise FileNotFoundError(f"config not found: {cfg_path}")
    cfg_dir = cfg_path.parent
    with open(cfg_path) as fh:
        cfg = yaml.safe_load(fh)

    # genome build (hg19-only for this release, via genomes.genome_build)
    cfg["genome"] = genome_build(cfg.get("genome"))

    workflow = workflow_override or cfg.get("workflow") or "benchmark"
    if workflow not in ("benchmark", "pairwise"):
        raise ValueError(
            f"unknown workflow: {workflow!r} (expected benchmark | pairwise)")
    cfg["workflow"] = workflow

    # workflow ↔ control consistency
    if workflow == "benchmark" and not cfg.get("control"):
        raise ValueError("workflow=benchmark requires 'control'")
    if workflow == "pairwise" and cfg.get("control"):
        raise ValueError(
            "workflow=pairwise must not set 'control'")

    root = Path(cfg.get("data_root", "."))
    if not root.is_absolute():
        root = (cfg_dir / root).resolve()
    cfg["_config_path"] = str(cfg_path)
    cfg["_root"] = root

    # validate cohort paths
    missing = []
    for cid, c in cfg["cohorts"].items():
        p = root / c["path"]
        if not p.is_dir():
            missing.append(f"{cid}: {p}")
        c["_abs"] = p
    if missing:
        raise FileNotFoundError("cohort path(s) missing:\n  " + "\n  ".join(missing))

    # optional sample metadata (P3); missing path is allowed — module skips
    meta = cfg.get("metadata")
    if meta:
        if not meta.get("path"):
            raise ValueError("metadata.path is required when 'metadata' is set")
        meta_path = Path(meta["path"])
        if not meta_path.is_absolute():
            meta_path = (root / meta_path)
        meta["_abs"] = meta_path
        meta.setdefault("sample_id", "sample")
        meta.setdefault("covariates", [])

    # comparison topology
    cfg["comparisons"] = build_comparisons(cfg)

    # presets (study-specific loci / paper panel); attach before palette so
    # preset files can be kept next to the config
    _setup_presets(cfg, cfg_dir)

    # dynamic palette: cohort id -> colour
    #   1. explicit output.figure.cohort_colors in config
    #   2. farthest-from-used tab10 fallback (keeps cohorts distinguishable)
    cohort_colors = (((cfg.get("output") or {}).get("figure") or {})
                     .get("cohort_colors") or {})
    palette: dict[str, str] = {}
    for cid in cfg["cohorts"]:
        col = cohort_colors.get(cid) or _pick_fallback_colour(
            list(palette.values()))
        palette[cid] = col
    cfg["_palette"] = palette

    # results/<run_id>/{tables,figures} — never inside the installed package
    _setup_output(cfg, cfg_dir, outdir)

    return cfg


def assign_comparison_kinds(cfg: dict[str, Any],
                            data: dict[str, dict],
                            paired_frac: float = 0.8) -> None:
    """Fill kind='unknown' comparisons from sample-ID overlap.

    paired when |shared samples| >= paired_frac * min(n_a, n_b).
    Requires data[cid]['calls'] (loaded by io_gistic).
    """
    for cmp_ in cfg["comparisons"]:
        if cmp_.get("kind") != "unknown":
            continue
        ca, cb = cmp_["a"], cmp_["b"]
        sa = set(data[ca]["calls"]["sample"].unique())
        sb = set(data[cb]["calls"]["sample"].unique())
        if not sa or not sb:
            cmp_["kind"] = "unpaired"
            continue
        frac = len(sa & sb) / min(len(sa), len(sb))
        cmp_["kind"] = "paired" if frac >= paired_frac else "unpaired"
        cmp_["_sample_overlap_frac"] = frac


def cohort_ids(cfg: dict[str, Any]) -> list[str]:
    return list(cfg["cohorts"])


def out_table(cfg: dict[str, Any], name: str) -> Path:
    return cfg["output"]["_tables"] / name


def out_fig(cfg: dict[str, Any], name: str) -> Path:
    return cfg["output"]["_figures"] / name
