"""Load preset (paper- / study-specific) configuration files.

Presets keep study-specific constants out of the source code:

* **loci sets**  — named genomic loci plus a list of "bookmark" cytobands that
  figures annotate (originally ``PAPER_LOCI`` / ``PAPER_TABLE2``).
* **reference panels** — numeric/qualitative reference values a cohort set is
  benchmarked against (originally the hard-coded ``PAPER`` dict of module G).

Both file formats are plain YAML; paths are resolved relative to the config
file, consistent with ``data_root``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


def _resolve(path: str | Path, base: Path) -> Path:
    p = Path(path)
    return p if p.is_absolute() else (base / p)


def load_loci_set(path: str | Path, base: Path) -> dict[str, Any]:
    """Load a loci-set preset.

    Expected schema::

        id: we_ea
        loci:                # named loci to test / annotate
          - {locus: 16p13.13, chrom: chr16, start: 10500000, end: 12600000,
             type: Amp, expected_cohort: EA}
        bands:               # cytobands bookmarked on figures
          - {band: 1p36.23, star: false}
    """
    p = _resolve(path, base)
    if not p.exists():
        raise FileNotFoundError(f"loci-set preset not found: {p}")
    with open(p) as fh:
        raw = yaml.safe_load(fh) or {}
    loci = []
    for r in raw.get("loci") or []:
        try:
            loci.append((str(r["locus"]), str(r["chrom"]),
                         int(r["start"]), int(r["end"]),
                         str(r["type"]), str(r.get("expected_cohort", ""))))
        except (KeyError, TypeError, ValueError) as e:
            raise ValueError(f"{p}: bad loci entry {r!r}: {e}") from e
    bands = []
    for r in raw.get("bands") or []:
        try:
            bands.append((str(r["band"]), bool(r.get("star", False))))
        except (TypeError, ValueError) as e:
            raise ValueError(f"{p}: bad bands entry {r!r}: {e}") from e
    return dict(id=str(raw.get("id", p.stem)), loci=loci, bands=bands,
                path=str(p))


def load_paper_panel(path: str | Path, base: Path) -> dict[str, Any]:
    """Load a reference panel for module G (paper bench).

    Expected schema::

        id: we_lung
        genome_pct: {amp: 15.5, del: 25.3}
        top_amp:     [{band: 5p15, min_freq_pct: 60.0}]
        top_del:     [{band: 3p14, min_freq_pct: 50.0}]
        amp_cytobands: ["7p11", "8q24"]
        del_cytobands: ["3p14", "8p21"]
        broad_amp_arms: ["1q", "5p"]
        focal_amp_arms: ["11q", "14q"]
        large_del_arms: ["3p", "6q"]
        windows:                # hg windows used when a cytoband prefix misses
          5p15: {chrom: chr5, start: 1500000, end: 10000000}
        drivers: ["EGFR", "MYC"]
        driver_groups: {E2F: ["E2F1", "E2F2", "E2F3"]}
    """
    p = _resolve(path, base)
    if not p.exists():
        raise FileNotFoundError(f"reference-panel preset not found: {p}")
    with open(p) as fh:
        raw = yaml.safe_load(fh) or {}

    gp = raw.get("genome_pct") or {}
    for key in ("amp", "del"):
        if key not in gp:
            raise ValueError(
                f"{p}: 'genome_pct.{key}' is required (%% of genome altered)")
    windows = {}
    for band, w in (raw.get("windows") or {}).items():
        windows[str(band)] = (str(w["chrom"]), int(w["start"]), int(w["end"]))

    panel = dict(
        id=str(raw.get("id", p.stem)),
        amp_genome_pct=float(gp["amp"]),
        del_genome_pct=float(gp["del"]),
        top_amp=[(str(r["band"]), float(r["min_freq_pct"]), str(r["band"]))
                 for r in (raw.get("top_amp") or [])],
        top_del=[(str(r["band"]), float(r["min_freq_pct"]), str(r["band"]))
                 for r in (raw.get("top_del") or [])],
        amp_cytobands=[str(x) for x in (raw.get("amp_cytobands") or [])],
        del_cytobands=[str(x) for x in (raw.get("del_cytobands") or [])],
        broad_amp_arms=[str(x) for x in (raw.get("broad_amp_arms") or [])],
        focal_amp_arms=[str(x) for x in (raw.get("focal_amp_arms") or [])],
        large_del_arms=[str(x) for x in (raw.get("large_del_arms") or [])],
        windows=windows,
        drivers=[str(x) for x in (raw.get("drivers") or [])],
        driver_groups={str(k): [str(v_) for v_ in v]
                       for k, v in (raw.get("driver_groups") or {}).items()},
        path=str(p),
    )
    return panel


def band_names(cfg: dict[str, Any] | None) -> set[str]:
    """Names of the bookmarked bands from the configured loci preset."""
    return {n for n, _s in from_cfg(cfg)["bands"]}


def loci_names(cfg: dict[str, Any] | None) -> list[str]:
    """Names of the loci from the configured loci preset."""
    return [n for n, *_ in from_cfg(cfg)["loci"]]


def from_cfg(cfg: dict[str, Any] | None) -> dict[str, Any]:
    """Return the preset dict recorded on ``cfg`` (empty when none)."""
    if not cfg:
        return dict(id="", loci=[], bands=[], path=None)
    return dict(cfg.get("_preset_loci") or
                dict(id="", loci=[], bands=[], path=None))


def loci_flag(loci: list[tuple], chrom: str, start: int, end: int) -> str:
    """Name of the first preset locus overlapping ``[start, end)``; else ''."""
    for name, c, lo, hi, _dirn, _who in loci:
        if chrom == c and end > lo and start < hi:
            return name
    return ""
