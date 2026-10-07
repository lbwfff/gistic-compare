"""Run manifest: a machine-readable record of each analysis run.

Every invocation writes ``results/<run_id>/run_manifest.json`` containing the
resolved configuration snapshot, package versions, and per-module status —
so any table in the run can be traced back to exactly what produced it.
"""

from __future__ import annotations

import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from . import __version__

_DEPENDENCIES = ("pandas", "numpy", "scipy", "matplotlib", "seaborn",
                 "statsmodels", "h5py", "pyyaml")


def _package_version(name: str) -> str | None:
    try:
        from importlib.metadata import version
        return version(name)
    except Exception:
        return None


def git_commit() -> str | None:
    """Short commit hash when the package is run from its git checkout."""
    try:
        out = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, timeout=5, check=True)
        return out.stdout.strip() or None
    except Exception:
        return None


def _sanitize(obj: Any) -> Any:
    """Make a config-like structure JSON-serialisable, preserving content."""
    if obj is None or isinstance(obj, (bool, int, float, str)):
        return obj
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, (list, tuple)):
        return [_sanitize(x) for x in obj]
    if isinstance(obj, set):
        return sorted(str(x) for x in obj)
    if isinstance(obj, dict):
        return {str(k): _sanitize(v) for k, v in obj.items()}
    import numpy as np
    if isinstance(obj, np.generic):
        return obj.item()
    return str(obj)


def build_manifest(cfg: dict[str, Any],
                   modules: list[dict[str, Any]] | None = None,
                   total_seconds: float | None = None,
                   data_summary: dict[str, Any] | None = None) -> dict[str, Any]:
    manifest: dict[str, Any] = dict(
        schema_version=1,
        tool="gistic-compare",
        tool_version=__version__,
        git_commit=git_commit(),
        timestamp=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        python_version=sys.version.split()[0],
        platform=platform.platform(),
        dependencies={d: _package_version(d) for d in _DEPENDENCIES},
        workflow=cfg.get("workflow"),
        run_id=cfg.get("_run_id"),
        genome=cfg.get("genome"),
        config_path=cfg.get("_config_path"),
        config=_sanitize({k: v for k, v in cfg.items()
                          if k not in ("output",)}),
        output=dict(
            tables=str(cfg["output"]["_tables"]),
            figures=str(cfg["output"]["_figures"]),
        ),
        modules=modules or [],
        total_seconds=(round(total_seconds, 1)
                       if total_seconds is not None else None),
    )
    if data_summary:
        manifest["cohorts"] = data_summary
    return manifest


def write_manifest(cfg: dict[str, Any],
                   modules: list[dict[str, Any]] | None = None,
                   total_seconds: float | None = None,
                   data_summary: dict[str, Any] | None = None) -> Path:
    """Write ``run_manifest.json`` into the run directory; return its path."""
    out_dir = Path(cfg["output"]["_tables"]).parent
    path = out_dir / "run_manifest.json"
    manifest = build_manifest(cfg, modules, total_seconds, data_summary)
    path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    return path
