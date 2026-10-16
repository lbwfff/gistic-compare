"""Module ABS (pairwise): absolute copy-number comparisons (ASCAT tables).

Phase 1+2 of ``docs/plan_absolute_cn.md``.  The R wrapper
``R/absolute_cn.R`` (Phase 2) runs ASCAT 3.x on each cohort's marker-
level logR/BAF (SNP6 Birdseed chain) and emits the harmonized tables

  absolute_cn_<cid>.tsv    sample, purity, ploidy_abs, wgd_status,
                           wgd_frac, loh_frac, cnloh_frac, ...,
                           flag_low_purity, wgd_flag_b, wgd_disagreement,
                           germline_mode, het_source
  absolute_seg_<cid>.tsv   sample, chrom, start, end, cn_total, cn_major,
                           cn_minor

The Python module reads them (by default next to the cohort's GISTIC
outputs; ``absolute_cn.tables_dir`` overrides) and runs the full
between-cohort suite:

  WGD              2x2 Fisher (WGD+ vs WGD-) with odds ratio + CI; the
                   mean status (0/1/2) is reported, not tested (a FACETS
                   snapshot cannot count *multiple* WGD rounds reliably)
  purity           Mann-Whitney + rank-biserial
  ploidy_abs       same
  LOH / CN-LOH     same
  CN-space burden  P9's three-level event counting transplanted into
                   copy-number space: focal = allelic segments beyond the
                   neutral band (+/- ``cn_focal_band`` copies of CN 2) not
                   covered by broad events; broad arm event = length-
                   weighted mean CN over > ``cn_min_frac`` of the arm
                   beyond ``cn_arm_gain`` / below ``cn_arm_del``; both
                   arms same sign and |delta CN| < ``cn_match_diff`` merge
                   to chromosome level.  Amplitude weight |C|: 2 copies
                   beyond the high threshold (+/-``cn_high``), else 1.
  purity-adjusted peaks  per-peak ``altered ~ cohort`` vs ``altered ~
                   cohort + purity`` logistic (P3 machinery; purity is the
                   covariate only absolute CN can provide)

Low-purity QC: samples flagged ``flag_low_purity`` (purity <
``min_purity``, default 0.30 — R05 §4.5) are excluded from a sensitivity
rerun of the headline WGD test; the delta is recorded.

ASCAT itself is **not** part of this package's downstream analysis: it is
run on a server from this repository's ``R/absolute_cn.R`` (input prep:
``scripts/prepare_abs_inputs.py``), and only its two harmonized tables
cross the machine boundary.  The run book and the frozen table schema are
in ``docs/ascat_handoff.md``; the tables are validated against that schema
on read (a mismatch skips ABS with a pointer to the document instead of
silently computing on mislabelled columns).

Pipeline mode (opt-in, server-side convenience): with
``analysis.absolute_cn.auto_run: true`` and the tables missing, the R
wrapper is invoked here (subprocess, logged to ``ABS_r_pipeline.log``)
when ``raw_dir/manifest`` exists and Rscript + the ASCAT package are
available.  **Default false** — a missing table otherwise produces
``ABS_skipped.txt`` with instructions.

Germline modes (recorded in the output columns
``germline_mode``/``het_source``): matched (measured germline BAF), or
tumour-only via the samples' own Birdseed genotype calls (``het_calls``
manifest column) or ASCAT's built-in inference — tumour-only makes
``loh_frac``/``cnloh_frac`` lower-bound estimates (somatic-LOH loci drop
out of the het set).

Outputs: ABS_wgd.tsv, ABS_purity_ploidy.tsv, ABS_cn_scores.tsv,
         ABS_cn_tests.tsv, ABS_purity_adjusted_peaks.tsv,
         ABS_loh_del_cross.tsv, ABS_summary.tsv

Two LOH-specific analyses beyond the per-sample fractions, both of which
need the allele-specific segments (not just the per-sample metrics):

  WGD-stratified LOH  the pooled LOH test mixes a WGD+ stratum (minor
                      CN >= 1 is common, so fewer LOH segments survive)
                      with a WGD- stratum; its sign follows the cohorts'
                      WGD composition, so LOH is also compared *within*
                      each WGD stratum (`<comp> [WGD+ only]` /
                      `[WGD- only]` rows in the metrics table).
  LOH x Del cross     every (sample x Del peak) pair is classified
                      del+LOH / del-only / LOH-only / neither, splitting
                      LOH-only into copy-neutral (no dosage change —
                      structurally invisible to GISTIC's del track) and
                      CN-loss (below GISTIC's amplitude call); see
                      ``_loh_del_cross``.
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from .config import CHROMS, out_table
from .pairwise_burden import _arm_intervals, _load_bands
from .pairwise_model import build_union_keys
from .pairwise_regional import _regional_cfg

_DEFAULT_R_SCRIPT = (Path(__file__).resolve().parents[2] / "R"
                     / "absolute_cn.R")

# the frozen schema of the two harmonized tables (docs/ascat_handoff.md);
# required = what the statistics below actually read, optional = metadata
# the runner also emits (recorded in ABS_summary when present)
_CN_REQUIRED = ("sample", "purity", "ploidy_abs", "wgd_status", "wgd_frac",
                "loh_frac", "cnloh_frac", "flag_low_purity")
_SEG_REQUIRED = ("sample", "chrom", "start", "end", "cn_total", "cn_major",
                 "cn_minor")
_CN_OPTIONAL = ("fit_ok", "note", "germline_mode", "het_source",
                "n_segments_allelic", "wgd_flag_b", "wgd_disagreement",
                "non_aberrant")
_SCHEMA_DOC = "docs/ascat_handoff.md"


class TableSchemaError(RuntimeError):
    """Harmonized absolute tables do not match the frozen schema."""


def _abs_cfg(cfg: dict[str, Any], rcfg: dict) -> dict:
    a = (cfg.get("analysis") or {}).get("absolute_cn") or {}
    return dict(
        tool=str(a.get("tool", "ascat")),
        auto_run=bool(a.get("auto_run", False)),
        loh_cross_core=str(a.get("loh_cross_core", "narrow")),
        loh_cross_cover=float(a.get("loh_cross_cover", 0.5)),
        tables_dir=(a.get("tables_dir") or None),
        raw_dir=str(a.get("raw_dir", "raw_snp6")),
        manifest=str(a.get("manifest", "absolute_manifest.tsv")),
        r_script=str(a.get("r_script", _DEFAULT_R_SCRIPT)),
        penalty=float(a.get("penalty", 70)),
        gamma=float(a.get("gamma", 0.55)),
        max_ploidy=float(a.get("max_ploidy", 5.5)),
        gc_file=str(a.get("gc_file") or ""),
        cache_dir=str(a.get("cache_dir") or ""),
        jobs=max(1, int(a.get("jobs", 4))),
        r_timeout_min_per_sample=float(a.get("r_timeout_min_per_sample",
                                             90.0)),
        cn_thr_wgd=float(a.get("cn_thr_wgd", 3.0)),
        loh_frac_wgd=float(a.get("loh_frac_wgd", 0.30)),
        wgd_class_frac=float(a.get("wgd_class_frac", 0.30)),
        wgd_multi_frac=float(a.get("wgd_multi_frac", 0.10)),
        min_purity=float(a.get("min_purity", 0.30)),
        ploidy_wgd=float(a.get("ploidy_wgd", 3.7)),
        cn_arm_gain=float(a.get("cn_arm_gain", 3.5)),
        cn_arm_del=float(a.get("cn_arm_del", 1.5)),
        cn_focal_band=float(a.get("cn_focal_band", 1.0)),
        cn_high=float(a.get("cn_high", 2.0)),
        cn_min_frac=float(a.get("cn_min_frac", 0.70)),
        cn_match_diff=float(a.get("cn_match_diff", 0.5)),
        fdr=float(a.get("fdr", rcfg["fdr"])),
        seed=int(a.get("seed", rcfg.get("seed", 42))),
    )


# ---------------------------------------------------------------------------
# table-mode input
# ---------------------------------------------------------------------------

def _norm_chrom(series: pd.Series) -> pd.Series:
    """'1' / '1.0' / 'chr1' -> 'chr1' (the package's chromosome naming)."""
    s = series.astype(str).str.strip()
    plain = s.str.fullmatch(r"\d+(\.0+)?")
    return s.mask(plain, "chr" + s.str.split(".").str[0])


def _schema_check(df: pd.DataFrame, required: tuple[str, ...],
                  path: Path, what: str) -> None:
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise TableSchemaError(
            f"{what} table {path.name} is missing required column(s) "
            f"{', '.join(missing)} — the schema is frozen; ask the ASCAT run "
            f"for {', '.join(required)} (see {_SCHEMA_DOC})")
    bad_chrom = sorted(set(df["chrom"].dropna().unique()) - set(CHROMS)) \
        if "chrom" in df.columns else []
    if bad_chrom:                       # '_norm_chrom' ran first, so this is
        raise TableSchemaError(        # a genome/naming mismatch, not 'chrN'
            f"{what} table {path.name} has chromosome label(s) "
            f"{', '.join(map(str, bad_chrom[:5]))} outside the configured "
            f"genome ({CHROMS[0]}..{CHROMS[-1]}) — check the build and the "
            f"naming convention in {_SCHEMA_DOC}")


def _read_abs_table(path: Path, required: tuple[str, ...],
                    what: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        raise TableSchemaError(f"{what} table {path} is empty")
    for col in ("sample", "chrom"):
        if col in df.columns:
            df[col] = df[col].astype(str).str.strip()
    if "chrom" in df.columns:
        df["chrom"] = _norm_chrom(df["chrom"])
    _schema_check(df, required, path, what)
    for col in required:
        if col in ("sample", "chrom", "note"):
            continue
        num = pd.to_numeric(df[col], errors="coerce")
        if num.isna().all():
            raise TableSchemaError(
                f"{what} table {path.name}: column '{col}' is not numeric "
                f"— see {_SCHEMA_DOC}")
        df[col] = num
    return df


def _load_abs_tables(cfg: dict[str, Any], data: dict[str, dict],
                     acfg: dict) -> tuple[dict, dict, list[str]]:
    """Per-cohort absolute tables.  Default location: the cohort directory
    (next to the GISTIC outputs); ``absolute_cn.tables_dir`` overrides.

    Both tables are validated against the frozen schema (see the module
    docstring); a missing table is reported in ``missing``, a malformed one
    raises :class:`TableSchemaError`."""
    tables_dir = (Path(acfg["tables_dir"]).expanduser()
                  if acfg["tables_dir"] else None)
    cn, seg, missing = {}, {}, []
    for cid, c in cfg["cohorts"].items():
        base = tables_dir if tables_dir is not None else c["_abs"]
        p_cn = base / f"absolute_cn_{cid}.tsv"
        p_seg = base / f"absolute_seg_{cid}.tsv"
        if p_cn.exists() and p_seg.exists():
            cn[cid] = _read_abs_table(p_cn, _CN_REQUIRED, "absolute CN")
            seg[cid] = _read_abs_table(p_seg, _SEG_REQUIRED, "allelic segment")
            absent = [c for c in _CN_OPTIONAL if c not in cn[cid].columns]
            if absent:
                print(f"  [ABS] note: {cid} absolute CN table lacks "
                      f"optional column(s) {', '.join(absent)}")
        else:
            missing.append(cid)
    return cn, seg, missing


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _fdr(p: np.ndarray) -> np.ndarray:
    out = np.full_like(p, np.nan, dtype=float)
    ok = np.isfinite(p)
    if ok.any():
        from statsmodels.stats.multitest import multipletests
        out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out


def _rank_biserial(xa: np.ndarray, xb: np.ndarray) -> float:
    try:
        u1 = stats.mannwhitneyu(xa, xb, alternative="two-sided")[0]
    except ValueError:
        return np.nan
    return 2.0 * u1 / (len(xa) * len(xb)) - 1.0


def _wgd_test(pool: pd.DataFrame, ca: str, cb: str, comp: str) -> dict | None:
    d = pool.dropna(subset=["wgd_status"])
    a, b = d[d["cohort"] == ca], d[d["cohort"] == cb]
    ka, kb = int((a["wgd_status"] > 0).sum()), int((b["wgd_status"] > 0).sum())
    na, nb = len(a), len(b)
    if na == 0 or nb == 0:
        return None
    table = [[ka, na - ka], [kb, nb - kb]]
    or_, p = stats.fisher_exact(table, alternative="two-sided")
    lo = hi = np.nan
    if np.isfinite(or_):
        try:
            from statsmodels.stats.contingency_tables import Table2x2
            t2 = Table2x2(np.asarray(table))
            lo, hi = t2.oddsratio_confint()
        except Exception:
            pass
    else:
        # degenerate 2x2 (a zero cell): Haldane-Anscombe 0.5 correction
        t2c = np.asarray(table, dtype=float) + 0.5
        or_ = (t2c[0, 0] * t2c[1, 1]) / (t2c[0, 1] * t2c[1, 0])
        try:
            from statsmodels.stats.contingency_tables import Table2x2
            t2 = Table2x2(t2c)
            lo, hi = t2.oddsratio_confint()
        except Exception:
            pass
    try:
        p_mw = float(stats.mannwhitneyu(a["wgd_status"].to_numpy(float),
                                        b["wgd_status"].to_numpy(float),
                                        alternative="two-sided")[1])
    except ValueError:
        p_mw = np.nan
    return dict(comparison=comp, n_a=na, n_b=nb, wgd_n_a=ka, wgd_n_b=kb,
                wgd_frac_a=round(ka / na, 4), wgd_frac_b=round(kb / nb, 4),
                odds_ratio=(round(float(or_), 4) if np.isfinite(or_)
                            else np.nan),
                or_lo=(round(float(lo), 4) if np.isfinite(lo) else np.nan),
                or_hi=(round(float(hi), 4) if np.isfinite(hi) else np.nan),
                p_fisher=float(p), p_mw_status=p_mw,
                mean_status_a=round(float(a["wgd_status"].mean()), 4),
                mean_status_b=round(float(b["wgd_status"].mean()), 4))


def _rscript() -> str | None:
    return shutil.which("Rscript")


def _ascat_available(rscript: str) -> bool:
    return _r_pkg_available(rscript, "ASCAT")


def _r_pkg_available(rscript: str, pkg: str) -> bool:
    try:
        r = subprocess.run(
            [rscript, "-e",
             f'cat(requireNamespace("{pkg}", quietly=TRUE))'],
            capture_output=True, text=True, timeout=120)
    except (OSError, subprocess.TimeoutExpired):
        return False
    return "TRUE" in r.stdout


def _run_r_pipeline(cfg: dict[str, Any], acfg: dict,
                    missing: list[str]) -> tuple[bool, str]:
    """Run the ASCAT wrapper (``R/absolute_cn.R``) for cohorts whose
    harmonized tables are missing — opt-in via ``auto_run`` (see
    ``docs/ascat_handoff.md``); in normal operation the tables come from
    the server-side run.  Returns (ok, message)."""
    rscript = _rscript()
    if rscript is None:
        return False, ("Rscript not on PATH — install R or place harmonized "
                       "absolute_cn_<cid>.tsv / absolute_seg_<cid>.tsv "
                       "tables per cohort")
    script = Path(acfg["r_script"]).expanduser()
    if not script.exists():
        return False, (f"R wrapper not found: {script} — set "
                       "analysis.absolute_cn.r_script or place harmonized "
                       "tables per cohort")
    if not _ascat_available(rscript):
        return False, ("R package 'ASCAT' not available — install via "
                       "remotes::install_github('VanLoo-lab/ascat') (the "
                       "package directory is ascat-master/ASCAT, and its "
                       "unused 'copynumber' Import may be removed) or place "
                       "harmonized tables per cohort")
    manifest = Path(cfg["_root"]) / acfg["raw_dir"] / acfg["manifest"]
    if not manifest.exists():
        return False, (f"manifest not found: {manifest} — pipeline mode "
                       "needs analysis.absolute_cn.raw_dir/manifest "
                       "(columns: sample, cohort, logr, baf [, "
                       "germline_logr, germline_baf | het_calls]) or "
                       "harmonized tables per cohort")
    tdir = Path(cfg["output"]["_tables"])
    jobs = max(1, int(acfg.get("jobs", 4)))
    timeout_min = float(acfg.get("r_timeout_min_per_sample", 90.0))
    log: list[str] = []
    for cid in missing:
        c = cfg["cohorts"][cid]
        out_dir = (Path(acfg["tables_dir"]).expanduser()
                   if acfg["tables_dir"] else c["_abs"])
        man_all = pd.read_csv(manifest, sep="\t")
        rows = man_all[man_all["cohort"].astype(str).str.lower()
                       == cid.lower()]
        if rows.empty:
            return False, (f"manifest has no rows for cohort {cid} "
                           f"({manifest})")
        # round-robin split into per-worker manifests, written next to the
        # source manifest so relative paths and the markers grid resolve
        n_chunks = min(jobs, len(rows))
        parts = [rows.iloc[i::n_chunks] for i in range(n_chunks)]
        chunk_dir = out_dir / "ASCAT_chunks"
        chunk_dir.mkdir(parents=True, exist_ok=True)
        log.append(f"## cohort {cid}: {len(rows)} sample(s) on "
                   f"{len(parts)} worker(s)\n")
        procs = []
        for i, part in enumerate(parts, 1):
            cm = manifest.parent / f"{manifest.stem}__{cid}_chunk{i}.tsv"
            part.to_csv(cm, sep="\t", index=False)
            chunk_out = chunk_dir / f"chunk{i}"
            cmd = [rscript, str(script), "--manifest", str(manifest.parent
                                                           / cm.name),
                   "--cohort", cid, "--out-dir", str(chunk_out),
                   "--penalty", str(acfg["penalty"]),
                   "--gamma", str(acfg["gamma"]),
                   "--seed", str(acfg["seed"]),
                   "--cn-thr-wgd", str(acfg["cn_thr_wgd"]),
                   "--loh-frac-wgd", str(acfg["loh_frac_wgd"]),
                   "--min-purity", str(acfg["min_purity"]),
                   "--ploidy-wgd", str(acfg["ploidy_wgd"]),
                   "--max-ploidy", str(acfg["max_ploidy"])]
            if acfg["gc_file"]:
                gc_p = Path(acfg["gc_file"])
                if not gc_p.is_absolute():
                    gc_p = (Path(cfg["_root"]) / gc_p).resolve()
                cmd += ["--gc-file", str(gc_p)]
            if acfg["cache_dir"]:
                cache_p = Path(acfg["cache_dir"])
                if not cache_p.is_absolute():
                    cache_p = (Path(cfg["_root"]) / cache_p).resolve()
                cache_p.mkdir(parents=True, exist_ok=True)
                cmd += ["--cache-rds", str(cache_p)]
            log.append(f"$ [{' '.join(cm.name.split('_chunk')[-1:])}] "
                       f"{len(part)} samples: {' '.join(cmd[1:])}\n")
            timeout = 600.0 + timeout_min * 60.0 * len(part)
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, text=True)
            procs.append((i, proc, timeout, len(part), chunk_out, cm))
        # ---- wait for all workers ----
        failed: list[str] = []
        cn_parts: list[pd.DataFrame] = []
        seg_parts: list[pd.DataFrame] = []
        for i, proc, timeout, n_s, chunk_out, cm in procs:
            tag = f"worker {i} ({n_s} samples)"
            try:
                out, _ = proc.communicate(timeout=timeout)
            except subprocess.TimeoutExpired:
                proc.kill()
                proc.communicate()
                out = ""
                failed.append(f"{tag}: TIMEOUT after {timeout:.0f}s")
                log.append(f"## {tag}: TIMEOUT\n")
                continue
            log.append(f"## {tag} (rc={proc.returncode})\n{out}\n")
            if proc.returncode != 0:
                failed.append(f"{tag}: rc={proc.returncode}")
                continue
            cn_f = chunk_out / f"absolute_cn_{cid}.tsv"
            seg_f = chunk_out / f"absolute_seg_{cid}.tsv"
            if cn_f.exists():
                cn_parts.append(pd.read_csv(cn_f, sep="\t"))
            if seg_f.exists():
                seg_parts.append(pd.read_csv(seg_f, sep="\t"))
        (tdir / "ABS_r_pipeline.log").write_text("\n".join(log),
                                                 encoding="utf-8")
        if not cn_parts:
            return False, (f"ASCAT pipeline failed for {cid} — see "
                           f"{tdir / 'ABS_r_pipeline.log'}")
        cn = pd.concat(cn_parts, ignore_index=True)
        cn.to_csv(out_dir / f"absolute_cn_{cid}.tsv", sep="\t", index=False)
        seg = (pd.concat(seg_parts, ignore_index=True) if seg_parts else
               pd.DataFrame(columns=["sample", "chrom", "start", "end",
                                     "cn_total", "cn_major", "cn_minor"]))
        seg.to_csv(out_dir / f"absolute_seg_{cid}.tsv", sep="\t",
                   index=False)
        # clean the per-worker work areas on success
        shutil.rmtree(chunk_dir, ignore_errors=True)
        for cm in procs and [p[5] for p in procs]:
            cm.unlink(missing_ok=True)
        if failed:
            log.append(f"## merged with failed workers: {failed}\n")
            (tdir / "ABS_r_pipeline.log").write_text("\n".join(log),
                                                     encoding="utf-8")
        n_ok = int((cn["fit_ok"] == 1).sum()) if "fit_ok" in cn else 0
        print(f"  [ABS] {cid}: {n_ok}/{len(cn)} samples fitted "
              f"({len(parts)} workers"
              + (f"; {len(failed)} worker failure(s)" if failed else "")
              + ")")
    (tdir / "ABS_r_pipeline.log").write_text("\n".join(log),
                                             encoding="utf-8")
    return True, (f"ran ASCAT on {len(missing)} cohort(s) "
                  f"({jobs} workers; manifest {manifest})")


# ---------------------------------------------------------------------------
# CN-space three-level burden (P9 transplanted into copy-number units)
# ---------------------------------------------------------------------------

def _arm_events_cn(seg_allelic: pd.DataFrame, arms: pd.DataFrame,
                   acfg: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Broad arm/chromosome events in copy-number space."""
    gain, dell = acfg["cn_arm_gain"], acfg["cn_arm_del"]
    min_frac = acfg["cn_min_frac"]
    samples = sorted(seg_allelic["sample"].unique())
    sidx = {s: i for i, s in enumerate(samples)}
    arm_chrom = arms["chrom"].to_numpy()
    arm_start = arms["start"].to_numpy()
    arm_end = arms["end"].to_numpy()
    arm_len = np.maximum(arm_end - arm_start, 1.0).astype(float)
    cov = np.zeros((len(arms), len(samples)))
    vsum = np.zeros((len(arms), len(samples)))
    for c in CHROMS:
        sel = np.flatnonzero(arm_chrom == c)
        if len(sel) == 0:
            continue
        sub = seg_allelic[seg_allelic["chrom"] == c]
        if sub.empty:
            continue
        a_s = arm_start[sel].astype(float)
        a_e = arm_end[sel].astype(float)
        for r in sub.itertuples():
            lo = np.maximum(a_s, r.start)
            hi = np.minimum(a_e, r.end)
            ov = np.clip(hi - lo, 0.0, None)
            if not ov.any():
                continue
            j = sidx[r.sample]
            cov[sel, j] += ov
            vsum[sel, j] += ov * r.cn_total
    with np.errstate(invalid="ignore", divide="ignore"):
        mean_cn = np.where(cov > 0, vsum / np.maximum(cov, 1.0), np.nan)
        frac = cov / arm_len[:, None]
    arm_rows = []
    for ai, arm_row in enumerate(arms.itertuples()):
        for si, s in enumerate(samples):
            m, f = mean_cn[ai, si], frac[ai, si]
            if not np.isfinite(m) or f <= min_frac:
                continue
            if m >= gain:
                direction = "Amp"
            elif m <= dell:
                direction = "Del"
            else:
                continue
            arm_rows.append(dict(
                sample=s, arm=str(arm_row.arm), chrom=arm_row.chrom,
                start=int(arm_row.start), end=int(arm_row.end),
                direction=direction, mean_cn=float(m), frac=float(f)))
    arm_ev = pd.DataFrame(arm_rows, columns=[
        "sample", "arm", "chrom", "start", "end", "direction", "mean_cn",
        "frac"])

    # chromosome-level merge (same sign, mean CN within match)
    chrom_rows = []
    arm_ev = arm_ev.reset_index(drop=True)
    arm_ev["_id"] = arm_ev.index
    merged_ids: list[int] = []
    for (sample, chrom, dirn), g in arm_ev.groupby(
            ["sample", "chrom", "direction"], sort=False):
        if len(g) < 2:
            continue
        p = g[g["arm"].str.endswith("p")]
        q = g[g["arm"].str.endswith("q")]
        if p.empty or q.empty:
            continue
        rp, rq = p.iloc[0], q.iloc[0]
        if abs(rp["mean_cn"] - rq["mean_cn"]) < acfg["cn_match_diff"]:
            w = rp["frac"] + rq["frac"]
            chrom_rows.append(dict(
                sample=sample, chrom=chrom,
                start=int(min(rp["start"], rq["start"])),
                end=int(max(rp["end"], rq["end"])),
                direction=dirn,
                mean_cn=float((rp["mean_cn"] * rp["frac"]
                               + rq["mean_cn"] * rq["frac"])
                              / max(w, 1e-12)),
                arms=f"{rp['arm']}+{rq['arm']}"))
            merged_ids += [int(rp["_id"]), int(rq["_id"])]
    arm_out = arm_ev[~arm_ev["_id"].isin(merged_ids)].drop(columns="_id")
    chrom_ev = pd.DataFrame(chrom_rows, columns=[
        "sample", "chrom", "start", "end", "direction", "mean_cn", "arms"])
    return arm_out, chrom_ev


def _focal_events_cn(seg_one: pd.DataFrame, broad_spans: list[dict],
                     acfg: dict) -> list[dict]:
    """Focal copy-number events: allelic segments beyond the neutral band
    (|cn_total - 2| > ``cn_focal_band``) not covered by broad events."""
    band, high = acfg["cn_focal_band"], acfg["cn_high"]
    events = []
    for r in seg_one.itertuples():
        dist = abs(r.cn_total - 2.0)
        if dist <= band:
            continue
        C = (2.0 if dist >= high else 1.0) * (1 if r.cn_total > 2 else -1)
        covered, length = 0.0, max(r.end - r.start, 1)
        ev_dir = "Amp" if r.cn_total > 2 else "Del"
        for sp in broad_spans:
            if sp["chrom"] != r.chrom or sp["direction"] != ev_dir:
                continue
            ov = min(r.end, sp["end"]) - max(r.start, sp["start"])
            if ov > 0:
                covered += ov
        if covered / length > 0.5:
            continue
        events.append(dict(chrom=r.chrom, start=int(r.start),
                           end=int(r.end), C=C))
    return events


def _cn_scores(cfg: dict[str, Any], seg: dict[str, pd.DataFrame],
               arms: pd.DataFrame, acfg: dict) -> pd.DataFrame:
    score_frames = []
    for cid, segc in seg.items():
        arm_ev, chrom_ev = _arm_events_cn(segc, arms, acfg)
        broad_spans = [dict(sample=r["sample"], chrom=r["chrom"],
                            direction=r["direction"], start=r["start"],
                            end=r["end"])
                       for r in arm_ev.to_dict("records")]
        broad_spans += [dict(sample=r["sample"], chrom=r["chrom"],
                             direction=r["direction"], start=r["start"],
                             end=r["end"])
                        for r in chrom_ev.to_dict("records")]
        dist_arm = (arm_ev["mean_cn"] - 2.0).abs()
        arm_ev["C"] = np.sign(arm_ev["mean_cn"] - 2.0) * np.where(
            dist_arm >= acfg["cn_high"], 2.0, 1.0)
        dist_ch = (chrom_ev["mean_cn"] - 2.0).abs()
        chrom_ev["C"] = np.sign(chrom_ev["mean_cn"] - 2.0) * np.where(
            dist_ch >= acfg["cn_high"], 2.0, 1.0)
        samples = sorted(segc["sample"].unique())
        focal_rows = []
        for s in samples:
            seg_one = segc[segc["sample"] == s]
            spans = [sp for sp in broad_spans if sp["sample"] == s]
            ev = _focal_events_cn(seg_one, spans, acfg)
            focal_rows.append(dict(
                sample=s, n_focal=len(ev),
                s_focal=float(np.abs([e["C"] for e in ev]).sum())
                if ev else 0.0))
        focal_map = {r["sample"]: r["s_focal"] for r in focal_rows}
        nfocal_map = {r["sample"]: r["n_focal"] for r in focal_rows}
        sc = pd.DataFrame({"sample": samples})
        sc["S_focal_cn"] = sc["sample"].map(focal_map).fillna(0.0)
        sc["n_focal_events"] = sc["sample"].map(nfocal_map).fillna(0).astype(int)
        sc["S_arm_cn"] = sc["sample"].map(
            arm_ev.groupby("sample")["C"].apply(
                lambda x: np.abs(x).sum())).fillna(0.0)
        sc["S_chrom_cn"] = sc["sample"].map(
            chrom_ev.groupby("sample")["C"].apply(
                lambda x: np.abs(x).sum())).fillna(0.0)
        sc.insert(0, "cohort", cid)
        score_frames.append(sc)
    scores = pd.concat(score_frames, ignore_index=True)
    # pooled rank normalization per comparison
    for cmp_ in cfg["comparisons"]:
        ca, cb = cmp_["a"], cmp_["b"]
        if cmp_.get("kind") == "paired":
            continue
        pool = scores["cohort"].isin([ca, cb])
        for lvl in ("S_focal_cn", "S_arm_cn", "S_chrom_cn"):
            scores.loc[pool, f"{lvl}_star"] = (
                scores.loc[pool, lvl].rank(method="average")
                / max(int(pool.sum()), 1))
    return scores


# ---------------------------------------------------------------------------
# purity-adjusted per-peak models
# ---------------------------------------------------------------------------

def _purity_adjusted_peaks(cfg: dict[str, Any], data: dict[str, dict],
                           cn_tables: dict) -> pd.DataFrame:
    """Per-peak altered ~ cohort vs altered ~ cohort + purity (P3
    machinery; purity is the covariate only absolute CN can provide)."""
    import statsmodels.api as sm

    cn_thr = float(((cfg.get("analysis") or {}).get("permtest") or {})
                   .get("cn_thr", 0.25))
    keys = build_union_keys(cfg, data, cfg["analysis"]["q_thresholds"][0])
    if keys.empty:
        return pd.DataFrame()

    purity: dict[str, pd.Series] = {}
    for cid, t in cn_tables.items():
        purity[cid] = t.set_index("sample")["purity"]

    rows = []
    for cmp_ in cfg["comparisons"]:
        ca, cb = cmp_["a"], cmp_["b"]
        if cmp_.get("kind") == "paired":
            continue
        comp = f"{ca}_vs_{cb}"
        built = _binary_matrix_cfg(cfg, data, ca, cb, keys, cn_thr)
        if built is None:
            continue
        X, samples, meta = built
        cohorts = [ca if s in set(data[ca]["peak_cn"]["sample"]) else cb
                   for s in samples]
        pur = np.array([
            (purity[ca].get(s, np.nan) if c == ca
             else purity[cb].get(s, np.nan))
            for s, c in zip(samples, cohorts)], dtype=float)
        keep = np.isfinite(pur)
        if keep.sum() < 10:
            continue
        Xo = X[keep]
        is_bo = (np.array(cohorts)[keep] == cb).astype(float)
        puro = pur[keep]

        def _fit(y: np.ndarray, x: np.ndarray) -> tuple[float, float]:
            if len(np.unique(y)) < 2 or len(np.unique(x)) < 2:
                return np.nan, np.nan
            try:
                fit = sm.GLM(y, np.column_stack([np.ones(len(y)), x]),
                             family=sm.families.Binomial()).fit(disp=0)
                return float(np.exp(fit.params[1])), float(fit.pvalues[1])
            except Exception:
                return np.nan, np.nan

        for j, r in meta.iterrows():
            y = Xo[:, j].astype(float)
            crude_or, crude_p = _fit(y, is_bo)
            if not np.isfinite(crude_or):
                continue
            adj_or, adj_p = _fit(y, np.column_stack([is_bo, puro]))
            rows.append(dict(
                comparison=comp, peak_key=r["peak_key"], ptype=r["ptype"],
                cytoband=r["cytoband"], n=int(keep.sum()),
                crude_or=round(crude_or, 4),
                crude_p=round(crude_p, 5),
                adj_or_purity=round(adj_or, 4),
                adj_p_purity=round(adj_p, 5),
            ))
    out = pd.DataFrame(rows)
    if not out.empty and "adj_p_purity" in out.columns:
        out["adj_q_purity"] = _fdr(out["adj_p_purity"].to_numpy(float))
    return out


def _binary_matrix_cfg(cfg: dict[str, Any], data: dict[str, dict],
                       ca: str, cb: str, keys: pd.DataFrame,
                       cn_thr: float):
    from .pairwise_permutation import _binary_matrix as _bm
    return _bm(cfg, data, ca, cb, keys, cn_thr)


# ---------------------------------------------------------------------------
# LOH x GISTIC-Del cross (mechanism, not frequency)
# ---------------------------------------------------------------------------

def _loh_del_cross(cfg: dict[str, Any], data: dict[str, dict],
                   seg_tables: dict[str, pd.DataFrame],
                   cn_tables: dict[str, pd.DataFrame],
                   acfg: dict) -> pd.DataFrame:
    """Classify every (sample x Del peak) pair as del+LOH / del-only /
    LOH-only / neither, per cohort.

    GISTIC calls deletion from the **total** copy number at a peak, so a
    locus can lose an allele without any loss of dosage (copy-neutral LOH,
    ``cn_total == 2 & cn_minor == 0``) and stay invisible to the del track.
    Intersecting the allele-specific segments with GISTIC's own per-sample
    del calls shows how much of the LOH at those loci the del frequency
    misses.  A pair counts as LOH when at least ``loh_cross_cover`` of the
    peak's core (``loh_cross_core``: narrow Peak Limits by default) is
    covered by ``cn_minor == 0`` segments; the LOH-only class is split into
    CN-neutral (no dosage change — structurally invisible to GISTIC) and
    CN-loss (dosage lost but below GISTIC's amplitude call) counts.

    The four classes partition the **covered** samples (``n_covered``;
    segments tile the genome, so this equals ``n_samples`` unless a cohort
    has a gap at that locus), and the fractions use ``n_samples`` — the
    cohort's samples with a successful fit, i.e. GISTIC's own denominator
    (``freq_altered``, reported as ``gistic_freq``) unless a fit failed.
    """
    q_main = float(((cfg.get("analysis") or {}).get("q_thresholds")
                    or [0.25])[0])
    core = "wide" if str(acfg.get("loh_cross_core", "narrow")) == "wide" \
        else "narrow"
    cover = float(acfg.get("loh_cross_cover", 0.5))
    fam = {c for cmp_ in cfg["comparisons"] for c in (cmp_["a"], cmp_["b"])}
    rows = []
    for cid, seg in (seg_tables or {}).items():
        if cid not in fam or seg is None or seg.empty:
            continue
        peaks = (data.get(cid) or {}).get("peaks")
        if peaks is None or peaks.empty or "type" not in peaks.columns:
            continue
        cn = (cn_tables or {}).get(cid)
        if cn is not None and "fit_ok" in cn.columns:
            ok = set(cn.loc[pd.to_numeric(cn["fit_ok"], errors="coerce") == 1,
                            "sample"])
            if ok:
                seg = seg[seg["sample"].isin(ok)]
        samples = set(seg["sample"].unique())
        if not samples:
            continue
        calls = data[cid].get("calls")
        del_pairs = set()
        if calls is not None and len(calls):
            cl = calls[pd.to_numeric(calls["call"], errors="coerce") > 0]
            del_pairs = set(zip(cl["peak_id"], cl["sample"]))
        comp = next((f"{c['a']}_vs_{c['b']}" for c in cfg["comparisons"]
                     if cid in (c["a"], c["b"])), cid)
        pk = peaks[(peaks["type"] == "Del")
                   & (pd.to_numeric(peaks["q"], errors="coerce") <= q_main)]
        for _, p in pk.iterrows():
            gc = "peak_" if core == "narrow" else "wide_"
            chrom = str(p.get(f"{gc}chrom", "") or "")
            cs = pd.to_numeric(p.get(f"{gc}start"), errors="coerce")
            ce = pd.to_numeric(p.get(f"{gc}end"), errors="coerce")
            if not chrom or not (np.isfinite(cs) and np.isfinite(ce)) \
                    or ce <= cs:
                continue
            n = dict(del_loh=0, del_only=0, loh_only=0, neither=0,
                     n_loh=0, n_cnloh=0, loh_only_cnloh=0,
                     loh_only_cnloss=0, covered=0)
            for smp, g in seg[seg["chrom"] == chrom].groupby("sample"):
                ov = (np.minimum(g["end"].to_numpy(float), ce)
                      - np.maximum(g["start"].to_numpy(float), cs)).clip(0)
                ln = float(ov.sum())
                if ln <= 0:
                    continue
                mn = g["cn_minor"].to_numpy(float)
                tt = g["cn_total"].to_numpy(float)
                lo = float(ov[mn == 0].sum()) >= cover * ln
                cn_loh = lo and float(ov[(mn == 0) & (tt == 2)].sum()) \
                    >= cover * ln
                cn_loss = lo and float(ov[(mn == 0) & (tt < 2)].sum()) \
                    >= cover * ln
                is_del = (p["peak_id"], smp) in del_pairs
                n["covered"] += 1
                n["n_loh"] += int(lo)
                n["n_cnloh"] += int(cn_loh)
                if lo and is_del:
                    n["del_loh"] += 1
                elif lo:
                    n["loh_only"] += 1
                    n["loh_only_cnloh"] += int(cn_loh)
                    n["loh_only_cnloss"] += int(cn_loss)
                elif is_del:
                    n["del_only"] += 1
                else:
                    n["neither"] += 1
            if not n["covered"]:
                continue
            nsmp = len(samples)
            rows.append(dict(
                comparison=comp, cohort=cid, peak_id=p["peak_id"],
                ptype="Del", cytoband=str(p.get("cytoband", "")).strip(),
                chrom=chrom, start=int(cs), end=int(ce),
                kb=round((float(ce) - float(cs)) / 1e3, 1), core=core,
                gistic_freq=round(float(p["freq_altered"]), 4),
                n_samples=nsmp, n_covered=n["covered"],
                n_loh=n["n_loh"], n_cnloh=n["n_cnloh"],
                del_loh=n["del_loh"], del_only=n["del_only"],
                loh_only=n["loh_only"], neither=n["neither"],
                loh_only_cnloh=n["loh_only_cnloh"],
                loh_only_cnloss=n["loh_only_cnloss"],
                loh_frac_any=round(n["n_loh"] / nsmp, 4),
                cnloh_frac_any=round(n["n_cnloh"] / nsmp, 4),
                loh_only_frac=round(n["loh_only"] / nsmp, 4)))
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values(["comparison", "cohort", "loh_only_frac"],
                           ascending=[True, True, False]).reset_index(drop=True)


def _cross_agg(cross: pd.DataFrame, cid: str, col: str) -> float:
    """Cohort-level totals of the LOH x Del cross (peaks pooled), for the
    one-line-per-comparison summary."""
    if cross is None or getattr(cross, "empty", True):
        return np.nan
    d = cross[cross["cohort"] == cid]
    n = float(d["n_samples"].sum())
    if d.empty or n <= 0:
        return np.nan
    src = "loh_only" if col == "loh_only_frac" else col
    if src not in d.columns:
        return np.nan
    return round(float(d[src].sum()) / n, 4)


# ---------------------------------------------------------------------------
# module entry
# ---------------------------------------------------------------------------

_ABS_METRICS = [
    ("purity", "purity"),
    ("ploidy_abs", "absolute ploidy"),
    ("loh_frac", "LOH fraction"),
    ("cnloh_frac", "copy-neutral LOH"),
    ("wgd_frac", "WGD fraction"),
]
# metrics that additionally get WGD-stratified rows: WGD raises minor CN
# genome-wide and so depresses loh_frac, and the cohort difference in WGD
# composition can dominate (even reverse) the pooled test — see
# docs/stats.md.  Purity/ploidy are not stratified this way: wgd_status is
# defined from them, so within-stratum tests would be conditioned on the
# outcome.
_ABS_STRAT_METRICS = (("loh_frac", "LOH fraction"),
                      ("cnloh_frac", "copy-neutral LOH"))
_ABS_CN_METRICS = ("S_focal_cn_star", "S_arm_cn_star", "S_chrom_cn_star")


def _wgd_class_columns(seg_tables: dict, cn_tables: dict,
                       acfg: dict) -> dict:
    """Per-sample WGD-multiplicity metrics + QC tier, from segment tables.

    strict_frac: autosomal length fraction with nMinor >= 2 — parsimonious
    whole-genome-doubling evidence ((2,1)-type segments are ambiguous: a
    single-allele gain produces them without any WGD).  f_min3p / f_t8p:
    minor >= 3 / total >= 8 fractions — the >=2x rungs of the doubling
    ladder.  qc_tier: purity >= 0.99 fits are non-informative (purity and
    ploidy are unidentifiable: flat "likely normal" genomes or
    balanced-elevated degenerate fits).
    """
    for cid, seg in seg_tables.items():
        cn = cn_tables.get(cid)
        if cn is None or seg is None or seg.empty:
            continue
        sg = seg.copy()
        sg["L"] = (sg["end"] - sg["start"]).clip(lower=1)
        rows = {}
        for s, g in sg.groupby("sample"):
            Ls, L = g["L"], g["L"].sum()
            mn, tt = g["cn_minor"], g["cn_total"]
            rows[s] = dict(
                strict_frac=Ls[mn >= 2].sum() / L,
                f_min3p=Ls[mn >= 3].sum() / L,
                f_t8p=Ls[tt >= 8].sum() / L)
        if not rows:
            continue
        cn = cn.merge(pd.DataFrame.from_dict(rows, orient="index")
                      .rename_axis("sample").reset_index(),
                      on="sample", how="left")
        cn["wgd_class"] = np.where(
            cn["strict_frac"] < acfg["wgd_class_frac"], "0x",
            np.where((cn["f_min3p"] >= acfg["wgd_multi_frac"])
                     | (cn["f_t8p"] >= acfg["wgd_multi_frac"]), ">=2x", "1x"))
        cn["qc_tier"] = np.where(cn["purity"] >= 0.99,
                                 "non-informative", "informative")
        cn_tables[cid] = cn
    return cn_tables


def run_module_abs(cfg: dict[str, Any],
                   data: dict[str, dict]) -> dict[str, Any]:
    rcfg = _regional_cfg(cfg)
    acfg = _abs_cfg(cfg, rcfg)
    out_skip = out_table(cfg, "ABS_skipped.txt")

    def skip(msg: str) -> dict[str, Any]:
        Path(out_skip).write_text(msg + "\n", encoding="utf-8")
        print(f"  [ABS] skipped: {msg}")
        empty = pd.DataFrame()
        return dict(cn_tables={}, scores=empty, wgd=empty, metrics=empty,
                    cn_tests=empty, purity_peaks=empty, summary=empty,
                    skipped=msg)

    if acfg["tool"] != "ascat":
        return skip(f"unknown tool: {acfg['tool']!r} — only 'ascat' is "
                    "implemented (a FACETS count-based path is reserved "
                    "for future sequencing input, not implemented)")

    cn_tables, seg_tables, missing = {}, {}, []
    try:
        cn_tables, seg_tables, missing = _load_abs_tables(cfg, data, acfg)
        if missing and acfg["auto_run"]:
            ok, msg = _run_r_pipeline(cfg, acfg, missing)
            if not ok:
                return skip(msg)
            print(f"  [ABS] R/ASCAT pipeline: {msg}")
            cn_tables, seg_tables, missing = _load_abs_tables(cfg, data, acfg)
            if missing:
                return skip("pipeline ran but tables still missing for "
                            + ", ".join(missing))
    except TableSchemaError as e:
        return skip(str(e))
    if missing:
        return skip(
            f"no absolute tables for {', '.join(missing)} — ASCAT runs on "
            "the server (R/absolute_cn.R, see docs/ascat_handoff.md); place "
            "absolute_cn_<cid>.tsv + absolute_seg_<cid>.tsv next to each "
            "cohort's GISTIC outputs (or set analysis.absolute_cn."
            "tables_dir), or set analysis.absolute_cn.auto_run: true to run "
            "the R wrapper here")

    bands, _ = _load_bands(cfg, rcfg, data)
    if bands.empty:
        return skip("no cytoband table available")
    arms = _arm_intervals(bands)
    cn_tables = _wgd_class_columns(seg_tables, cn_tables, acfg)

    # ---- per-sample absolute metrics + CN-space scores ----
    scores = _cn_scores(cfg, seg_tables, arms, acfg)
    # attach the per-sample metric columns from the harmonized tables
    metric_cols = ["purity", "ploidy_abs", "wgd_status", "wgd_frac",
                   "loh_frac", "cnloh_frac", "flag_low_purity"]
    for col in metric_cols:
        of = {}
        for cid, t in cn_tables.items():
            for r in t.itertuples():
                of[r.sample] = getattr(r, col)
        scores[col] = scores["sample"].map(of)
    low_purity = set(scores.loc[scores["flag_low_purity"].astype(bool),
                                "sample"])

    # ---- per-comparison tests ----
    wgd_rows, metric_rows, cn_rows = [], [], []
    for cmp_ in cfg["comparisons"]:
        ca, cb = cmp_["a"], cmp_["b"]
        if cmp_.get("kind") == "paired":
            print(f"  [ABS] skipped paired comparison {ca}_vs_{cb}")
            continue
        comp = f"{ca}_vs_{cb}"
        pool = scores[scores["cohort"].isin([ca, cb])]

        wg = _wgd_test(pool, ca, cb, comp)
        if wg:
            wgd_rows.append(wg)
            if low_purity:
                sub = pool[~pool["sample"].isin(low_purity)]
                wg2 = _wgd_test(sub, ca, cb, f"{comp} [high-purity only]")
                if wg2:
                    wgd_rows.append(wg2)

        for metric, label in _ABS_METRICS:
            d = pool.dropna(subset=[metric])
            ya = d.loc[d["cohort"] == ca, metric].to_numpy(float)
            yb = d.loc[d["cohort"] == cb, metric].to_numpy(float)
            if len(ya) < 3 or len(yb) < 3:
                continue
            try:
                p_mw = float(stats.mannwhitneyu(ya, yb,
                                                alternative="two-sided")[1])
            except ValueError:
                p_mw = np.nan
            metric_rows.append(dict(
                comparison=comp, metric=metric, label=label, n_a=len(ya),
                n_b=len(yb), p_mw=p_mw, r_rb=round(_rank_biserial(ya, yb), 4),
                median_a=round(float(np.median(ya)), 4),
                median_b=round(float(np.median(yb)), 4),
                cohort_higher=(ca if np.nanmedian(ya) > np.nanmedian(yb)
                               else cb)))

        # LOH within WGD strata: the pooled row above mixes a WGD+ stratum
        # (minor CN >= 1 is common, LOH is suppressed) with a WGD- stratum,
        # so its sign/size follows the cohorts' WGD composition
        for metric, label in _ABS_STRAT_METRICS:
            for sval, slab in ((1, "WGD+ only"), (0, "WGD- only")):
                sub = pool[pool["wgd_status"] == sval]
                sub = sub.dropna(subset=[metric])
                ya = sub.loc[sub["cohort"] == ca, metric].to_numpy(float)
                yb = sub.loc[sub["cohort"] == cb, metric].to_numpy(float)
                if len(ya) < 3 or len(yb) < 3:
                    continue
                try:
                    p_mw = float(stats.mannwhitneyu(
                        ya, yb, alternative="two-sided")[1])
                except ValueError:
                    p_mw = np.nan
                metric_rows.append(dict(
                    comparison=f"{comp} [{slab}]", metric=metric,
                    label=label, n_a=len(ya), n_b=len(yb), p_mw=p_mw,
                    r_rb=round(_rank_biserial(ya, yb), 4),
                    median_a=round(float(np.median(ya)), 4),
                    median_b=round(float(np.median(yb)), 4),
                    cohort_higher=(ca if np.nanmedian(ya) > np.nanmedian(yb)
                                   else cb)))

        for metric in _ABS_CN_METRICS:
            d = pool.dropna(subset=[metric])
            ya = d.loc[d["cohort"] == ca, metric].to_numpy(float)
            yb = d.loc[d["cohort"] == cb, metric].to_numpy(float)
            if len(ya) < 3 or len(yb) < 3:
                continue
            try:
                p_mw = float(stats.mannwhitneyu(ya, yb,
                                                alternative="two-sided")[1])
            except ValueError:
                p_mw = np.nan
            cn_rows.append(dict(
                comparison=comp, metric=metric, n_a=len(ya), n_b=len(yb),
                p_mw=p_mw, r_rb=round(_rank_biserial(ya, yb), 4),
                median_a=round(float(np.median(ya)), 4),
                median_b=round(float(np.median(yb)), 4),
                cohort_higher=(ca if np.nanmedian(ya) > np.nanmedian(yb)
                               else cb)))

    wgd = pd.DataFrame(wgd_rows)
    metrics = pd.DataFrame(metric_rows)
    cn_tests = pd.DataFrame(cn_rows)

    purity_peaks = _purity_adjusted_peaks(cfg, data, cn_tables)
    loh_del = _loh_del_cross(cfg, data, seg_tables, cn_tables, acfg)

    sum_rows = []
    # metric_rows need >=3 samples per cohort (Mann-Whitney); small pilots
    # produce an empty metrics frame without columns
    if len(wgd) and "comparison" in metrics.columns:
        main = wgd[~wgd["comparison"].str.contains("high-purity", na=False)]

        def _mode_of(cid: str, col: str) -> str:
            t = cn_tables.get(cid)
            if t is None or col not in t.columns:
                return ""
            v = t[col].dropna()
            return str(v.iloc[0]) if len(v) else ""

        for r in main.itertuples():
            ca, cb = r.comparison.split("_vs_")
            mp = metrics[(metrics["comparison"] == r.comparison)]
            sum_rows.append(dict(
                comparison=r.comparison, wgd_frac_a=r.wgd_frac_a,
                wgd_frac_b=r.wgd_frac_b, wgd_or=r.odds_ratio,
                wgd_p=r.p_fisher,
                purity_p=(mp.loc[mp["metric"] == "purity", "p_mw"].iloc[0]
                          if len(mp) else np.nan),
                ploidy_p=(mp.loc[mp["metric"] == "ploidy_abs", "p_mw"].iloc[0]
                          if len(mp) else np.nan),
                n_low_purity=len(low_purity),
                loh_only_frac_a=_cross_agg(loh_del, ca, "loh_only_frac"),
                loh_only_frac_b=_cross_agg(loh_del, cb, "loh_only_frac"),
                loh_only_cnloh_frac_a=_cross_agg(loh_del, ca,
                                                 "loh_only_cnloh"),
                loh_only_cnloh_frac_b=_cross_agg(loh_del, cb,
                                                 "loh_only_cnloh"),
                germline_mode_a=_mode_of(ca, "germline_mode"),
                germline_mode_b=_mode_of(cb, "germline_mode"),
                het_source_a=_mode_of(ca, "het_source"),
                het_source_b=_mode_of(cb, "het_source"),
            ))
    summary = pd.DataFrame(sum_rows)

    wgd.to_csv(out_table(cfg, "ABS_wgd.tsv"), sep="\t", index=False)
    metrics.to_csv(out_table(cfg, "ABS_purity_ploidy.tsv"), sep="\t",
                   index=False)
    scores.to_csv(out_table(cfg, "ABS_cn_scores.tsv"), sep="\t", index=False)
    cn_tests.to_csv(out_table(cfg, "ABS_cn_tests.tsv"), sep="\t",
                    index=False)
    purity_peaks.to_csv(out_table(cfg, "ABS_purity_adjusted_peaks.tsv"),
                        sep="\t", index=False)
    loh_del.to_csv(out_table(cfg, "ABS_loh_del_cross.tsv"), sep="\t",
                   index=False)
    summary.to_csv(out_table(cfg, "ABS_summary.tsv"), sep="\t", index=False)
    if out_skip.exists():
        out_skip.unlink()

    if len(wgd):
        for r in wgd[~wgd["comparison"].str.contains(
                "high-purity", na=False)].itertuples():
            print(f"  [ABS] {r.comparison}: WGD+ {r.wgd_frac_a:.0%} vs "
                  f"{r.wgd_frac_b:.0%}  OR={r.odds_ratio} "
                  f"p={r.p_fisher:.3g}")
        if low_purity:
            print(f"  [ABS] {len(low_purity)} low-purity sample(s) excluded "
                  "from the sensitivity rerun")
    return dict(cn_tables=cn_tables, scores=scores, wgd=wgd, metrics=metrics,
                cn_tests=cn_tests, purity_peaks=purity_peaks,
                loh_del=loh_del, summary=summary, skipped=None, acfg=acfg)
