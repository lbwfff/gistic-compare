"""Confounder sensitivity analysis for the pairwise workflow (module ``CF``).

This module is deliberately **independent** of P4/P5: it only *reads* their
saved tables and never changes them. It answers a bounded question:

    Given only per-cohort marginal distributions of confounders (sex, age,
    smoking, driver mutations), how much of each observed between-cohort
    difference could those confounders explain?

Because no sample-level confounder labels exist, the effect is not
identifiable: without the confounder->outcome association inside each cohort
you cannot *adjust*. What you *can* do is bound the bias. For a binary
confounder with cohort exposure difference ``dp`` the confounding bias is at
most ``dp * delta`` where ``delta`` is the (unknown) confounder-outcome risk
difference, and the minimum ``delta`` needed to fully explain an observed
difference ``D`` is ``delta* = D / dp`` (an E-value in the confounder's own
scale). For age we use ``delta* = D / dMedian`` per year (the agreed
``dMedian ~= dMean`` approximation).

Outputs (all prefixed ``CF_``):
  CF_envelope.tsv    per comparison x confounder: exposure difference dp
  CF_sensitivity.tsv significant/paper findings x confounder: delta*, bounds
  CF_evalue_rr.tsv   Tier-1 generic E-value (RR scale) per finding
  CF_evalue_or.tsv   Tier-1 generic E-value (OR scale) per finding
  CF_summary.tsv     aggregate robustness per source x comparison

Caveats (also in README):
  * margins only -> bounding, not adjustment;
  * multi-confounder combination is an additive, conservative approximation
    (interactions and confounder-correlations unknown);
  * age uses median difference only;
  * EGFR/KRAS/TP53 may be mediators or colliders -> adjusting/bounding on them
    risks over-adjustment; they are reported for sensitivity only.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from .config import out_table
from .presets import band_names as preset_band_names

# default assumed confounder->outcome risk differences to sweep
DEFAULT_DELTA_GRID = [0.02, 0.05, 0.10, 0.20]
# per-year age effects to sweep (smaller scale)
DEFAULT_AGE_GRID = [0.001, 0.002, 0.005, 0.010]
# effects above these are deemed implausible -> finding "robust"
DEFAULT_MAX_BINARY = 0.20
DEFAULT_MAX_AGE = 0.01

# loci that CF always reports — resolved from the configured loci preset
# (analysis.regional.loci_set); empty when no preset is set


def _as_fraction(v: float, n: int | None) -> float:
    """Counts (v>1) -> v/n ; already a proportion -> as is."""
    v = float(v)
    if v > 1.0:
        if not n:
            raise ValueError("count given but cohort size unknown")
        return v / n
    return v


def _parse_cohort_confounders(c: dict[str, Any]) -> dict[str, Any]:
    sex = {k: float(v) for k, v in (c.get("sex") or {}).items()}
    n = c.get("n")
    if n is None and sex:
        n = int(sum(sex.values()))
    n = int(n) if n else None

    out: dict[str, dict] = {}
    if sex:
        tot = sum(sex.values())
        out["sex"] = dict(kind="binary", unit="fraction",
                          p_male=(sex.get("male", 0.0) / tot if tot else np.nan),
                          n=n, raw=sex)
    if c.get("age_median") is not None:
        out["age"] = dict(kind="continuous", unit="year",
                          median=float(c["age_median"]), n=n)
    sm = c.get("smoking")
    if sm:
        known = float(sm.get("smoker", 0.0)) + float(sm.get("nonsmoker", 0.0))
        if known > 0:
            out["smoking"] = dict(
                kind="binary", unit="fraction",
                p_smoker=float(sm.get("smoker", 0.0)) / known,
                n_known=known, n=n, raw=sm)
    muts = c.get("mutations") or {}
    for m, v in muts.items():
        out[f"mut_{m}"] = dict(
            kind="binary", unit="fraction",
            p=_as_fraction(v, n), n=n, raw=v)
    return out


def build_envelope(cfg: dict[str, Any]) -> pd.DataFrame:
    """One row per comparison x confounder with the cohort exposure gap dp."""
    conf = cfg.get("confounding") or {}
    cohorts = conf.get("cohorts") or {}
    parsed = {cid: _parse_cohort_confounders(c) for cid, c in cohorts.items()}
    rows = []
    for cmp_ in cfg["comparisons"]:
        ca, cb = cmp_["a"], cmp_["b"]
        pa, pb = parsed.get(ca, {}), parsed.get(cb, {})
        for name in pa:
            if name not in pb:
                continue
            da, db = pa[name], pb[name]
            if da["kind"] == "binary":
                key = next(k for k in ("p_male", "p_smoker", "p") if k in da)
                va, vb = da[key], db[key]
                dp = va - vb
                detail = (f"{ca}={va:.3f} vs {cb}={vb:.3f}")
            else:  # age, continuous
                va, vb = da["median"], db["median"]
                dp = va - vb
                detail = (f"{ca} med={va:g} vs {cb} med={vb:g}")
            rows.append(dict(comparison=f"{ca}_vs_{cb}", confounder=name,
                             kind=da["kind"], unit=da["unit"],
                             value_a=va, value_b=vb, delta_p=dp,
                             n_a=da.get("n"), n_b=db.get("n"),
                             detail=detail))
    return pd.DataFrame(rows)


def _evalue(x: float) -> float:
    """VanderWeele-Ding E-value for a risk/exposure ratio x (>0)."""
    if not np.isfinite(x) or x <= 0:
        return np.nan
    xr = x if x >= 1 else 1.0 / x
    return float(xr + np.sqrt(xr * (xr - 1.0)))


def _rr_or_cells(fa: float, fb: float, na: int, nb: int):
    """RR and OR from two frequencies + cohort sizes (Haldane +0.5)."""
    a = np.clip(fa, 0, 1); b = np.clip(fb, 0, 1)
    k1 = a * na + 0.5; n1 = na + 1.0
    k2 = b * nb + 0.5; n2 = nb + 1.0
    rr = (k1 / n1) / (k2 / n2)
    or_ = (k1 * (n2 - k2)) / ((n1 - k1) * k2)
    return float(rr), float(or_)


def _load_source(cfg: dict[str, Any], source: str) -> pd.DataFrame:
    """Concatenate window- + band-level differential tables for a module."""
    frames = []
    band_p = out_table(cfg, f"{source}_cytoband_diff.tsv")
    win_p = out_table(cfg, f"{source}_regional_diff.tsv")
    if band_p.exists():
        b = pd.read_csv(band_p, sep="\t")
        if not b.empty:
            b["level"] = "band"
            b["band"] = b["band"].astype(str)
            frames.append(b)
    if win_p.exists():
        w = pd.read_csv(win_p, sep="\t")
        if not w.empty:
            w["level"] = "window"
            w["band"] = w.get("paper_locus", pd.Series("", index=w.index))
            w["band"] = w["band"].fillna("").astype(str)
            frames.append(w)
    if not frames:
        return pd.DataFrame()
    keep = ["comparison", "level", "chrom", "start", "end", "band", "n_a", "n_b",
            "freq_amp_a", "freq_amp_b", "freq_del_a", "freq_del_b",
            "diff_amp", "diff_del",
            "p_amp", "p_del", "p_amp_fdr", "p_del_fdr", "p_mw_fdr"]
    df = pd.concat([f[[c for c in keep if c in f.columns]] for f in frames],
                   ignore_index=True)
    return df


def _conf_effect_columns(cfg: dict[str, Any]) -> tuple[list[float], list[float]]:
    conf = cfg.get("confounding") or {}
    grid = list(conf.get("delta_grid") or DEFAULT_DELTA_GRID)
    age_grid = list(conf.get("age_grid") or DEFAULT_AGE_GRID)
    return grid, age_grid


def _sensitivity_rows(cfg: dict[str, Any], source: str,
                      envelope: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per-finding x confounder bounds + Tier1 E-values. Returns (sens, evalue)."""
    conf = cfg.get("confounding") or {}
    grid, age_grid = _conf_effect_columns(cfg)
    max_bin = float(conf.get("plausible_max_binary", DEFAULT_MAX_BINARY))
    max_age = float(conf.get("plausible_max_age", DEFAULT_MAX_AGE))
    df = _load_source(cfg, source)
    if df.empty:
        return pd.DataFrame(), pd.DataFrame()
    # cohort sizes fall back to n_a/n_b already in the table
    sens_rows, ev_rows = [], []
    for cmp_ in cfg["comparisons"]:
        comp = f"{cmp_['a']}_vs_{cmp_['b']}"
        if comp not in set(df["comparison"]):
            continue
        env = envelope[envelope["comparison"] == comp]
        sub = df[df["comparison"] == comp]
        for _, r in sub.iterrows():
            d_amp = float(r["diff_amp"]); d_del = float(r["diff_del"])
            for direction, D, p, fdr, mw in (
                    ("Amp", d_amp, r["p_amp"], r["p_amp_fdr"], r["p_mw_fdr"]),
                    ("Del", d_del, r["p_del"], r["p_del_fdr"], r["p_mw_fdr"])):
                sig = (np.isfinite(fdr) and fdr < 0.05) or \
                      (np.isfinite(mw) and mw < 0.05)
                is_paper = str(r["band"]) in preset_band_names(cfg)
                if not (sig or is_paper):
                    continue
                freq_a = float(r["freq_amp_a"] if direction == "Amp"
                               else r["freq_del_a"])
                freq_b = float(r["freq_amp_b"] if direction == "Amp"
                               else r["freq_del_b"])
                na = int(r["n_a"]); nb = int(r["n_b"])
                rr, or_ = _rr_or_cells(freq_a, freq_b, na, nb)
                ev_rows.append(dict(
                    source=source, comparison=comp, level=r["level"],
                    chrom=r["chrom"], start=int(r["start"]), end=int(r["end"]),
                    band=r["band"], direction=direction, diff=D,
                    significant=bool(sig),
                    p=p, fdr=fdr, mw_fdr=mw, rr=rr, or_=or_,
                    evalue_rr=_evalue(rr),
                    evalue_or=_evalue(float(np.sqrt(or_))) if or_ > 0 else np.nan))
                if abs(D) < 1e-9:
                    continue
                for _, e in env.iterrows():
                    dp = float(e["delta_p"])
                    if e["kind"] == "binary":
                        if abs(dp) < 1e-6:
                            continue
                        req = D / dp
                        robust = abs(req) > max_bin
                        row = dict(
                            source=source, comparison=comp, level=r["level"],
                            chrom=r["chrom"], start=int(r["start"]),
                            end=int(r["end"]), band=r["band"],
                            direction=direction, diff=D, significant=bool(sig),
                            fdr=fdr, mw_fdr=mw,
                            confounder=e["confounder"], delta_p=dp,
                            required_effect=req, required_unit="risk_diff",
                            robust=bool(robust))
                        for g in grid:
                            row[f"bias_e{g:g}"] = abs(dp) * g
                            row[f"residual_e{g:g}"] = abs(D) - abs(dp) * g
                            row[f"explained_e{g:g}"] = abs(dp) * g >= abs(D)
                    else:  # age
                        if abs(dp) < 1e-6:
                            continue
                        req = D / dp
                        robust = abs(req) > max_age
                        row = dict(
                            source=source, comparison=comp, level=r["level"],
                            chrom=r["chrom"], start=int(r["start"]),
                            end=int(r["end"]), band=r["band"],
                            direction=direction, diff=D, significant=bool(sig),
                            fdr=fdr, mw_fdr=mw,
                            confounder=e["confounder"], delta_p=dp,
                            required_effect=req, required_unit="per_year",
                            robust=bool(robust))
                        for g in age_grid:
                            row[f"bias_e{g:g}"] = abs(dp) * g
                            row[f"residual_e{g:g}"] = abs(D) - abs(dp) * g
                            row[f"explained_e{g:g}"] = abs(dp) * g >= abs(D)
                    sens_rows.append(row)
    return pd.DataFrame(sens_rows), pd.DataFrame(ev_rows)


def _summarise(sens: pd.DataFrame, ev: pd.DataFrame,
               envelope: pd.DataFrame) -> pd.DataFrame:
    rows = []
    if not ev.empty:
        keys = ev.groupby(["source", "comparison"]).size()
    else:
        keys = pd.Series(dtype=int)
    for (source, comp), n_find in keys.items():
        evs = ev[(ev["source"] == source) & (ev["comparison"] == comp)]
        n_sig = int(evs["significant"].sum())
        for _, e in envelope[envelope["comparison"] == comp].iterrows():
            s = sens[(sens["source"] == source)
                     & (sens["comparison"] == comp)
                     & (sens["confounder"] == e["confounder"])]
            n_conf = len(s)
            rows.append(dict(
                source=source, comparison=comp, confounder=e["confounder"],
                kind=e["kind"], unit=e["unit"], delta_p=e["delta_p"],
                n_findings=int(n_find), n_significant=n_sig,
                n_with_confounder_bound=n_conf,
                n_robust=int(s["robust"].sum()) if n_conf else 0,
                median_abs_required_effect=(
                    float(np.nanmedian(np.abs(s["required_effect"])))
                    if n_conf else np.nan),
                min_abs_required_effect=(
                    float(np.nanmin(np.abs(s["required_effect"])))
                    if n_conf else np.nan),
            ))
    return pd.DataFrame(rows)


def run_module_confounding(cfg: dict[str, Any],
                           data: dict[str, dict]) -> dict[str, Any]:
    conf = cfg.get("confounding") or {}
    sources = list(conf.get("sources") or ["P4", "P5"])

    envelope = build_envelope(cfg)
    if envelope.empty:
        print("  [CF] no confounder marginals configured -- nothing to do")
        empty = pd.DataFrame()
        return dict(envelope=empty, sensitivity=empty, evalue_rr=empty,
                    evalue_or=empty, summary=empty, sources=sources)

    sens_all, ev_all = [], []
    for source in sources:
        s, e = _sensitivity_rows(cfg, source, envelope)
        if not s.empty:
            sens_all.append(s)
        if not e.empty:
            ev_all.append(e)
    sens = pd.concat(sens_all, ignore_index=True) if sens_all else pd.DataFrame()
    ev = pd.concat(ev_all, ignore_index=True) if ev_all else pd.DataFrame()

    ev_rr = (ev[["source", "comparison", "level", "chrom", "start", "end",
                 "band", "direction", "diff", "significant", "p", "fdr",
                 "mw_fdr", "rr", "evalue_rr"]].copy() if not ev.empty
             else pd.DataFrame())
    ev_or = (ev[["source", "comparison", "level", "chrom", "start", "end",
                 "band", "direction", "diff", "significant", "p", "fdr",
                 "mw_fdr", "or_", "evalue_or"]].copy() if not ev.empty
             else pd.DataFrame())
    summary = _summarise(sens, ev, envelope)

    envelope.to_csv(out_table(cfg, "CF_envelope.tsv"), sep="\t", index=False)
    sens.to_csv(out_table(cfg, "CF_sensitivity.tsv"), sep="\t", index=False)
    ev_rr.to_csv(out_table(cfg, "CF_evalue_rr.tsv"), sep="\t", index=False)
    ev_or.to_csv(out_table(cfg, "CF_evalue_or.tsv"), sep="\t", index=False)
    summary.to_csv(out_table(cfg, "CF_summary.tsv"), sep="\t", index=False)

    print(f"  [CF] envelope={len(envelope)} rows, "
          f"sensitivity={len(sens)}, evalue={len(ev_rr)}")
    return dict(envelope=envelope, sensitivity=sens, evalue_rr=ev_rr,
                evalue_or=ev_or, summary=summary, sources=sources)
