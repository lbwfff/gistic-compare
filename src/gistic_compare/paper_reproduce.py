"""Module G: benchmark experiment cohorts against a reference panel.

Only experiment cohorts are benchmarked (the control, if any, is excluded).

The reference values are **not** hard-coded: set

    paper_benchmark:
      enabled: true
      panel: presets/reference_panels/we_lung.yaml

in the config, and this module tests the experiment cohorts against those
values. ``presets/reference_panels/we_lung.yaml`` ships as an example (the WE
lung series used while developing this package); copy it and edit the numbers
for your own reference series.

Cytoband Descriptor is the band of the peak center; the same locus can be
labeled 1q21.3 vs 1q22 across runs (coords nearly identical). For panel
regions we therefore fall back to a genomic window when the prefix misses.

Outputs (results/<run_id>/tables/):
  G_paper_genome_pct.tsv   — amp/del genome % per experiment vs panel
  G_paper_top_freq.tsv     — peak frequency of panel-named top regions
  G_paper_cytoband.tsv     — presence of each panel cytoband/arm per experiment
  G_paper_driver.tsv       — presence of panel-named driver genes per experiment
  G_paper_summary.tsv      — one row per experiment: n metrics recovered
"""

from __future__ import annotations

import re
from typing import Any

import numpy as np
import pandas as pd

from .config import out_table
from .genomes import GENOME_LENGTH

_SKIP_HINT = (
    "[G] paper_benchmark not configured — add\n"
    "    paper_benchmark:\n"
    "      enabled: true\n"
    "      panel: presets/reference_panels/we_lung.yaml\n"
    "  to the config (copy the example and edit the numbers for your study).")


def panel(cfg: dict[str, Any]) -> dict[str, Any]:
    """Return the reference panel attached to `cfg`, or raise a clear error."""
    p = (cfg or {}).get("_paper_panel")
    if not p:
        raise RuntimeError(_SKIP_HINT)
    return p


# ---------------------------------------------------------------------------
# interval helpers
# ---------------------------------------------------------------------------

_WIDE_PEAK_RE = re.compile(r"^chr[0-9XYM]+:(\d+)-(\d+)")


def _union_bp(peaks: pd.DataFrame) -> int:
    """Total bp covered by union of intervals on chr1-22."""
    if peaks.empty:
        return 0
    total = 0
    for _, g in peaks.groupby("chrom"):
        ivs = sorted(zip(g["start"].to_numpy(int), g["end"].to_numpy(int)))
        cur_s, cur_e = ivs[0]
        for s, e in ivs[1:]:
            if s > cur_e:
                total += cur_e - cur_s
                cur_s, cur_e = s, e
            else:
                cur_e = max(cur_e, e)
        total += cur_e - cur_s
    return total


def _wide_peak_len(wp: str) -> int | None:
    m = _WIDE_PEAK_RE.match(str(wp).strip())
    if not m:
        return None
    return int(m.group(2)) - int(m.group(1))


def _match_region(peaks: pd.DataFrame, prefix: str,
                  panel: dict[str, Any]) -> tuple[pd.DataFrame, str]:
    """Match by cytoband prefix; fall back to panel genomic window.

    Returns (hits, method) where method in {"cytoband", "window", "none"}.
    """
    cb = peaks["cytoband"].astype(str).str.strip()
    hits = peaks[cb.str.startswith(prefix)]
    if len(hits):
        return hits, "cytoband"
    win = panel["windows"].get(prefix)
    if win:
        chrom, lo, hi = win
        m = ((peaks["chrom"] == chrom)
             & (peaks["end"] > lo) & (peaks["start"] < hi))
        hits = peaks[m]
        if len(hits):
            return hits, "window"
    return hits, "none"


def _max_freq(hits: pd.DataFrame) -> float:
    if hits.empty:
        return np.nan
    return float(hits["freq_altered"].max()) * 100.0  # percent


def _gene_present(genes: pd.DataFrame, gene: str, q_main: float,
                  max_len_mb: dict[str, float | None]) -> tuple[bool, str]:
    """True if gene appears in amp/del peak gene lists (q + length filtered).

    Returns (present, which) where which in {"Amp", "Del", "both", ""}.
    """
    g = genes[genes["gene"].astype(str).str.strip("[]") == gene]
    if g.empty:
        return False, ""
    found = set()
    for ptype in ("Amp", "Del"):
        sub = g[g["type"] == ptype]
        sub = sub[sub["q"] < q_main]
        mb = max_len_mb.get(ptype)
        if mb is not None and "wide_peak" in sub.columns:
            lengths = sub["wide_peak"].map(_wide_peak_len)
            sub = sub[lengths.notna() & (lengths <= mb * 1_000_000)]
        if len(sub):
            found.add(ptype)
    if found == {"Amp", "Del"}:
        return True, "both"
    if found:
        return True, next(iter(found))
    return False, ""


def run_module_g(cfg: dict[str, Any], data: dict[str, dict]) -> dict[str, Any]:
    panel_ = panel(cfg)
    q_main = cfg["analysis"]["q_thresholds"][0]
    gconf = cfg["analysis"].get("genes", {})
    legacy = gconf.get("max_wide_len_mb")
    max_len_mb = {
        "Amp": gconf.get("max_wide_len_mb_amp", legacy),
        "Del": gconf.get("max_wide_len_mb_del", legacy),
    }

    # only experiments — the control is not deemed panel-like
    control = cfg.get("control")
    experiments = [c for c in data if c != control]
    if not experiments:
        experiments = list(data)

    genome_rows = []
    freq_rows = []
    cyto_rows = []
    driver_rows = []

    amp_list = list(dict.fromkeys(
        panel_["amp_cytobands"] + panel_["broad_amp_arms"] + panel_["focal_amp_arms"]))
    del_list = list(dict.fromkeys(
        panel_["del_cytobands"] + panel_["large_del_arms"]))

    for cid in experiments:
        d = data[cid]
        peaks = d["peaks"]
        sig = peaks[peaks["q"] < q_main]
        amp = sig[sig["type"] == "Amp"]
        delp = sig[sig["type"] == "Del"]

        amp_pct = 100.0 * _union_bp(amp) / GENOME_LENGTH
        del_pct = 100.0 * _union_bp(delp) / GENOME_LENGTH
        genome_rows.append(dict(
            cohort=cid, label=d.get("label", cid), q_thresh=q_main,
            amp_genome_pct=amp_pct, paper_amp_genome_pct=panel_["amp_genome_pct"],
            delta_amp_pct=amp_pct - panel_["amp_genome_pct"],
            amp_match=abs(amp_pct - panel_["amp_genome_pct"]) <= 3.0,
            del_genome_pct=del_pct, paper_del_genome_pct=panel_["del_genome_pct"],
            delta_del_pct=del_pct - panel_["del_genome_pct"],
            del_match=abs(del_pct - panel_["del_genome_pct"]) <= 3.0,
            n_amp_peaks=len(amp), n_del_peaks=len(delp),
        ))

        # top-frequency named regions (prefix + window fallback)
        for entry in panel_["top_amp"] + panel_["top_del"]:
            prefix, min_f = entry[0], entry[1]
            label = entry[2] if len(entry) > 2 else prefix
            side = "Amp" if entry in panel_["top_amp"] else "Del"
            pool = amp if side == "Amp" else delp
            hits, method = _match_region(pool, prefix, panel_)
            f = _max_freq(hits)
            freq_rows.append(dict(
                cohort=cid, q_thresh=q_main, side=side, region=label,
                prefix=prefix, match_method=method,
                paper_min_freq_pct=min_f,
                observed_freq_pct=f, n_peaks=len(hits),
                best_cytoband=(hits.sort_values("q").iloc[0]["cytoband"]
                               if len(hits) else ""),
                freq_ok=bool(np.isfinite(f) and f >= min_f),
            ))

        # presence of each panel-claimed region
        for prefix in amp_list:
            hits, method = _match_region(amp, prefix, panel_)
            kind = ("broad_amp" if prefix in panel_["broad_amp_arms"]
                    else "focal_amp" if prefix in panel_["focal_amp_arms"]
                    else "amp")
            cyto_rows.append(dict(
                cohort=cid, q_thresh=q_main, side="Amp", prefix=prefix,
                paper_class=kind, match_method=method,
                present=not hits.empty, n_peaks=len(hits),
                best_cytoband=(hits.sort_values("q").iloc[0]["cytoband"]
                               if len(hits) else ""),
                best_q=(float(hits["q"].min()) if len(hits) else np.nan),
                max_freq_pct=_max_freq(hits),
            ))
        for prefix in del_list:
            hits, method = _match_region(delp, prefix, panel_)
            kind = "large_del" if prefix in panel_["large_del_arms"] else "del"
            cyto_rows.append(dict(
                cohort=cid, q_thresh=q_main, side="Del", prefix=prefix,
                paper_class=kind, match_method=method,
                present=not hits.empty, n_peaks=len(hits),
                best_cytoband=(hits.sort_values("q").iloc[0]["cytoband"]
                               if len(hits) else ""),
                best_q=(float(hits["q"].min()) if len(hits) else np.nan),
                max_freq_pct=_max_freq(hits),
            ))

        # panel-named driver genes (groups collapse to a single row)
        groups = panel_.get("driver_groups", {})
        for gene in panel_["drivers"]:
            members = groups.get(gene, [gene])
            found: set[str] = set()
            for m in members:
                ok, which = _gene_present(d["genes"], m, q_main, max_len_mb)
                if ok:
                    found.add(which)
            if found:
                has_amp = any("Amp" in w for w in found)
                has_del = any("Del" in w for w in found)
                which = ("both" if (has_amp and has_del)
                         else ("Amp" if has_amp else "Del"))
                present = True
            else:
                present, which = False, ""
            driver_rows.append(dict(
                cohort=cid, q_thresh=q_main, gene=gene,
                present=present, amp_or_del=which,
            ))

    genome = pd.DataFrame(genome_rows)
    freq = pd.DataFrame(freq_rows)
    cyto = pd.DataFrame(cyto_rows)
    driver = pd.DataFrame(driver_rows)

    # per-experiment recovery summary
    sum_rows = []
    for cid in experiments:
        g = genome[genome["cohort"] == cid].iloc[0]
        f = freq[freq["cohort"] == cid]
        c = cyto[cyto["cohort"] == cid]
        dv = driver[driver["cohort"] == cid]
        n_named = len(c)
        n_hit = int(c["present"].sum())
        core_mask = ((c["side"] == "Amp")
                     & c["paper_class"].isin(["amp", "focal_amp"]))
        sum_rows.append(dict(
            cohort=cid,
            amp_pct_ok=bool(g["amp_match"]),
            del_pct_ok=bool(g["del_match"]),
            amp_genome_pct=round(float(g["amp_genome_pct"]), 1),
            del_genome_pct=round(float(g["del_genome_pct"]), 1),
            top_freq_n_ok=int(f["freq_ok"].sum()),
            top_freq_n=len(f),
            named_regions_recovered=f"{n_hit}/{n_named}",
            named_region_recall=n_hit / n_named if n_named else np.nan,
            core_amp_ok=int(c.loc[core_mask, "present"].sum()),
            core_amp_n=int(core_mask.sum()),
            listed_del_ok=int(c.loc[c["side"] == "Del", "present"].sum()),
            listed_del_n=int((c["side"] == "Del").sum()),
            drivers_recovered=f"{int(dv['present'].sum())}/{len(dv)}",
            driver_recall=(float(dv["present"].mean()) if len(dv) else np.nan),
        ))
    summary = pd.DataFrame(sum_rows)

    genome.to_csv(out_table(cfg, "G_paper_genome_pct.tsv"), sep="\t", index=False)
    freq.to_csv(out_table(cfg, "G_paper_top_freq.tsv"), sep="\t", index=False)
    cyto.to_csv(out_table(cfg, "G_paper_cytoband.tsv"), sep="\t", index=False)
    driver.to_csv(out_table(cfg, "G_paper_driver.tsv"), sep="\t", index=False)
    summary.to_csv(out_table(cfg, "G_paper_summary.tsv"), sep="\t", index=False)

    return dict(genome=genome, top_freq=freq, cytoband=cyto,
                driver=driver, summary=summary, panel=panel_)
