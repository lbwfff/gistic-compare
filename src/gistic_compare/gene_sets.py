"""Module D: driver gene / cytoband set comparisons."""

from __future__ import annotations

import re
from typing import Any

import pandas as pd

from .config import out_table

_WIDE_PEAK_RE = re.compile(r"^chr[0-9XYM]+:(\d+)-(\d+)")


def _clean_gene(g: str, strip_brackets: bool, keep_mir: bool) -> str | None:
    g = g.strip()
    if not g:
        return None
    if strip_brackets:
        g = g.strip("[]")
    if not keep_mir and re.search(r"(?i)mir", g):
        return None
    if not g:
        return None
    return g


def _wide_peak_len(wp: str) -> int | None:
    """Length of 'chrN:start-end'; None if unparsable."""
    m = _WIDE_PEAK_RE.match(str(wp).strip())
    if not m:
        return None
    return int(m.group(2)) - int(m.group(1))


def run_module_d(cfg: dict[str, Any], data: dict[str, dict]) -> dict[str, Any]:
    ana = cfg["analysis"]
    q_main = ana["q_thresholds"][0]
    gconf = ana.get("genes", {})
    strip_br = gconf.get("strip_brackets", True)
    keep_mir = gconf.get("keep_mir", True)
    # asymmetric caps; fall back to legacy max_wide_len_mb if present
    legacy_mb = gconf.get("max_wide_len_mb")
    max_len_mb = {
        "Amp": gconf.get("max_wide_len_mb_amp", legacy_mb),
        "Del": gconf.get("max_wide_len_mb_del", legacy_mb),
    }
    max_len_bp = {
        ptype: (int(mb * 1_000_000) if mb else None)
        for ptype, mb in max_len_mb.items()
    }

    def _filter_sub(sub: pd.DataFrame, ptype: str) -> pd.DataFrame:
        """q < threshold and optional per-type wide-peak length cap."""
        out = sub[sub["q"] < q_main]
        cap = max_len_bp.get(ptype)
        if cap is not None and "wide_peak" in out.columns:
            lengths = out["wide_peak"].map(_wide_peak_len)
            out = out[lengths.notna() & (lengths <= cap)]
        return out

    # ---------- gene sets ----------
    gene_sets: dict[tuple[str, str], set[str]] = {}
    cytoband_sets: dict[tuple[str, str], set[str]] = {}
    gene_records = []
    for cid, d in data.items():
        gdf = d["genes"]
        for ptype in ("Amp", "Del"):
            sub = gdf[gdf["type"] == ptype]
            if sub.empty:
                gene_sets[(cid, ptype)] = set()
                cytoband_sets[(cid, ptype)] = set()
                continue
            sub = _filter_sub(sub, ptype)
            genes = set()
            for g in sub["gene"].unique():
                cg = _clean_gene(str(g), strip_br, keep_mir)
                if cg:
                    genes.add(cg)
            gene_sets[(cid, ptype)] = genes
            cytoband_sets[(cid, ptype)] = set(sub["cytoband"].unique())
            for g in sorted(genes):
                gene_records.append(dict(cohort=cid, type=ptype, gene=g))

    gene_long = pd.DataFrame(gene_records)

    # ---------- pairwise Jaccard of gene / cytoband sets ----------
    jaccard_rows = []
    for cmp_ in cfg["comparisons"]:
        ca, cb = cmp_["a"], cmp_["b"]
        for ptype in ("Amp", "Del"):
            ga, gb = gene_sets[(ca, ptype)], gene_sets[(cb, ptype)]
            ba, bb = cytoband_sets[(ca, ptype)], cytoband_sets[(cb, ptype)]
            for level, sa, sb in (("gene", ga, gb), ("cytoband", ba, bb)):
                inter = len(sa & sb)
                union = len(sa | sb)
                jaccard_rows.append(dict(
                    comparison=f"{ca}_vs_{cb}",
                    kind=cmp_.get("kind", "unpaired"),
                    level=level, ptype=ptype,
                    n_a=len(sa), n_b=len(sb),
                    n_shared=inter,
                    jaccard=inter / union if union else float("nan"),
                    only_a=len(sa - sb), only_b=len(sb - sa),
                ))
    jaccard_df = pd.DataFrame(jaccard_rows)

    # ---------- three-way gene presence matrix ----------
    cohorts = list(data.keys())
    all_genes = sorted(set().union(*gene_sets.values()) if gene_sets else set())
    presence = pd.DataFrame({"gene": all_genes})
    for cid in cohorts:
        presence[f"{cid}_amp"] = presence["gene"].isin(gene_sets.get((cid, "Amp"), set()))
        presence[f"{cid}_del"] = presence["gene"].isin(gene_sets.get((cid, "Del"), set()))
    presence["n_amp"] = presence[[f"{c}_amp" for c in cohorts]].sum(axis=1)
    presence["n_del"] = presence[[f"{c}_del" for c in cohorts]].sum(axis=1)
    presence["shared_amp_all"] = presence["n_amp"] == len(cohorts)
    presence["shared_del_all"] = presence["n_del"] == len(cohorts)

    # ---------- driver matrix ----------
    drivers = cfg["analysis"].get("drivers", []) or []
    drivers = list(dict.fromkeys(drivers))  # dedupe, keep order
    drv = pd.DataFrame({"gene": drivers})
    for cid in cohorts:
        drv[f"{cid}_amp"] = drv["gene"].isin(gene_sets.get((cid, "Amp"), set()))
        drv[f"{cid}_del"] = drv["gene"].isin(gene_sets.get((cid, "Del"), set()))
        # note: gene may appear in both if peak sets overlap types; keep as-is
    drv["n_cohorts_amp"] = drv[[f"{c}_amp" for c in cohorts]].sum(axis=1)
    drv["n_cohorts_del"] = drv[[f"{c}_del" for c in cohorts]].sum(axis=1)
    drv["detected_any"] = drv["n_cohorts_amp"] + drv["n_cohorts_del"] > 0

    # ---------- per-cohort peak → gene table (for reference) ----------
    peak_gene_rows = []
    for cid, d in data.items():
        gdf = d["genes"]
        for ptype in ("Amp", "Del"):
            sub = gdf[gdf["type"] == ptype]
            if sub.empty:
                continue
            sub = _filter_sub(sub, ptype)
            for (cb, q), grp in sub.groupby(["cytoband", "q"]):
                genes = sorted({cg for g in grp["gene"]
                                if (cg := _clean_gene(str(g), strip_br, keep_mir))})
                wp = grp["wide_peak"].iloc[0] if "wide_peak" in grp.columns else ""
                peak_gene_rows.append(dict(cohort=cid, type=ptype,
                                           cytoband=cb, q=q,
                                           wide_peak=wp,
                                           wide_len=_wide_peak_len(wp),
                                           n_genes=len(genes),
                                           genes=",".join(genes)))
    peak_gene = pd.DataFrame(peak_gene_rows)

    # ---------- set lists (long) ----------
    set_rows = []
    for (cid, ptype), genes in gene_sets.items():
        for g in sorted(genes):
            set_rows.append(dict(cohort=cid, type=ptype, level="gene", item=g))
    for (cid, ptype), cbs in cytoband_sets.items():
        for cb in sorted(cbs):
            set_rows.append(dict(cohort=cid, type=ptype, level="cytoband", item=cb))
    set_long = pd.DataFrame(set_rows)

    for name, df in (("gene_jaccard", jaccard_df), ("gene_presence", presence),
                     ("driver_matrix", drv), ("peak_genes", peak_gene),
                     ("gene_sets_long", set_long), ("gene_sets_by_cohort", gene_long)):
        df.to_csv(out_table(cfg, f"D_{name}.tsv"), sep="\t", index=False)

    return dict(jaccard=jaccard_df, presence=presence, drivers=drv,
                peak_genes=peak_gene, sets_long=set_long,
                gene_sets={k: sorted(v) for k, v in gene_sets.items()},
                cytoband_sets={k: sorted(v) for k, v in cytoband_sets.items()})
