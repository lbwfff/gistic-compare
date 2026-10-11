"""Parsers for GISTIC 2.0 output files.

All parsers return tidy pandas DataFrames with stable column names.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .config import CHROM_RANK, CHROMS

# ---------------------------------------------------------------------------
# coordinate helpers
# ---------------------------------------------------------------------------

_COORD_RE = re.compile(
    r"^(?P<chrom>chr[0-9XYM]+):(?P<start>\d+)-(?P<end>\d+)"
)
_PROBES_RE = re.compile(r"\(probes\s+\d+:\d+\)")


def parse_coord(s: str) -> tuple[str, int, int]:
    """Parse 'chr1:13884087-15253589(probes 5947:7043)' -> ('chr1', start, end)."""
    s = str(s).strip()
    s = _PROBES_RE.sub("", s).strip()
    m = _COORD_RE.match(s)
    if not m:
        raise ValueError(f"cannot parse coordinate: {s!r}")
    return m.group("chrom"), int(m.group("start")), int(m.group("end"))


def normalize_chrom(c: Any) -> str | None:
    """Normalize chromosome label to 'chrN'; return None if not autosomal/sex."""
    if c is None or (isinstance(c, float) and np.isnan(c)):
        return None
    s = str(c).strip()
    if not s:
        return None
    if not s.startswith("chr"):
        s = "chr" + s
    # strip .0 from numeric
    if s.endswith(".0") and s[3:-2].isdigit():
        s = "chr" + s[3:-2]
    if s in ("chrM", "chrMT"):
        return None
    return s


def autosomal_only(df: pd.DataFrame, chrom_col: str = "chrom") -> pd.DataFrame:
    return df[df[chrom_col].isin(CHROMS)].copy()


def sort_peaks(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["_r"] = out["chrom"].map(CHROM_RANK)
    out = out.sort_values(["type", "_r", "start", "end"], kind="mergesort")
    return out.drop(columns="_r").reset_index(drop=True)


# ---------------------------------------------------------------------------
# all_lesions.conf_90.txt
# ---------------------------------------------------------------------------

def load_all_lesions(path: str | Path, cohort: str,
                     return_cn: bool = False):
    """Parse all_lesions file.

    GISTIC writes each peak twice: a summary row with 0/1/2 calls and a
    ``- CN values`` row with continuous per-sample CN. Only summary rows are
    used for peaks/calls; CN-value rows are returned separately when
    ``return_cn`` is true (used by module P7).

    Returns
    -------
    peaks : DataFrame
        one row per peak: cohort, peak_id, type(Amp/Del), cytoband,
        wide_peak/peak/region (chrom,start,end), q, residual_q, n_altered_samples
    calls : DataFrame
        long format per-sample calls: peak_id, sample, call(0/1/2)
    cn : DataFrame  (only when return_cn=True)
        long format per-sample continuous CN: peak_id, sample, cn
    """
    path = Path(path)
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    df.columns = [c.strip() for c in df.columns]

    # sample columns = everything after the 9 metadata cols
    meta_cols = [
        "Unique Name", "Descriptor", "Wide Peak Limits", "Peak Limits",
        "Region Limits", "q values",
        "Residual q values after removing segments shared with higher peaks",
        "Broad or Focal", "Amplitude Threshold",
    ]
    # tolerate minor header variations
    sample_cols = [c for c in df.columns if c not in meta_cols]
    # drop accidental empty column names
    sample_cols = [c for c in sample_cols if c.strip() != ""]

    rows = []
    call_rows = []
    cn_rows = []
    uname_to_pid: dict[str, str] = {}
    n_skipped_cn = 0
    for i, r in df.iterrows():
        uname = str(r["Unique Name"]).strip()
        if "CN values" in uname:
            n_skipped_cn += 1
            if return_cn:
                base = uname.replace("- CN values", "").strip()
                pid = uname_to_pid.get(base)
                if pid is None:
                    continue
                for sc in sample_cols:
                    v = str(r[sc]).strip()
                    if v == "":
                        continue
                    try:
                        cn = float(v)
                    except ValueError:
                        continue
                    if np.isfinite(cn):
                        cn_rows.append((pid, sc, cn))
            continue
        if uname.lower().startswith("amplification"):
            ptype = "Amp"
        elif uname.lower().startswith("deletion"):
            ptype = "Del"
        else:
            continue

        wide_c, wide_s, wide_e = parse_coord(r["Wide Peak Limits"])
        pk_c, pk_s, pk_e = parse_coord(r["Peak Limits"])
        rg_c, rg_s, rg_e = parse_coord(r["Region Limits"])

        q = float(r["q values"])
        try:
            rq = float(r["Residual q values after removing segments shared with higher peaks"])
        except (TypeError, ValueError):
            rq = np.nan

        peak_id = f"{cohort}:{ptype}:{i}"
        cytoband = str(r["Descriptor"]).strip()
        uname_to_pid[uname] = peak_id

        rows.append(dict(
            cohort=cohort, peak_id=peak_id, uname=uname, type=ptype,
            cytoband=cytoband,
            wide_chrom=wide_c, wide_start=wide_s, wide_end=wide_e,
            peak_chrom=pk_c, peak_start=pk_s, peak_end=pk_e,
            region_chrom=rg_c, region_start=rg_s, region_end=rg_e,
            q=q, residual_q=rq,
            wide_len=wide_e - wide_s,
        ))

        for sc in sample_cols:
            v = str(r[sc]).strip()
            if v == "":
                continue
            try:
                call = int(float(v))
            except ValueError:
                continue
            # amplitude-threshold calls are 0/1/2; anything else is corrupt
            if call in (0, 1, 2):
                call_rows.append((peak_id, sc, call))

    if n_skipped_cn:
        print(f"  [io] {cohort}: skipped {n_skipped_cn} 'CN values' rows")
    if not rows:
        raise ValueError(f"no peak rows parsed from {path}")

    peaks = pd.DataFrame(rows)
    peaks["chrom"] = peaks["wide_chrom"]      # convenience alias for main field
    peaks["start"] = peaks["wide_start"]
    peaks["end"] = peaks["wide_end"]
    peaks = sort_peaks(peaks)

    calls = pd.DataFrame(call_rows, columns=["peak_id", "sample", "call"])
    n_alt = (
        calls.assign(altered=(calls["call"] > 0).astype(int))
        .groupby("peak_id", sort=False)["altered"].sum()
    )
    peaks["n_altered_samples"] = peaks["peak_id"].map(n_alt).fillna(0).astype(int)
    n_samples = calls["sample"].nunique()
    peaks["freq_altered"] = peaks["n_altered_samples"] / n_samples if n_samples else np.nan

    if return_cn:
        cn = pd.DataFrame(cn_rows, columns=["peak_id", "sample", "cn"])
        return peaks, calls, cn
    return peaks, calls


# ---------------------------------------------------------------------------
# broad_values_by_arm.txt  /  amp_qv / del_qv (arm-level q-values)
# ---------------------------------------------------------------------------

def _arm_key(arm: str) -> str | None:
    """Normalize an arm label ('1p', 'chr1p', '1 p') -> '1p'; None if invalid."""
    s = str(arm).strip().lower().replace("chr", "").replace(" ", "")
    if len(s) < 2 or not s[0].isdigit():
        return None
    if s[1] not in ("p", "q"):
        return None
    num = s[:len(s) - 1]
    if not num.isdigit():
        return None
    return f"{int(num)}{s[-1]}"


def load_broad_by_arm(path: str | Path, cohort: str) -> pd.DataFrame:
    """Parse broad_values_by_arm.txt (per-sample arm-level broad calls).

    Handles the two value conventions seen in the wild:
      * 0/1/2 — 1 = broad amplification, 2 = broad deletion
      * 0/1/-1 — 1 = broad amplification, -1 = broad deletion
    The convention is auto-detected from the value domain (any negative
    value selects the sign convention).  Rows/columns are auto-oriented.

    Returns a long DataFrame: cohort, sample, arm ('1p'...'22q'), type
    (Amp/Del) — one row per (sample, arm) that carries a broad event.
    """
    path = Path(path)
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    df.columns = [str(c).strip() for c in df.columns]
    first = df.columns[0]
    cols = [c for c in df.columns[1:]]
    col_is_arm = all(_arm_key(c) for c in cols) and len(cols) > 0
    if col_is_arm:
        # orientation A: rows = samples, columns = arms
        arms = [_arm_key(c) for c in cols]
        samples = df[first].astype(str).str.strip()
        mat = df[cols].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    else:
        # orientation B: rows = arms, columns = samples
        arms = [_arm_key(v) for v in df[first]]
        sample_cols = cols
        if not all(a for a in arms):
            return pd.DataFrame(columns=["cohort", "sample", "arm", "type"])
        samples = sample_cols
        mat = (df[sample_cols]
               .apply(pd.to_numeric, errors="coerce").to_numpy(float)).T
    # drop columns whose arm label could not be parsed, keeping mat aligned
    valid = [i for i, a in enumerate(arms) if a]
    if not valid or mat.size == 0:
        return pd.DataFrame(columns=["cohort", "sample", "arm", "type"])
    arms = [arms[i] for i in valid]
    mat = mat[:, valid]

    rows = []
    sign_mode = bool((mat < 0).any().any())
    for ai, arm in enumerate(arms):
        col = mat[:, ai]
        for si, s in enumerate(samples):
            v = col[si]
            if not np.isfinite(v) or v == 0:
                continue
            if sign_mode:
                t = "Amp" if v > 0 else "Del"
            else:
                t = "Amp" if v == 1 else ("Del" if v == 2 else None)
            if t is None:
                continue
            rows.append(dict(cohort=cohort, sample=str(s), arm=arm, type=t))
    return pd.DataFrame(rows, columns=["cohort", "sample", "arm", "type"])


def load_arm_qv(path: str | Path, cohort: str, ptype: str) -> pd.DataFrame:
    """Parse amp_qv / del_qv.conf_XX.txt (arm-level q-values), best effort.

    The exact layout varies between GISTIC versions, so this parser is
    tolerant: it expects one label column (arm names) plus one numeric
    q-value column, in either orientation.  Returns a long DataFrame
    (cohort, arm, q) or an empty frame when the file cannot be interpreted;
    the values are display-only (see docs/stats.md).
    """
    path = Path(path)
    try:
        df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    except Exception:
        return pd.DataFrame(columns=["cohort", "arm", "q"])
    if df.empty:
        return pd.DataFrame(columns=["cohort", "arm", "q"])
    df.columns = [str(c).strip() for c in df.columns]

    def _q_from_cols(sub: pd.DataFrame, label_col: str) -> pd.DataFrame | None:
        for c in sub.columns:
            if c == label_col:
                continue
            v = pd.to_numeric(sub[c], errors="coerce")
            if v.notna().sum() >= max(1, int(0.5 * len(sub))):
                arms = [_arm_key(x) for x in sub[label_col]]
                out = pd.DataFrame(dict(arm=arms, q=v.to_numpy(float)))
                out = out[out["arm"].notna() & out["q"].notna()]
                if len(out):
                    return out
        return None

    out = _q_from_cols(df, df.columns[0])
    if out is None:
        # transposed layout: arm labels in the first data row
        sub = df.T
        sub.columns = sub.iloc[0]
        sub = sub.reset_index()
        out = _q_from_cols(sub, sub.columns[0])
    if out is None:
        print(f"  [io] {cohort}: could not parse arm q-values from "
              f"{path.name} (ignored)")
        return pd.DataFrame(columns=["cohort", "arm", "q"])
    out.insert(0, "cohort", cohort)
    return out.reset_index(drop=True)


def apply_limit_field(peaks: pd.DataFrame, field: str = "wide_peak") -> pd.DataFrame:
    """Set chrom/start/end to the requested limit field."""
    m = {"wide_peak": ("wide_chrom", "wide_start", "wide_end"),
         "peak": ("peak_chrom", "peak_start", "peak_end"),
         "region": ("region_chrom", "region_start", "region_end")}
    if field not in m:
        raise ValueError(f"unknown peak_limit_field: {field}")
    c, s, e = m[field]
    out = peaks.copy()
    out["chrom"] = out[c]
    out["start"] = out[s]
    out["end"] = out[e]
    out["wide_len"] = out["end"] - out["start"]
    return out


# ---------------------------------------------------------------------------
# amp_genes / del_genes.conf_90.txt (transposed)
# ---------------------------------------------------------------------------

def load_genes(path: str | Path, cohort: str, ptype: str) -> pd.DataFrame:
    """Parse amp/del_genes.conf_90.txt.

    Returns long DataFrame: cohort, type, col_index, cytoband, q, residual_q,
    wide_peak(chr:start-end), gene
    """
    path = Path(path)
    raw = pd.read_csv(path, sep="\t", dtype=str, header=None, keep_default_na=False)
    if raw.shape[0] < 5:
        return pd.DataFrame(columns=[
            "cohort", "type", "col_index", "cytoband", "q", "residual_q",
            "wide_peak", "gene"])

    def row_vals(i: int) -> list[str]:
        return [str(v).strip() for v in raw.iloc[i, 1:]]

    cytobands = row_vals(0)
    qvals = row_vals(1)
    rqvals = row_vals(2)
    wide_peaks = row_vals(3)

    # gene rows: row 4 is 'genes in wide peak' (top gene per col);
    # subsequent rows are additional genes aligned by column.
    n_cols = len(cytobands)
    # build per-column gene lists
    gene_cols: list[list[str]] = [[] for _ in range(n_cols)]
    for ri in range(4, raw.shape[0]):
        label = str(raw.iloc[ri, 0]).strip()
        vals = row_vals(ri)
        for ci, v in enumerate(vals[:n_cols]):
            if v:
                gene_cols[ci].append(v)

    rows = []
    for ci in range(n_cols):
        cb = cytobands[ci]
        if not cb:
            continue
        try:
            q = float(qvals[ci])
        except ValueError:
            q = np.nan
        try:
            rq = float(rqvals[ci])
        except ValueError:
            rq = np.nan
        wp = wide_peaks[ci]
        for g in gene_cols[ci]:
            rows.append(dict(cohort=cohort, type=ptype, col_index=ci,
                             cytoband=cb, q=q, residual_q=rq,
                             wide_peak=wp, gene=g))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# scores.gistic
# ---------------------------------------------------------------------------

def load_scores(path: str | Path, cohort: str) -> pd.DataFrame:
    """Parse scores.gistic -> tidy DataFrame with numeric columns."""
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    df.columns = [c.strip() for c in df.columns]
    out = pd.DataFrame({
        "cohort": cohort,
        "type": df["Type"].str.strip(),
        "chrom": df["Chromosome"].map(normalize_chrom),
        "start": pd.to_numeric(df["Start"], errors="coerce").astype("int64"),
        "end": pd.to_numeric(df["End"], errors="coerce").astype("int64"),
        "neg_log10_q": pd.to_numeric(df["-log10(q-value)"], errors="coerce"),
        "g_score": pd.to_numeric(df["G-score"], errors="coerce"),
        "avg_amplitude": pd.to_numeric(df["average amplitude"], errors="coerce"),
        "frequency": pd.to_numeric(df["frequency"], errors="coerce"),
    })
    out = out[out["chrom"].isin(CHROMS)].copy()
    out = out.sort_values(["type", "chrom", "start"], key=lambda s: s.map(CHROM_RANK) if s.name == "chrom" else s)
    return out.reset_index(drop=True)


# ---------------------------------------------------------------------------
# seg files
# ---------------------------------------------------------------------------

def load_seg(path: str | Path, cohort: str) -> pd.DataFrame:
    """Parse tumors.seg (GSE has Num_Probes, TCGA does not)."""
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    df.columns = [c.strip() for c in df.columns]
    has_probes = "Num_Probes" in df.columns
    cols = {"Sample": "sample", "Chromosome": "chrom",
            "Start": "start", "End": "end", "Segment_Mean": "seg_mean"}
    if has_probes:
        cols["Num_Probes"] = "num_probes"
    out = pd.DataFrame({v: df[k] for k, v in cols.items()})
    out["cohort"] = cohort
    out["chrom"] = out["chrom"].map(normalize_chrom)
    for c in ("start", "end"):
        out[c] = pd.to_numeric(out[c], errors="coerce").astype("int64")
    out["seg_mean"] = pd.to_numeric(out["seg_mean"], errors="coerce")
    if has_probes:
        out["num_probes"] = pd.to_numeric(out["num_probes"], errors="coerce")
    else:
        out["num_probes"] = np.nan
    out = out[out["chrom"].isin(CHROMS)].copy()
    return out.reset_index(drop=True)


def load_sample_seg_counts(path: str | Path, cohort: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    df.columns = [c.strip() for c in df.columns]
    out = pd.DataFrame({
        "cohort": cohort,
        "sample": df["sample"].str.strip(),
        "segment_count": pd.to_numeric(df["segment_count"], errors="coerce"),
        "included": df["included"].str.strip(),
    })
    return out


# ---------------------------------------------------------------------------
# regions_track bed
# ---------------------------------------------------------------------------

def load_regions_bed(path: str | Path, cohort: str) -> pd.DataFrame:
    lines = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            lines.append(line.split("\t"))
    df = pd.DataFrame(lines)
    out = pd.DataFrame({
        "cohort": cohort,
        "chrom": df[0].map(normalize_chrom),
        "start": pd.to_numeric(df[1], errors="coerce").astype("int64"),
        "end": pd.to_numeric(df[2], errors="coerce").astype("int64"),
        "name": df[3] if df.shape[1] > 3 else "",
    })
    out["type"] = np.where(out["name"].str.contains("AP", case=False), "Amp",
                   np.where(out["name"].str.contains("DP", case=False), "Del", "NA"))
    # all_lesions AP=amp peak, DP=del peak — but bed uses Any-AP/Any-DP
    out.loc[out["name"].str.contains("-DP", case=False, na=False), "type"] = "Del"
    out.loc[out["name"].str.contains("-AP", case=False, na=False), "type"] = "Amp"
    out = out[out["chrom"].isin(CHROMS)].copy()
    return out.reset_index(drop=True)


# ---------------------------------------------------------------------------
# high-level loader
# ---------------------------------------------------------------------------

def load_cohort(cfg: dict[str, Any], cid: str) -> dict[str, Any]:
    """Load every GISTIC file for one cohort into a dict of DataFrames."""
    c = cfg["cohorts"][cid]
    d = Path(c["_abs"])
    prefix = c["prefix"]

    def p(suffix: str) -> Path:
        return d / f"{prefix}.{suffix}"

    peaks, calls, peak_cn = load_all_lesions(p("all_lesions.conf_90.txt"),
                                             cid, return_cn=True)
    limit_field = cfg["analysis"].get("peak_limit_field", "wide_peak")
    peaks = apply_limit_field(peaks, limit_field)

    genes_amp = load_genes(p("amp_genes.conf_90.txt"), cid, "Amp")
    genes_del = load_genes(p("del_genes.conf_90.txt"), cid, "Del")
    genes = pd.concat([genes_amp, genes_del], ignore_index=True)

    scores = load_scores(p("scores.gistic"), cid)
    seg = load_seg(p("tumors.seg"), cid)
    seg_counts = load_sample_seg_counts(p("sample_seg_counts.txt"), cid)
    bed = load_regions_bed(p("regions_track.conf_90.bed"), cid)

    # optional arm-level outputs (module P6); absent in many GISTIC runs
    broad = None
    broad_path = p("broad_values_by_arm.txt")
    if broad_path.exists():
        broad = load_broad_by_arm(broad_path, cid)
    arm_qv = {}
    for ptype, suffix in (("Amp", "amp_qv.conf_90.txt"),
                          ("Del", "del_qv.conf_90.txt")):
        qp = p(suffix)
        if qp.exists():
            arm_qv[ptype] = load_arm_qv(qp, cid, ptype)

    return dict(id=cid, label=c["label"], peaks=peaks, calls=calls,
                peak_cn=peak_cn, genes=genes, scores=scores, seg=seg,
                seg_counts=seg_counts, bed=bed, broad=broad, arm_qv=arm_qv)


def load_all(cfg: dict[str, Any]) -> dict[str, dict[str, Any]]:
    data = {cid: load_cohort(cfg, cid) for cid in cfg["cohorts"]}
    # expose per-cohort sample counts so figures can annotate n without
    # threading the full data dict through every plotting call
    cfg["_n_samples"] = {
        cid: int(d["seg"]["sample"].nunique()) for cid, d in data.items()
    }
    return data
