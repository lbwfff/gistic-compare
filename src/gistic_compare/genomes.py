"""Reference-genome constants (hg19).

Currently only hg19 (the GISTIC runs this package was developed against) is
supported.  Keep everything chromosome-related in this single module so a
future hg38 addition is a one-place change.
"""

from __future__ import annotations

from typing import Any

# hg19 autosome lengths (bp), chr1-22.  Mitochondrial and alt contigs are
# excluded on purpose: GISTIC analyses here cover the autosomes only.
_HG19_CHROM_SIZES: dict[str, int] = {
    "chr1": 249250621, "chr2": 243199373, "chr3": 198022430,
    "chr4": 191154276, "chr5": 180915260, "chr6": 171115067,
    "chr7": 159138663, "chr8": 146364022, "chr9": 141213431,
    "chr10": 135534747, "chr11": 135006516, "chr12": 133851895,
    "chr13": 115169878, "chr14": 107349540, "chr15": 102531392,
    "chr16": 90354753, "chr17": 81195210, "chr18": 78077248,
    "chr19": 59128983, "chr20": 63025520, "chr21": 48129895,
    "chr22": 51304566,
}

SUPPORTED_GENOMES: tuple[str, ...] = ("hg19",)


def genome_build(genome: str | None = "hg19") -> str:
    """Validate a genome build name; only 'hg19' is supported."""
    g = str(genome or "hg19").lower().strip()
    if g not in SUPPORTED_GENOMES:
        raise ValueError(
            f"unsupported genome: {genome!r} — this release only supports "
            f"{', '.join(SUPPORTED_GENOMES)}. Run your GISTIC analysis on hg19, "
            "or set genome: hg19 in the config.")
    return g


def get_genome(genome: str | None = "hg19") -> dict[str, Any]:
    """Return the chromosome constants dict for a genome build."""
    return {"build": genome_build(genome), "chrom_sizes": dict(_HG19_CHROM_SIZES)}


# Stable names re-exported by gistic_compare.config for backward
# compatibility with the original project layout.
HG19_CHROM_SIZES = _HG19_CHROM_SIZES
CHROMS: list[str] = list(HG19_CHROM_SIZES)
CHROM_RANK: dict[str, int] = {c: i for i, c in enumerate(CHROMS)}

GENOME_LENGTH = sum(HG19_CHROM_SIZES.values())
