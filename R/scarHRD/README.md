# scarHRD (vendored)

Source: https://github.com/sztup/scarHRD (MIT License, (c) 2020 Zsofia
Sztupinszki).  Vendored unmodified — we drive it by `source()`ing these
files (the `sequenza` import is only needed for the seqz input path,
which we do not use).

Reference: Sztupinszki et al., Nat Commun 2018 (PMCID: PMC6028448) —
migration of the SNP-array HRD measures (Abkevich 2012, Popova 2012) to
seg inputs; threshold 42 from Telli et al. 2016, Clin Cancer Res.

Used by `R/hrd_scores.R` to compute per-sample HRD-LOH / TAI / LST /
HRD-sum from the ABS allele-specific segment tables
(`absolute_seg_<cohort>.tsv`).
