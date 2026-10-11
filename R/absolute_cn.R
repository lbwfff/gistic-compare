#!/usr/bin/env Rscript
# absolute_cn.R — ASCAT-based absolute copy-number pipeline (module ABS,
# Phase 2 wrapper; see docs/plan_absolute_cn.md §4.3).
#
# Reads a manifest (sample, cohort, logr, baf [, germline_logr,
# germline_baf] [, het_calls]) whose paths are relative to the manifest
# directory, runs ASCAT 3.x on every sample of one cohort, applies the
# WGD heuristics from plan §5, and writes the harmonized tables consumed
# by the Python module `gistic_compare/pairwise_absolute.py`:
#
#   <out-dir>/absolute_cn_<cohort>.tsv    one row per sample
#   <out-dir>/absolute_seg_<cohort>.tsv   one row per allelic segment
#
# Marker files are tab-separated with columns (chr, pos, value); all
# samples of a cohort must share the same chr/pos grid (pre-harmonized,
# e.g. from Birdseed logR/BAF on the SNP6 marker list).
#
# Germline handling (the germline data's only role in ASCAT's array
# pipeline is het-SNP selection — the copy-number fit consumes tumour
# logR/BAF exclusively):
#   matched     germline_logr/germline_baf columns (measured BAF);
#               takes precedence when both are present
#   het_calls   the samples' OWN Birdseed genotype calls (AA/AB/BB, or
#               0/1/2): SNP6 arrays measure the germline genotype at
#               most loci, so the tumour array itself is a valid
#               germline-genotype proxy — passed to ascat.aspcf via
#               ascat.gg.  Tumour-only caveat (plan §9.1): loci with
#               somatic LOH look homozygous in the tumour calls and drop
#               out of the het set, so loh_frac/cnloh_frac are
#               lower-bound estimates (recorded in germline_mode).
#   inferred    no germline and no het_calls: ASCAT's built-in
#               predictGermlineGenotypes(platform="AffySNP6")
#
# Per-sample failure tolerance (plan §4.3): a failing sample is reported
# in the per-sample table (fit_ok=0, note) instead of interrupting the
# cohort.
#
# WGD heuristics (plan §5, thresholds are CLI flags):
#   A (primary): autosomal length fraction with minor CN >= 1 AND
#                total CN >= --cn-thr-wgd  >= --loh-frac-wgd  -> WGD+
#   B (reference): ploidy >= --ploidy-wgd (classical ploidy rule)
#   wgd_status is binary (0/1); a single snapshot cannot count multiple
#   WGD rounds reliably (R05 §4.6 caveat).  wgd_disagreement=1 marks
#   samples where A and B disagree.
#
# Usage:
#   Rscript absolute_cn.R --manifest m.tsv --cohort cue --out-dir /path \
#     [--cn-thr-wgd 3.0] [--loh-frac-wgd 0.30] [--min-purity 0.30] \
#     [--ploidy-wgd 3.7] [--penalty 70] [--gamma 0.55] [--seed 42]

suppressWarnings(suppressMessages(library(ASCAT)))

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[i + 1])
  default
}
num <- function(flag, default) as.numeric(getarg(flag, as.character(default)))

manifest_path <- getarg("--manifest")
cohort <- getarg("--cohort")
out_dir <- getarg("--out-dir", ".")
cn_thr_wgd <- num("--cn-thr-wgd", 3.0)
loh_frac_wgd <- num("--loh-frac-wgd", 0.30)
min_purity <- num("--min-purity", 0.30)
ploidy_wgd <- num("--ploidy-wgd", 3.7)
penalty <- num("--penalty", 70)
gamma <- num("--gamma", 0.55)
seed <- as.integer(getarg("--seed", "42"))

if (is.null(manifest_path) || is.null(cohort)) {
  stop("usage: --manifest <tsv> --cohort <id> --out-dir <dir> [flags]")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

man <- read.delim(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
need <- c("sample", "cohort", "logr", "baf")
missing_cols <- setdiff(need, names(man))
if (length(missing_cols)) {
  stop("manifest missing columns: ", paste(missing_cols, collapse = ", "))
}
man$cohort <- tolower(trimws(as.character(man$cohort)))
man <- man[man$cohort == tolower(cohort), ]
if (nrow(man) == 0) stop("no manifest rows for cohort ", cohort)
base_dir <- dirname(normalizePath(manifest_path))
resolve <- function(p) {
  if (is.null(p) || length(p) == 0 || is.na(p) || p == "") return(NA_character_)
  pth <- file.path(base_dir, p)
  if (file.exists(pth)) normalizePath(pth) else NA_character_
}

cat(sprintf("[abs] cohort %s: %d sample(s)\n", cohort, nrow(man)))

# ---------------------------------------------------------------------------
# per-sample marker files (chr, pos, value) -> ASCAT wide-format TSVs
# ---------------------------------------------------------------------------
read_marker_file <- function(p, colname) {
  d <- read.delim(p, stringsAsFactors = FALSE, check.names = FALSE)
  lc <- tolower(names(d))
  i_chr <- which(lc %in% c("chr", "chrom", "chromosome"))
  i_pos <- which(lc %in% c("pos", "position", "start"))
  i_val <- which(lc %in% c(colname, "value", "logr", "baf"))
  if (!length(i_chr) || !length(i_pos) || !length(i_val)) {
    stop(sprintf("file %s lacks chr/pos/%s columns (got: %s)",
                 basename(p), colname, paste(names(d), collapse = ",")))
  }
  d[, c(i_chr[1], i_pos[1], i_val[1])]
}

tumor_logr <- lapply(man$logr, function(p) read_marker_file(resolve(p), "logr"))
tumor_baf <- lapply(man$baf, function(p) read_marker_file(resolve(p), "baf"))
grid_ok <- vapply(seq_len(nrow(man)), function(i)
  nrow(tumor_logr[[i]]) > 0 && nrow(tumor_baf[[i]]) > 0, TRUE)
if (!all(grid_ok)) {
  stop("empty marker files for: ", paste(man$sample[!grid_ok], collapse = ", "))
}
grid <- tumor_logr[[1]][, 1:2]
for (i in seq_along(tumor_logr)) {
  g <- tumor_logr[[i]][, 1:2]
  if (!identical(dim(g), dim(grid)) ||
      !isTRUE(all.equal(as.data.frame(g), as.data.frame(grid)))) {
    stop("marker grids differ across samples (first mismatch near sample ",
         man$sample[i], ") — pre-harmonize marker files to a common grid")
  }
}

work <- file.path(out_dir, "ASCAT_input")
dir.create(work, recursive = TRUE, showWarnings = FALSE)
write_wide <- function(frames, path) {
  marker_ids <- sprintf("M%d", seq_len(nrow(grid)))
  vals <- do.call(cbind, frames)
  wide <- cbind(marker_ids, grid, vals)
  names(wide) <- c("", "chrs", "pos", man$sample)
  write.table(wide, path, sep = "\t", quote = FALSE, row.names = FALSE,
              na = "NA")
}
write_wide(lapply(tumor_logr, `[[`, 3), file.path(work, "Tumor_LogR.txt"))
write_wide(lapply(tumor_baf, `[[`, 3), file.path(work, "Tumor_BAF.txt"))

# ---------------------------------------------------------------------------
# germline mode
#   matched      germline_logr + germline_baf columns (measured BAF)
#   tumor_only   the samples' OWN Birdseed genotype calls (het_calls
#                columns; SNP6 arrays measure the germline genotype at
#                most loci) — passed to ascat.aspcf via ascat.gg
#   tumor_only   last resort: ASCAT's built-in genotype inference
#                (predictGermlineGenotypes, AffySNP6 platform priors)
# A matched germline takes precedence over het_calls; both partially
# present (some samples only) is an error.
# ---------------------------------------------------------------------------
has_het <- "het_calls" %in% names(man) &&
  all(!is.na(man$het_calls) & man$het_calls != "")
has_germline <- all(c("germline_logr", "germline_baf") %in% names(man)) &&
  all(!is.na(man$germline_logr)) && all(!is.na(man$germline_baf))
if (has_germline && has_het) {
  warning("manifest has both germline and het_calls columns — the matched ",
          "germline takes precedence; het_calls ignored")
  has_het <- FALSE
}
germline_mode <- "tumor_only"
het_source <- "inferred"
if (has_germline) {
  germline_mode <- "matched"
  het_source <- "germline_baf"
} else if (has_het) {
  het_source <- "birdseed_calls"
}

if (has_germline) {
  germline_logr <- lapply(man$germline_logr,
                          function(p) read_marker_file(resolve(p), "logr"))
  germline_baf <- lapply(man$germline_baf,
                         function(p) read_marker_file(resolve(p), "baf"))
  write_wide(lapply(germline_logr, `[[`, 3),
             file.path(work, "Germline_LogR.txt"))
  write_wide(lapply(germline_baf, `[[`, 3),
             file.path(work, "Germline_BAF.txt"))
  cat(sprintf("[abs] germline matched for all %d sample(s)\n", nrow(man)))
} else if (has_het) {
  cat(sprintf("[abs] tumor-only mode: het SNPs from the samples' own ",
              "Birdseed calls (%d sample(s))\n", nrow(man)))
} else {
  cat("[abs] no germline columns and no het_calls -> ASCAT built-in ",
      "tumor-only genotype inference\n")
}

ascat.data <- ascat.loadData(
  Tumor_LogR_file = file.path(work, "Tumor_LogR.txt"),
  Tumor_BAF_file = file.path(work, "Tumor_BAF.txt"),
  if (has_germline) file.path(work, "Germline_LogR.txt") else NULL,
  if (has_germline) file.path(work, "Germline_BAF.txt") else NULL,
  chrs = c(1:22))

ascat.gg <- NULL
if (!has_germline && has_het) {
  # per-sample Birdseed call files (chr, pos, call): AA/BB (or 0/2, hom)
  # = germline homozygous, AB (or 1, het) = heterozygous; NA = unknown
  read_call_file <- function(p) {
    d <- read.delim(p, stringsAsFactors = FALSE, check.names = FALSE)
    lc <- tolower(names(d))
    i_chr <- which(lc %in% c("chr", "chrom", "chromosome"))
    i_pos <- which(lc %in% c("pos", "position", "start"))
    i_val <- which(lc %in% c("call", "genotype", "gt", "het_calls", "value"))
    if (!length(i_chr) || !length(i_pos) || !length(i_val)) {
      stop(sprintf("call file %s lacks chr/pos/call columns (got: %s)",
                   basename(p), paste(names(d), collapse = ",")))
    }
    d[, c(i_chr[1], i_pos[1], i_val[1])]
  }
  calls <- lapply(man$het_calls, function(p) read_call_file(resolve(p)))
  for (i in seq_along(calls)) {
    g <- calls[[i]][, 1:2]
    if (!identical(dim(g), dim(grid)) ||
        !isTRUE(all.equal(as.data.frame(g), as.data.frame(grid)))) {
      stop("call grids differ from the marker grid (first mismatch near ",
           "sample ", man$sample[i], ") — het_calls must be on the same ",
           "marker list as logr/baf")
    }
  }
  raw <- do.call(cbind, lapply(calls, function(d)
    tolower(as.character(d[[3]]))))
  het_mat <- matrix(NA, nrow = nrow(grid), ncol = length(calls),
                    dimnames = list(sprintf("M%d", seq_len(nrow(grid))),
                                    man$sample))
  het_mat[raw %in% c("ab", "1", "het")] <- FALSE
  het_mat[raw %in% c("aa", "bb", "0", "2", "hom", "homz")] <- TRUE
  ascat.gg <- list(germlinegenotypes = !het_mat)
}

if (has_germline) {
  ascat.data <- ascat.aspcf(ascat.data, penalty = penalty, seed = seed,
                            out.dir = tempdir())
} else if (!is.null(ascat.gg)) {
  ascat.data <- ascat.aspcf(ascat.data, ascat.gg = ascat.gg,
                            penalty = penalty, seed = seed,
                            out.dir = tempdir())
} else {
  gg <- ascat.predictGermlineGenotypes(ascat.data, platform = "AffySNP6",
                                       img.dir = tempdir())
  ascat.data <- ascat.aspcf(ascat.data, ascat.gg = gg, penalty = penalty,
                            seed = seed, out.dir = tempdir())
}
ascat.output <- ascat.runAscat(ascat.data, gamma = gamma,
                               img.dir = tempdir())

# ---------------------------------------------------------------------------
# per-sample harmonized metrics + WGD heuristics
# ---------------------------------------------------------------------------
segs_all <- ascat.output$segments
if (!is.null(segs_all) && nrow(segs_all) > 0) {
  segs_all$chr <- as.character(segs_all$chr)
} else {
  segs_all <- data.frame(sample = character(), chr = character(),
                         startpos = numeric(), endpos = numeric(),
                         nMajor = numeric(), nMinor = numeric())
}
purity_of <- ascat.output$purity
ploidy_of <- ascat.output$ploidy
nonaberrant <- unique(as.character(ascat.output$nonaberrantarrays))
failed <- unique(as.character(ascat.output$failedarrays))

autos <- function(seg) seg[seg$chr %in% as.character(1:22), , drop = FALSE]
wgd_metrics <- function(seg) {
  seg <- autos(seg)
  if (nrow(seg) == 0) {
    return(list(wgd_frac = NA_real_, loh_frac = NA_real_,
                cnloh_frac = NA_real_, n_seg = 0L))
  }
  L <- sum(pmax(seg$endpos - seg$startpos, 1))
  total <- seg$nMajor + seg$nMinor
  dup <- seg$nMinor >= 1 & total >= cn_thr_wgd
  loh <- seg$nMinor == 0 & total > 0
  cnloh <- seg$nMinor == 0 & total == 2
  list(wgd_frac = sum(pmax(seg$endpos - seg$startpos, 1)[dup]) / L,
       loh_frac = sum(pmax(seg$endpos - seg$startpos, 1)[loh]) / L,
       cnloh_frac = sum(pmax(seg$endpos - seg$startpos, 1)[cnloh]) / L,
       n_seg = nrow(seg))
}

cn_rows <- list()
seg_rows <- list()
for (i in seq_len(nrow(man))) {
  s <- man$sample[i]
  purity <- if (s %in% names(purity_of)) as.numeric(purity_of[[s]]) else NA_real_
  ploidy <- if (s %in% names(ploidy_of)) as.numeric(ploidy_of[[s]]) else NA_real_
  seg <- segs_all[segs_all$sample == s, , drop = FALSE]
  ok <- (s %in% failed) == FALSE && nrow(seg) > 0 &&
    is.finite(purity) && is.finite(ploidy)
  note <- if (s %in% nonaberrant) "non-aberrant" else ""
  if (s %in% failed) note <- paste(note, "failedarray")
  if (ok) {
    m <- wgd_metrics(seg)
    seg_a <- autos(seg)
    seg_rows[[length(seg_rows) + 1]] <- data.frame(
      sample = s, chrom = paste0("chr", seg_a$chr),
      start = as.integer(seg_a$startpos), end = as.integer(seg_a$endpos),
      cn_total = as.integer(seg_a$nMajor + seg_a$nMinor),
      cn_major = as.integer(seg_a$nMajor),
      cn_minor = as.integer(seg_a$nMinor),
      stringsAsFactors = FALSE)
    wgd_a <- !is.na(m$wgd_frac) && m$wgd_frac >= loh_frac_wgd
    wgd_b <- !is.na(ploidy) && ploidy >= ploidy_wgd
    cn_rows[[length(cn_rows) + 1]] <- data.frame(
      sample = s,
      purity = round(purity, 4),
      ploidy_abs = round(ploidy, 3),
      wgd_status = as.integer(wgd_a),
      wgd_frac = round(m$wgd_frac, 4),
      loh_frac = round(m$loh_frac, 4),
      cnloh_frac = round(m$cnloh_frac, 4),
      n_segments_allelic = m$n_seg,
      flag_low_purity = (!is.na(purity) && purity < min_purity),
      wgd_flag_b = as.integer(wgd_b),
      wgd_disagreement = as.integer(wgd_a != wgd_b),
      non_aberrant = as.integer(note == "non-aberrant"),
      germline_mode = germline_mode, het_source = het_source,
      fit_ok = 1L, note = note, stringsAsFactors = FALSE)
  } else {
    if (note == "") note <- "no segments or fit failed"
    cn_rows[[length(cn_rows) + 1]] <- data.frame(
      sample = s, purity = NA_real_, ploidy_abs = NA_real_,
      wgd_status = NA_integer_, wgd_frac = NA_real_, loh_frac = NA_real_,
      cnloh_frac = NA_real_, n_segments_allelic = 0L,
      flag_low_purity = NA, wgd_flag_b = NA_integer_,
      wgd_disagreement = NA_integer_, non_aberrant = NA_integer_,
      germline_mode = germline_mode, het_source = het_source,
      fit_ok = 0L, note = trimws(note), stringsAsFactors = FALSE)
    cat(sprintf("[abs] WARNING sample %s failed: %s\n", s, note))
  }
}

cn <- do.call(rbind, cn_rows)
seg <- if (length(seg_rows)) do.call(rbind, seg_rows) else
  data.frame(sample = character(), chrom = character(), start = integer(),
             end = integer(), cn_total = integer(), cn_major = integer(),
             cn_minor = integer())

p_cn <- file.path(out_dir, sprintf("absolute_cn_%s.tsv", cohort))
p_seg <- file.path(out_dir, sprintf("absolute_seg_%s.tsv", cohort))
write.table(cn, p_cn, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
if (nrow(seg)) {
  write.table(seg, p_seg, sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  writeLines(c("sample\tchrom\tstart\tend\tcn_total\tcn_major\tcn_minor"),
             p_seg)
}
n_ok <- sum(cn$fit_ok == 1)
n_wgd <- sum(cn$wgd_status == 1, na.rm = TRUE)
n_dis <- sum(cn$wgd_disagreement == 1, na.rm = TRUE)
n_low <- sum(cn$flag_low_purity %in% TRUE, na.rm = TRUE)
cat(sprintf("[abs] wrote %s (%d/%d fit ok; WGD+ %d; disagreement %d; "
            , p_cn, n_ok, nrow(cn), n_wgd, n_dis))
cat(sprintf("low purity %d)\n", n_low))
cat(sprintf("[abs] wrote %s (%d segments)\n", p_seg, nrow(seg)))
