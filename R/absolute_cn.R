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
# GC content reference for logR correction (ASCAT 3 standard SNP6 step;
# without it the grid search fails to reach the GoF threshold on
# HAPMAP-normalised logR — verified on the CUE/CHB pilot):
#   --gc-file GC_AffySNP6_102015.txt   (Zenodo record 14008443; keyed by
#   probe name — the markers grid carries the probe ids)
#
# Usage:
#   Rscript absolute_cn.R --manifest m.tsv --cohort cue --out-dir /path \
#     [--cn-thr-wgd 3.0] [--loh-frac-wgd 0.30] [--min-purity 0.30] \
#     [--ploidy-wgd 3.7] [--penalty 70] [--gamma 0.55] [--seed 42] \
#     [--gc-file <GC_AffySNP6_102015.txt>]

suppressWarnings(suppressMessages(library(ASCAT)))
data.table::setDTthreads(1L)  # before any fread: macOS fork safety (see
                              # GSEdata/ASCAT/ASCAT_SNP6_优化报告.md §5.6)

# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[i + 1])
  default
}
num <- function(flag, default) as.numeric(getarg(flag, as.character(default)))

# optional ASCAT fast stack (R/ascat_opt/): bit-identical optimisations of
# the ASCAT 3.2.0 hot paths — Rcpp aspcfpart, vectorised helpers, fread
# loader, sample-level mclapply parallelism.  Source: our own benchmarking
# (GSEdata/ASCAT/ASCAT_SNP6_优化报告.md, identical()-verified).  Falls back
# to stock ASCAT when the files or Rcpp are unavailable.
# ---------------------------------------------------------------------------
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                        value = TRUE)[1])
opt_dir <- getarg("--ascat-opt-dir",
                  file.path(dirname(script_path), "ascat_opt"))
use_fast_ascat <- FALSE
if (all(file.exists(file.path(opt_dir, c("ascat_fast.R", "aspcf_cpp.cpp",
                                         "ascat_io_fast.R",
                                         "ascat_parallel.R")))) &&
    requireNamespace("Rcpp", quietly = TRUE)) {
  tryCatch({
    Rcpp::sourceCpp(file.path(opt_dir, "aspcf_cpp.cpp"))
    source(file.path(opt_dir, "ascat_fast.R"), local = TRUE)
    source(file.path(opt_dir, "ascat_io_fast.R"), local = TRUE)
    source(file.path(opt_dir, "ascat_parallel.R"), local = TRUE)
    patch_ascat()
    use_fast_ascat <- TRUE
    cat("[abs] ASCAT fast stack loaded (patch_ascat)\n")
  }, error = function(e) {
    cat("[abs] ASCAT fast stack failed to load (",
        conditionMessage(e), ") — using stock ASCAT\n")
  })
} else {
  cat("[abs] ASCAT fast stack not found at ", opt_dir,
      " — using stock ASCAT\n", sep = "")
}


manifest_path <- getarg("--manifest")
cohort <- getarg("--cohort")
out_dir <- getarg("--out-dir", ".")
cn_thr_wgd <- num("--cn-thr-wgd", 3.0)
loh_frac_wgd <- num("--loh-frac-wgd", 0.30)
min_purity <- num("--min-purity", 0.30)
ploidy_wgd <- num("--ploidy-wgd", 3.7)
penalty <- num("--penalty", 70)
gamma <- num("--gamma", 0.55)
max_ploidy <- num("--max-ploidy", 5.5)
seed <- as.integer(getarg("--seed", "42"))
gc_file <- getarg("--gc-file")
ncores <- as.integer(getarg("--ncores", "4"))
cat(sprintf("[dbg] parsed max_ploidy=%s gamma=%s\n", max_ploidy, gamma))

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

# BAF placeholder semantics (official ASCAT_fromCELfiles.R): PennCNV marks
# CN probes (no genotype) with BAF=2 — the official script converts those
# to NA before ASCAT.  BAF NA is therefore the EXPECTED state on CN probes
# and ASCAT handles it natively (they simply don't enter BAF segmentation).
for (i in seq_along(tumor_baf)) {
  v <- tumor_baf[[i]][[3]]
  v[!is.na(v) & v == 2] <- NA
  tumor_baf[[i]][[3]] <- v
}

# NA hygiene (ASCAT 优化报告 §5.1): a single NA probe in logR makes getMad()
# return NA for its whole 1000-probe window and fastAspcf then silently
# SKIPS that window (degenerate segmentation, no warning).  Drop rows with
# NA logR here, once, across all samples.  BAF NAs are NOT a reason to drop
# a row (official CN-probe state — see above).
na_keep <- stats::complete.cases(do.call(
  cbind, lapply(tumor_logr, `[[`, 3)))
if (!all(na_keep)) {
  cat(sprintf("[abs] dropping %d marker row(s) with NA logR values\n",
              sum(!na_keep)))
  grid <- grid[na_keep, , drop = FALSE]
  tumor_logr <- lapply(tumor_logr, function(d) d[na_keep, , drop = FALSE])
  tumor_baf <- lapply(tumor_baf, function(d) d[na_keep, , drop = FALSE])
}

work <- file.path(out_dir, "ASCAT_input")
dir.create(work, recursive = TRUE, showWarnings = FALSE)

# probe-name grid — the ASCAT GC-content reference file is keyed by probe
# name, so the wide-format matrices must carry real probe ids as rownames
markers_path <- file.path(dirname(manifest_path), cohort,
                          "markers.tsv.gz")
if (!file.exists(markers_path)) {
  stop("markers grid not found: ", markers_path, " — re-run ",
       "scripts/prepare_abs_inputs.py")
}
markers <- read.delim(markers_path)
if (!all(na_keep)) markers <- markers[na_keep, , drop = FALSE]
if (nrow(markers) != nrow(grid) || any(markers$chr != grid[[1]]) ||
    any(markers$pos != grid[[2]])) {
  stop("markers.tsv.gz grid differs from the per-sample marker grid")
}

# Official SNP6 protocol (ASCAT ExampleData/ASCAT_fromCELfiles.R): PennCNV
# derives CN probes and SNP probes on different baselines, so logR must be
# centred SEPARATELY per probe class, per sample — otherwise the two
# interleaved probe classes form a systematic sawtooth (measured offsets
# +0.01..+0.11 log2 in this cohort) that pollutes the noise estimate, the
# segmentation and the (rho, psi) distance surface.
is_cn_probe <- startsWith(as.character(markers$name), "CN")
for (i in seq_along(tumor_logr)) {
  v <- tumor_logr[[i]][[3]]
  v[is_cn_probe] <- v[is_cn_probe] - mean(v[is_cn_probe], na.rm = TRUE)
  v[!is_cn_probe] <- v[!is_cn_probe] - mean(v[!is_cn_probe], na.rm = TRUE)
  tumor_logr[[i]][[3]] <- v
}
cat(sprintf(paste0("[abs] official SNP6 CN/SNP probe-class centring applied ",
                   "(%d CN / %d SNP probes)\n"),
            sum(is_cn_probe), sum(!is_cn_probe)))

write_wide <- function(frames, path) {
  vals <- do.call(cbind, frames)
  wide <- cbind(markers$name, grid, vals)
  names(wide) <- c("", "chrs", "pos", man$sample)
  write.table(wide, path, sep = "\t", quote = FALSE, row.names = FALSE,
              na = "NA")
}

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

# ---------------------------------------------------------------------------
# two-level RDS cache (content-keyed; the pipeline is deterministic given
# the inputs):  d1 = post loadData + GC correction, d2 = post aspcf.
# A d2 hit skips loadData/GC/aspcf entirely — downstream-parameter
# experiments (gamma, max_ploidy, manual-fit logic) then cost minutes
# instead of a full recompute.
# ---------------------------------------------------------------------------
cache_dir <- getarg("--cache-rds")
use_cache <- !is.null(cache_dir) && nzchar(cache_dir)
d1_hit <- d2_hit <- FALSE
d1_path <- d2_path <- NULL
if (use_cache) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  in_files <- c(markers_path,
                resolve(man$logr), resolve(man$baf),
                if (has_het) resolve(man$het_calls) else character(0),
                if (has_germline) resolve(man$germline_logr) else character(0),
                if (has_germline) resolve(man$germline_baf) else character(0),
                if (!is.null(gc_file) && nzchar(gc_file))
                  resolve(gc_file) else character(0))
  in_files <- in_files[!is.na(in_files)]
  keyfile <- tempfile(fileext = ".key")
  # pipeline version invalidates caches produced before the official SNP6
  # probe-class centring + BAF-placeholder-as-NA semantics were adopted
  writeLines(paste(c("abs-cn", "v3-gg-namekeyed", cohort, penalty, seed,
                     unname(tools::md5sum(in_files)))), keyfile)
  key_hash <- unname(tools::md5sum(keyfile))
  unlink(keyfile)
  d1_path <- file.path(cache_dir, paste0("d1_", key_hash, ".rds"))
  d2_path <- file.path(cache_dir, paste0("d2_", key_hash, ".rds"))
  d1_hit <- file.exists(d1_path)
  d2_hit <- file.exists(d2_path)
  if (d2_hit) {
    cat("[abs] cache hit: d2 (post aspcf) — skipping loadData/GC/aspcf\n")
    ascat.data <- readRDS(d2_path)
  } else if (d1_hit) {
    cat("[abs] cache hit: d1 (post loadData+GC) — skipping loadData/GC\n")
    ascat.data <- readRDS(d1_path)
  }
}

if (!d2_hit) {
write_wide(lapply(tumor_logr, `[[`, 3), file.path(work, "Tumor_LogR.txt"))
write_wide(lapply(tumor_baf, `[[`, 3), file.path(work, "Tumor_BAF.txt"))

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
  cat(paste0("[abs] tumor-only mode: het SNPs from the samples' own ",
             "Birdseed calls (", nrow(man), " sample(s))\n"))
} else {
  cat("[abs] no germline columns and no het_calls -> ASCAT built-in ",
      "tumor-only genotype inference\n")
}

tmark <- function(label, t0) cat(sprintf("[abs][t] %s: %.1f s\n", label,
                                as.numeric(difftime(Sys.time(), t0,
                                                    units = "secs"))))
if (d1_hit) {
  ascat.data <- readRDS(d1_path)
} else {
t0 <- Sys.time()
ascat.data <- ascat.loadData(
  Tumor_LogR_file = file.path(work, "Tumor_LogR.txt"),
  Tumor_BAF_file = file.path(work, "Tumor_BAF.txt"),
  if (has_germline) file.path(work, "Germline_LogR.txt") else NULL,
  if (has_germline) file.path(work, "Germline_BAF.txt") else NULL,
  chrs = c(1:22))
tmark("loadData", t0)

# GC correction (the standard SNP6 step; without it ASCAT's grid search
# fails to reach the goodness-of-fit threshold on HAPMAP-normalised
# HapMap-referenced logR — verified on the CUE/CHB pilot)
if (!is.null(gc_file) && nzchar(gc_file)) {
  gc_path <- if (file.exists(gc_file)) normalizePath(gc_file) else
    resolve(gc_file)
  if (is.na(gc_path)) stop("GC content file not found: ", gc_file)
  cat(sprintf("[abs] GC correction: %s\n", basename(gc_path)))
  t0 <- Sys.time()
  ascat.data <- ascat.correctLogR(ascat.data, GCcontentfile = gc_path)
  tmark("GC correction", t0)
}
if (use_cache) saveRDS(ascat.data, d1_path)
}

}


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
  calls <- lapply(man$het_calls, function(p) {
    d <- read_call_file(resolve(p))
    d[na_keep, , drop = FALSE]          # same NA-row filter as logR/BAF
  })
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
  # KEY FIX: ascat.loadData reorders probes by genomic position (stock and
  # -fast alike), so a het mask built positionally on the manifest/grid
  # order is SCRAMBLED relative to the loaded matrices — het flags land on
  # hom/uncalled probes, the BAF PCF then segments garbage (b≈0.105
  # everywhere), and every downstream fit degenerates.  Build the mask
  # keyed by (chr, pos) here and map it onto the LOADED object's row order
  # by probe name further down.
  het_flag <- matrix(TRUE, nrow = nrow(markers), ncol = length(calls),
                     dimnames = list(markers$name, man$sample))  # TRUE = hom
  key_m <- paste(markers$chr, markers$pos)
  for (i in seq_along(calls)) {
    hit <- match(paste(calls[[i]][[1]], calls[[i]][[2]]), key_m)
    if (anyNA(hit)) stop("call file has rows outside the marker grid (sample ",
                         man$sample[i], ")")
    cv <- tolower(as.character(calls[[i]][[3]]))
    het_flag[hit, i] <- !(cv %in% c("ab", "1", "het"))  # NA call -> hom
  }
  # ASCAT's contract: germlinegenotypes TRUE = germline homozygous, and the
  # matrix rows must follow the LOADED object's (reordered) row order.
  # Uncalled (NN) markers keep TRUE — excluding them from the het set and
  # treating them as homozygous stretches, as in the manual ASCAT SNP-array
  # workflow (NA here would break the hom-stretch logic and the fit
  # acceptance — verified on the CUE/CHB pilot)
  rn <- rownames(ascat.data$Tumor_BAF)
  gg <- matrix(TRUE, nrow = length(rn), ncol = ncol(het_flag),
               dimnames = list(rn, man$sample))
  for (i in seq_len(ncol(het_flag))) {
    # intersect: ascat.correctLogR keeps only probes present in the GC
    # reference file (official behaviour — our grid loses ~4.3k probes
    # there), so a handful of het probes may have no row in the object
    hn <- intersect(rownames(het_flag)[het_flag[, i] == FALSE], rn)
    gg[hn, i] <- FALSE
  }
  ascat.gg <- list(germlinegenotypes = gg)
}

pcf <- function(gg = NULL) {
  args <- list(ASCATobj = ascat.data, ncores = ncores,
               penalty = penalty, seed = seed)
  if (!is.null(gg)) args$ascat.gg <- gg
  do.call(if (use_fast_ascat) ascat.aspcf_parallel
          else ascat.aspcf_sequential, args)
}
if (!d2_hit) {
t0 <- Sys.time()
if (has_germline) {
  ascat.data <- pcf()
} else if (!is.null(ascat.gg)) {
  ascat.data <- pcf(ascat.gg)
} else {
  gg <- ascat.predictGermlineGenotypes(ascat.data, platform = "AffySNP6",
                                       img.dir = tempdir())
  ascat.data <- pcf(list(germlinegenotypes = gg$germlinegenotypes,
                         failedarrays = gg$failedarrays))
}
tmark("aspcf", t0)
if (use_cache) saveRDS(ascat.data, d2_path)
}

ras <- function(rho = NULL, psi = NULL) {
  args <- list(ASCATobj = ascat.data, ncores = ncores, gamma = gamma,
               max_ploidy = max_ploidy, plot = FALSE, img.dir = tempdir())
  if (!is.null(rho)) { args$rho_manual = rho; args$psi_manual = psi }
  do.call(if (use_fast_ascat) ascat.runAscat_parallel
          else ascat.runAscat_sequential, args)
}
t0 <- Sys.time()
ascat.output <- ras()
tmark("runAscat (automatic)", t0)

# ---------------------------------------------------------------------------
# fallback: samples that fail the automatic grid search (ASCAT's
# percentzero / local-minimum filters reject e.g. near-diploid or
# low-purity samples) are re-fitted at their grid-optimal (rho, psi) —
# the documented ASCAT manual mode
# ---------------------------------------------------------------------------
failed <- unique(as.character(ascat.output$failedarrays))
auto_ok <- setdiff(names(ascat.output$purity), failed)
retry <- setdiff(man$sample, auto_ok)
rho_vec <- rep(NA_real_, nrow(man)); psi_vec <- rep(NA_real_, nrow(man))
names(rho_vec) <- names(psi_vec) <- man$sample
no_plots <- function() {
  noop <- function(...) invisible(NULL)
  for (fn in c("ascat.plotSunrise", "ascat.plotNonRounded",
               "ascat.plotAscatProfile")) {
    try(assignInNamespace(fn, noop, ns = "ASCAT"), silent = TRUE)
  }
}
subset1 <- function(s_name) {
  keep = which(ascat.data$samples == s_name)
  o = ascat.data
  o$Tumor_LogR = ascat.data$Tumor_LogR[, keep, drop = FALSE]
  o$Tumor_BAF = ascat.data$Tumor_BAF[, keep, drop = FALSE]
  o$Tumor_LogR_segmented =
    ascat.data$Tumor_LogR_segmented[, keep, drop = FALSE]
  o$Tumor_BAF_segmented = ascat.data$Tumor_BAF_segmented[keep]
  o$samples = ascat.data$samples[keep]
  o$gender = if (!is.null(ascat.data$gender)) ascat.data$gender[keep] else NULL
  o$failedarrays = character(0)
  o
}
t0 <- Sys.time()
if (length(retry)) {
  for (s_name in retry) {
    s_pos = match(s_name, ascat.data$samples)
    lrrsegm = ascat.data$Tumor_LogR_segmented[, s_name]
    bafsegm = ascat.data$Tumor_BAF_segmented[[s_pos]][, , drop = FALSE]
    if (is.null(bafsegm) || nrow(bafsegm) == 0) next
    rr = lrrsegm[rownames(bafsegm)]
    bb = as.numeric(bafsegm)
    ap = !(ascat.data$SNPpos[rownames(bafsegm), 1] %in% c("X", "Y"))
    if (!any(ap)) next
    s_segs = make_segments(rr[ap], bb[ap])
    dmat = ASCAT:::create_distance_matrix(s_segs, gamma,
                                          min_ploidy = 1.5,
                                          max_ploidy = max_ploidy,
                                          min_purity = 0.1,
                                          max_purity = 1.05)
    # scan the distance matrix with runASCAT's own model, applying its
    # sanity bounds (recalculated ploidy window, rho >= 0.2) but NOT the
    # percentzero / strict-local-minimum filters that reject this cohort
    psi_pos = as.numeric(rownames(dmat))
    rho_pos = as.numeric(colnames(dmat))
    L = as.numeric(s_segs[, "length"])
    sb = as.numeric(s_segs[, "b"]); sr = as.numeric(s_segs[, "r"])
    best_m = Inf; best_rho = NA_real_; best_psi = NA_real_
    for (i in seq_along(psi_pos)) {
      psi_i = psi_pos[i]
      for (j in seq_along(rho_pos)) {
        rho_j = rho_pos[j]
        if (rho_j < 0.2 || rho_j > 1) next
        D = (1 - rho_j) * 2 + rho_j * psi_i
        nA = (rho_j - 1 - (sb - 1) * 2^(sr / gamma) * D) / rho_j
        nB = (rho_j - 1 + sb * 2^(sr / gamma) * D) / rho_j
        ploidy_j = sum((nA + nB) * L) / sum(L)
        if (ploidy_j < 1.5 || ploidy_j > max_ploidy) next
        if (dmat[i, j] < best_m) {
          best_m = dmat[i, j]
          best_rho = rho_j; best_psi = psi_i
        }
      }
    }
    if (is.na(best_rho)) next
    rho_vec[s_name] = best_rho
    psi_vec[s_name] = best_psi
    cat(sprintf("[abs] manual fit %s: rho=%.2f psi=%.2f\n", s_name,
                rho_vec[s_name], psi_vec[s_name]))
  }
  if (any(!is.na(rho_vec[retry]))) {
    # one manual fit per failed sample, parallel over failed samples
    # (scalar rho_manual/psi_manual — runAscat requires length == N)
    fits = parallel::mclapply(retry, function(s_name) {
      if (is.na(rho_vec[s_name])) return(NULL)
      no_plots()
      f = if (use_fast_ascat) ascat.runAscat_fast else ASCAT::ascat.runAscat
      f(subset1(s_name), gamma = gamma, max_ploidy = max_ploidy,
        rho_manual = rho_vec[s_name], psi_manual = psi_vec[s_name])
    }, mc.cores = max(1L, min(ncores, length(retry))))
    names(fits) = retry
    for (s_name in names(fits)) {
      r = fits[[s_name]]
      if (is.null(r) || is.null(r$purity) || is.na(r$purity[[1]])) next
      ascat.output$purity[s_name] = r$purity[[1]]
      ascat.output$aberrantcellfraction[s_name] = r$purity[[1]]
      ascat.output$ploidy[s_name] = r$ploidy[[1]]
      ascat.output$psi[s_name] = r$psi[[1]]
      ascat.output$goodnessOfFit[s_name] = r$goodnessOfFit[[1]]
      ascat.output$segments = rbind(ascat.output$segments, r$segments)
      ascat.output$failedarrays =
        setdiff(ascat.output$failedarrays, s_name)
    }
  }
  tmark("manual-fit fallback", t0)
}
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
  if (s %in% names(rho_vec) && !is.na(rho_vec[s])) {
    note <- paste(note, "manual-fit")
  }
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
