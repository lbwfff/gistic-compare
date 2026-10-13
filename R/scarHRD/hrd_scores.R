#!/usr/bin/env Rscript
# scarHRD runner for gistic-compare ABS outputs.
#
# Input : allele-specific segmentation table, 8 tab-separated columns:
#         sample, chrom, start, end, totalCN, A, B, ploidy
#         (chrom may carry the "chr" prefix; one row per segment, multiple
#         samples per file are supported)
# Output: <outdir>/all_HRDresults.tsv — per sample:
#         HRD (=HRD-LOH count), Telomeric AI, LST, HRD-sum
#
# Usage: Rscript R/hrd_scores.R <input.tsv> <outdir> [reference]
#        reference: grch37 (= hg19, default) or grch38

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: hrd_scores.R <input.tsv> <outdir> [reference]")
input <- args[1]
outdir <- args[2]
reference <- if (length(args) >= 3) args[3] else "grch37"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

scar_dir <- dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
for (f in c("shrink.seg.ai.R", "shrink.seg.ai.wrapper.R", "preprocess.hrd.R",
            "calc.hrd.R", "calc.ai_new.R", "calc.lst.R", "scar_score.R"))
  source(file.path(scar_dir, f))
e <- new.env(); load(file.path(scar_dir, "sysdata.rda"), envir = e)
list2env(as.list(e), envir = globalenv())

seg <- read.table(input, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
t0 <- Sys.time()
out <- list(); errs <- character(0)
for (s in unique(seg[, 1])) {
  f <- file.path(outdir, paste0("seg1_", s, ".txt"))
  write.table(seg[seg[, 1] == s, ], f, sep = "\t", row.names = FALSE,
              quote = FALSE)
  r <- tryCatch(scar_score(seg = f, reference = reference,
                           chr.in.names = TRUE, seqz = FALSE,
                           outputdir = outdir),
                error = function(e) {
                  errs <<- c(errs, paste(s, conditionMessage(e))); NULL })
  if (!is.null(r)) out[[s]] <- r
  unlink(f)
}
cat(sprintf("scarHRD: %d ok, %d errors, %.1f s\n", length(out),
            length(errs), as.numeric(difftime(Sys.time(), t0,
                                              units = "secs"))))
if (length(errs)) cat(paste(errs, collapse = "\n"), "\n")
if (length(out)) {
  m <- do.call(rbind, out)
  rownames(m) <- names(out)
  write.table(m, file.path(outdir, "all_HRDresults.tsv"), sep = "\t",
              col.names = NA, row.names = TRUE, quote = FALSE)
}
