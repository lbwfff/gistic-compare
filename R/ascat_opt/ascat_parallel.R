## =============================================================================
##  Sample-level parallelism for ASCAT on SNP6.0 arrays.
##
##  ascat.aspcf() and ascat.runAscat() are both loops over samples; nothing is
##  shared between samples except the SNP positions.  Splitting the ASCAT object
##  into per-sample sub-objects and running them in forked workers gives a
##  near-linear speedup, and is fully deterministic (no RNG involved for XX
##  samples; for XY the seed is passed through explicitly).
##
##  Requires: ascat_fast.R sourced (for the *_fast workers) - optional.
## =============================================================================

## On macOS, forking a process whose data.table OpenMP thread pool has been
## started makes the child abort inside the Objective-C runtime as soon as it
## touches Quartz (png()).  Force data.table back to one thread before forking.
.guard_fork_safety <- function() {
  if (requireNamespace("data.table", quietly = TRUE) && data.table::getDTthreads() > 1L) {
    data.table::setDTthreads(1L)
    message("data.table threads set to 1 to keep fork()-based parallelism safe on macOS")
  }
  invisible(TRUE)
}

## Split an ASCAT object into a one-sample sub-object.
.ascat_subset <- function(bc, s, segmented = FALSE) {
  b <- bc
  b$Tumor_LogR <- bc$Tumor_LogR[, s, drop = FALSE]
  b$Tumor_BAF  <- bc$Tumor_BAF[, s, drop = FALSE]
  if (!is.null(bc$Germline_LogR)) b$Germline_LogR <- bc$Germline_LogR[, s, drop = FALSE]
  if (!is.null(bc$Germline_BAF))  b$Germline_BAF  <- bc$Germline_BAF[, s, drop = FALSE]
  b$samples <- bc$samples[s]
  b$gender  <- bc$gender[s]
  b$failedarrays <- NULL
  if (segmented) {
    b$Tumor_LogR_segmented <- bc$Tumor_LogR_segmented[, s, drop = FALSE]
    b$Tumor_BAF_segmented  <- list(bc$Tumor_BAF_segmented[[s]])
  } else {
    b$Tumor_LogR_segmented <- NULL
    b$Tumor_BAF_segmented  <- NULL
  }
  b
}

#' Allele-specific segmentation of all samples in parallel.
#' @param ASCATobj an ASCAT object (from ascat.loadData / ascat.correctLogR)
#' @param ncores number of forked workers (default: all but one core)
#' @param use_fast use ascat.aspcf_fast (default TRUE); set FALSE to keep upstream code
#' @param ascat.gg optional germline genotypes, as in ascat.aspcf.  Each worker
#'   receives only its own sample's column: ascat.aspcf indexes the mask as
#'   `gg[, sample]` where `sample` counts within the sub-object it is given, so
#'   passing the full multi-sample mask would silently give every worker the
#'   FIRST sample's genotypes.
#' @param ... passed on to ascat.aspcf (e.g. penalty)
ascat.aspcf_parallel <- function(ASCATobj, ncores = max(1L, parallel::detectCores() - 1L),
                                 use_fast = TRUE, ascat.gg = NULL, ...) {
  .guard_fork_safety()
  S <- length(ASCATobj$samples)
  ncores <- max(1L, min(as.integer(ncores), S))
  worker <- function(s) {
    f <- if (use_fast) ascat.aspcf_fast else ASCAT::ascat.aspcf
    gg <- ascat.gg
    if (!is.null(gg)) {
      if (!is.null(gg$germlinegenotypes) && ncol(gg$germlinegenotypes) >= s)
        gg$germlinegenotypes <- gg$germlinegenotypes[, s, drop = FALSE]
      if (!is.null(gg$failedarrays) && length(gg$failedarrays) >= s)
        gg$failedarrays <- gg$failedarrays[s]
    }
    r <- f(.ascat_subset(ASCATobj, s), out.dir = NA, ascat.gg = gg, ...)
    ## Return bare numeric payloads.  The probe-ID attributes are by far the
    ## largest part of the result (a 900k-element named vector costs ~50 MB to
    ## serialise back to the parent), so strip them and ship a compact integer
    ## index instead; the parent re-attaches the identical names.
    bm <- r$Tumor_BAF_segmented[[1]]
    list(lr = unname(r$Tumor_LogR_segmented[, 1]),
         baf = unname(bm),
         bafidx = match(rownames(bm), rownames(ASCATobj$Tumor_BAF)))
  }
  res <- if (ncores > 1L) parallel::mclapply(seq_len(S), worker, mc.cores = ncores)
         else lapply(seq_len(S), worker)
  bad <- vapply(res, function(x) inherits(x, "try-error") || is.null(x), logical(1))
  if (any(bad)) stop("ascat.aspcf_parallel failed for samples: ",
                     paste(ASCATobj$samples[bad], collapse = ", "))

  L <- nrow(ASCATobj$Tumor_LogR)
  lrseg <- matrix(NA_real_, L, S, dimnames = list(rownames(ASCATobj$Tumor_LogR), ASCATobj$samples))
  bafseg <- vector("list", S)
  for (s in seq_len(S)) {
    lrseg[, s] <- res[[s]]$lr
    bm <- matrix(res[[s]]$baf, ncol = 1)
    rownames(bm) <- rownames(ASCATobj$Tumor_BAF)[res[[s]]$bafidx]
    bafseg[[s]] <- bm
  }
  ASCATobj$Tumor_LogR_segmented <- lrseg
  ASCATobj$Tumor_BAF_segmented  <- bafseg
  ASCATobj
}

#' Turn the three ASCAT figure functions into no-ops (in this process only).
#' ~2/3 of ascat.runAscat's runtime is spent drawing 900k-point profiles, which
#' is pure overhead if you only need the segments / purity / ploidy.
.no_ascat_plots <- function() {
  noop <- function(...) invisible(NULL)
  .patch_binding("ascat.plotSunrise", noop)
  .patch_binding("ascat.plotNonRounded", noop)
  .patch_binding("ascat.plotAscatProfile", noop)
  invisible(TRUE)
}

#' Purity / ploidy / allele-specific copy number of all samples in parallel.
#' @inheritParams ascat.aspcf_parallel
#' @param img.dir,img.prefix as in ascat.runAscat
#' @param plot draw the sunrise/ASCAT/raw profiles?  FALSE skips the figure
#'   drawing (the PNG files are still created but left empty).  Default TRUE.
#' @param ... passed on to ascat.runAscat
ascat.runAscat_parallel <- function(ASCATobj, ncores = max(1L, parallel::detectCores() - 1L),
                                    use_fast = TRUE, img.dir = ".", img.prefix = "",
                                    plot = TRUE, ...) {
  .guard_fork_safety()
  S <- length(ASCATobj$samples)
  ncores <- max(1L, min(as.integer(ncores), S))
  worker <- function(s) {
    if (!plot) .no_ascat_plots()
    f <- if (use_fast) ascat.runAscat_fast else ASCAT::ascat.runAscat
    r <- f(.ascat_subset(ASCATobj, s, segmented = TRUE), img.dir = img.dir, img.prefix = img.prefix, ...)
    ## ship nA/nB without their (very large) probe-ID attributes
    r$nA <- if (is.null(r$nA)) NULL else unname(r$nA)
    r$nB <- if (is.null(r$nB)) NULL else unname(r$nB)
    r
  }
  res <- if (ncores > 1L) parallel::mclapply(seq_len(S), worker, mc.cores = ncores)
         else lapply(seq_len(S), worker)

  purity <- vapply(res, function(x) if (is.null(x$purity)) NA_real_ else x$purity, numeric(1))
  psi    <- vapply(res, function(x) if (is.null(x$psi)) NA_real_ else x$psi, numeric(1))
  ploidy <- vapply(res, function(x) if (is.null(x$ploidy)) NA_real_ else x$ploidy, numeric(1))
  gof    <- vapply(res, function(x) if (is.null(x$goodnessOfFit)) NA_real_ else x$goodnessOfFit, numeric(1))
  names(purity) <- names(psi) <- names(ploidy) <- names(gof) <- ASCATobj$samples

  nA <- nB <- NULL
  good <- which(!is.na(purity))
  if (length(good) > 0) {
    nA <- do.call(cbind, lapply(res[good], function(x) x$nA))
    nB <- do.call(cbind, lapply(res[good], function(x) x$nB))
    rownames(nA) <- rownames(nB) <- rownames(ASCATobj$Tumor_LogR)
    colnames(nA) <- colnames(nB) <- ASCATobj$samples[good]
  }
  segs     <- do.call(rbind, lapply(res[good], function(x) x$segments))
  segs_raw <- do.call(rbind, lapply(res[good], function(x) x$segments_raw))
  dm <- lapply(res, function(x) x$distance_matrix)
  names(dm) <- ASCATobj$samples

  list(nA = nA, nB = nB, purity = purity, aberrantcellfraction = purity, ploidy = ploidy,
       psi = psi, goodnessOfFit = gof,
       failedarrays = names(purity)[is.na(purity)],
       nonaberrantarrays = unlist(lapply(res, function(x) x$nonaberrantarrays)),
       segments = segs, segments_raw = segs_raw, distance_matrix = dm)
}

## Convenience: run the whole array pipeline in parallel and return everything.
#' @param ncores workers for segmentation and purity estimation
ascat.pipeline_parallel <- function(bc, ncores = max(1L, parallel::detectCores() - 1L),
                                    out.dir = ".", img.prefix = "", ...) {
  t <- Sys.time()
  bc <- ascat.aspcf_parallel(bc, ncores = ncores, ...)
  cat(sprintf("  aspcf (parallel, %d cores): %.1f s\n", ncores,
              as.numeric(difftime(Sys.time(), t, units = "secs"))))
  t <- Sys.time()
  res <- ascat.runAscat_parallel(bc, ncores = ncores, img.dir = out.dir, img.prefix = img.prefix)
  cat(sprintf("  runAscat (parallel, %d cores): %.1f s\n", ncores,
              as.numeric(difftime(Sys.time(), t, units = "secs"))))
  list(ascat.bc = bc, ascat.output = res)
}
