## =============================================================================
##  Fast SNP-array input for ASCAT.
##
##  ascat.loadData() uses read.table() on 900k-row text files; on a SNP6.0
##  project that is one of the largest single costs.  fread() is ~10x faster and
##  produces the same object (verified identical with identical()).
## =============================================================================

#' Drop-in, faster replacement for ASCAT::ascat.loadData
#' @inheritParams ASCAT::ascat.loadData
#' @export
ascat.loadData_fast <- function(Tumor_LogR_file, Tumor_BAF_file, Germline_LogR_file = NULL,
                                Germline_BAF_file = NULL, chrs = c(1:22, "X", "Y"), gender = NULL,
                                sexchromosomes = c("X", "Y"), genomeVersion = NULL, isTargetedSeq = FALSE) {
  stopifnot(length(isTargetedSeq) == 1 && isTargetedSeq %in% c(TRUE, FALSE))
  ## data.table's OpenMP thread pool, once created, makes a later fork()
  ## (parallel::mclapply) abort on macOS inside the Objective-C runtime as soon
  ## as the child touches Quartz.  Disable it before the first fread; note that
  ## setting it later does not undo an already-created pool.
  data.table::setDTthreads(1L)
  rd <- function(f) {
    ## nThread = 1 is deliberate: data.table's OpenMP threads make a later
    ## fork() (parallel::mclapply) abort on macOS with
    ##   "+[Swift.__SharedStringStorage initialize] ... when fork() was called".
    ## Single-threaded fread is still ~5x faster than read.table.
    dt <- data.table::fread(f, sep = "\t", header = TRUE, showProgress = FALSE,
                            na.strings = c("NA", "NaN", "", "Inf", "-Inf"),
                            data.table = FALSE, check.names = FALSE, nThread = 1)
    ## upstream ascat.loadData returns data.frames, not matrices - keep that
    m <- dt[, -(1:3), drop = FALSE]
    for (cc in seq_len(ncol(m))) m[[cc]] <- as.numeric(m[[cc]])
    rownames(m) <- as.character(dt[[1]])
    ## column layout in the file is ID, Chr, Position, samples...
    pos <- dt[, 2:3, drop = FALSE]
    rownames(pos) <- as.character(dt[[1]])
    list(mat = m, pos = pos)
  }
  tl <- rd(Tumor_LogR_file); tb <- rd(Tumor_BAF_file)
  Tumor_LogR <- tl$mat; Tumor_BAF <- tb$mat
  SNPpos <- tl$pos
  SNPpos <- SNPpos[SNPpos[[1]] %in% chrs, , drop = FALSE]
  chrs <- intersect(chrs, unique(SNPpos[[1]]))

  Tumor_LogR <- Tumor_LogR[rownames(SNPpos), , drop = FALSE]
  Tumor_BAF  <- Tumor_BAF[rownames(SNPpos), , drop = FALSE]

  Germline_LogR <- Germline_BAF <- NULL
  if (!is.null(Germline_LogR_file)) {
    gl <- rd(Germline_LogR_file); gb <- rd(Germline_BAF_file)
    Germline_LogR <- gl$mat[rownames(SNPpos), , drop = FALSE]
    Germline_BAF  <- gb$mat[rownames(SNPpos), , drop = FALSE]
  }

  ## order probes by chromosome, then by position (stable, as upstream)
  SNPorder <- integer(nrow(SNPpos))
  ch <- list(); last <- 0
  for (i in seq_along(chrs)) {
    sel <- which(SNPpos[[1]] == chrs[i])
    o <- sel[order(SNPpos[[2]][sel], method = "radix")]
    ch[[i]] <- (last + 1):(last + length(o))
    SNPorder[ch[[i]]] <- o
    last <- last + length(o)
  }
  SNPpos <- SNPpos[SNPorder, , drop = FALSE]
  Tumor_LogR <- Tumor_LogR[SNPorder, , drop = FALSE]
  Tumor_BAF  <- Tumor_BAF[SNPorder, , drop = FALSE]
  if (!is.null(Germline_LogR)) {
    Germline_LogR <- Germline_LogR[SNPorder, , drop = FALSE]
    Germline_BAF  <- Germline_BAF[SNPorder, , drop = FALSE]
  }
  SNPpos[[1]] <- as.character(SNPpos[[1]])
  SNPpos[[2]] <- as.numeric(SNPpos[[2]])
  rownames(SNPpos) <- rownames(Tumor_LogR)

  chr <- if (!isTargetedSeq) ASCAT:::split_genome(SNPpos) else ch
  if (is.null(gender)) gender <- rep("XX", ncol(Tumor_LogR))
  X_nonPAR <- NULL
  if (!is.null(genomeVersion)) {
    X_nonPAR <- switch(genomeVersion,
                       hg19 = c(2699521, 154931043),
                       hg38 = c(2781480, 155701382),
                       CHM13 = c(2394411, 153925834),
                       stop("genomeVersion must be 'hg19', 'hg38' or 'CHM13'."))
  }
  list(Tumor_LogR = Tumor_LogR, Tumor_BAF = Tumor_BAF,
       Tumor_LogR_segmented = NULL, Tumor_BAF_segmented = NULL,
       Germline_LogR = Germline_LogR, Germline_BAF = Germline_BAF,
       SNPpos = SNPpos, ch = ch, chr = chr, chrs = chrs,
       samples = colnames(Tumor_LogR), gender = gender,
       sexchromosomes = sexchromosomes, X_nonPAR = X_nonPAR,
       isTargetedSeq = isTargetedSeq, failedarrays = NULL)
}
