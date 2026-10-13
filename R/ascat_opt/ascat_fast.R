## =============================================================================
##  ASCAT hot-spot replacements for SNP6.0 array data
##  (ASCAT 3.2.0 / https://github.com/VanLoo-lab/ascat)
##
##  Every function below implements EXACTLY the same algorithm and the same
##  statistics as the corresponding upstream function -- only the implementation
##  changes.  R interpreter loops are replaced by vectorised / compiled code.
##
##  Usage A (recommended -- no change to your existing script):
##      library(ASCAT); library(Rcpp)
##      sourceCpp("aspcf_cpp.cpp")
##      source("ascat_fast.R")
##      patch_ascat()                 # patches the ASCAT namespace in this session
##      ... your existing ascat.loadData / ascat.aspcf / ascat.runAscat script ...
##
##  Usage B (explicit calls):
##      bc <- ascat.aspcf_fast(bc)
##      res <- ascat.runAscat_fast(bc)
##
##  What is replaced, and why:
##    ASCAT:::aspcfpart                          -> compiled C++ (48% of aspcf time)
##    c() growth inside ascat.aspcf              -> preallocated vectors (38%)
##    ASCAT:::predictGermlineHomozygousStretches -> vectorised (5%)
##    ASCAT::make_segments                       -> vectorised
##    ASCAT:::create_distance_matrix             -> vectorised over the rho grid
## =============================================================================

## ---------------------------------------------------------------------------
## 1. make_segments: upstream walks every heterozygous probe in an R for-loop
##    (~300k iterations on a SNP6.0 array).  Consecutive-run detection is
##    vectorised instead.  Same output, including the "length" run counts.
## ---------------------------------------------------------------------------
make_segments_fast <- function(r, b) {
  m <- cbind(r, b)
  m <- m[stats::complete.cases(m), , drop = FALSE]
  n <- nrow(m)
  if (n == 0L) {
    out <- matrix(numeric(0), ncol = 3L)
    colnames(out) <- c("r", "b", "length")
    return(out)
  }
  same <- c(FALSE, (m[-1L, 1L] == m[-n, 1L]) & (m[-1L, 2L] == m[-n, 2L]))
  rl <- rle(cumsum(!same))
  starts <- cumsum(c(1L, rl$lengths[-length(rl$lengths)]))
  out <- cbind(m[starts, 1L], m[starts, 2L], as.numeric(rl$lengths))
  colnames(out) <- c("r", "b", "length")
  out
}

## ---------------------------------------------------------------------------
## 2. create_distance_matrix: upstream is a doubly-nested R loop over
##    101 ploidies x 96 purities = 9,696 interpreter iterations.
##    Key observation: which allele is the minor one does NOT depend on
##    (rho, psi), so it is decided once and the whole rho grid is then
##    evaluated in one vectorised operation per ploidy (101 iterations).
## ---------------------------------------------------------------------------
create_distance_matrix_fast <- function(segments, gamma, min_ploidy = NULL, max_ploidy = NULL,
                                        min_purity = NULL, max_purity = NULL) {
  s <- segments
  if (is.null(min_ploidy) || is.null(max_ploidy)) {
    psi_pos <- seq(1, 6, 0.05)
  } else {
    psi_pos <- seq(min_ploidy - 0.5, max_ploidy + 0.5, 0.05)
  }
  if (is.null(min_purity) || is.null(max_purity)) {
    rho_pos <- seq(0.1, 1.05, 0.01)
  } else {
    rho_pos <- seq(round(min_purity, 2), round(max_purity, 2), 0.01)
  }
  d <- matrix(nrow = length(psi_pos), ncol = length(rho_pos))
  rownames(d) <- psi_pos
  colnames(d) <- rho_pos

  bb <- s[, "b"]; rr <- s[, "r"]; len <- s[, "length"]
  q <- 2^(rr / gamma)
  cw <- ifelse(bb == 0.5, 0.05, 1)
  ## sum(nA) < sum(nB)  <=>  sum((1-2b)*q) < 0, independent of rho and psi
  useA <- sum((1 - 2 * bb) * q, na.rm = TRUE) < 0
  qb <- (bb - 1) * q        # (b-1)*q, exactly as upstream associates it
  qbb <- bb * q             # b*q
  nE <- length(qb)
  for (i in seq_along(psi_pos)) {
    psi <- psi_pos[i]
    ## keep upstream's exact association: (rho-1 - (b-1)*q*C0)/rho
    C0 <- (1 - rho_pos) * 2 + rho_pos * psi
    nA <- (rep(rho_pos - 1, each = nE) - outer(qb, C0)) / rep(rho_pos, each = nE)
    nB <- (rep(rho_pos - 1, each = nE) + outer(qbb, C0)) / rep(rho_pos, each = nE)
    nm <- if (useA) nA else nB
    ## upstream: sum(abs(nMinor - pmax(round(nMinor),0))^2 * length * ifelse(b==0.5,0.05,1))
    A <- (nm - pmax(round(nm), 0))^2
    d[i, ] <- colSums((A * len) * cw, na.rm = TRUE)
  }
  d
}

## ---------------------------------------------------------------------------
## 3. predictGermlineHomozygousStretches: upstream is a per-probe R loop that
##    also grows its index vector with c() (quadratic on long homozygous runs).
##    Exact behaviour preserved, including:
##      * NA probes are skipped and do NOT break a run;
##      * a trailing TRUE-run is only reported when the very last probe of the
##        chunk is itself homozygous;
##      * the minimum run length is in TRUE probes, not in genomic span.
## ---------------------------------------------------------------------------
predictGermlineHomozygousStretches_fast <- function(chr, hom) {
  HomoStretches <- vector("list", ncol(hom))
  for (sam in seq_len(ncol(hom))) {
    homsam <- hom[, sam]
    perchom <- sum(homsam, na.rm = TRUE) / sum(!is.na(homsam))
    homthres <- ceiling(log(0.001, perchom))
    allhprobes <- NULL
    for (chrke in seq_along(chr)) {
      idx <- chr[[chrke]]
      h <- homsam[idx]
      n <- length(h)
      if (n == 0L) next
      keep <- which(!is.na(h))
      if (length(keep) == 0L) next
      v <- h[keep]
      nv <- length(v)
      rid <- cumsum(c(TRUE, v[-1L] != v[-nv]))
      rl <- rle(rid)
      ends_pos <- cumsum(rl$lengths)
      starts_pos <- ends_pos - rl$lengths + 1L
      vals <- v[starts_pos]
      emit <- vals & (seq_along(vals) < length(vals) | keep[ends_pos] == n)
      if (any(emit)) {
        sel <- which(emit)
        ok <- rl$lengths[sel] >= homthres
        if (any(ok)) {
          sel <- sel[ok]
          allhprobes <- rbind(allhprobes,
                              cbind(chrke, idx[keep[starts_pos[sel]]], idx[keep[ends_pos[sel]]]))
        }
      }
    }
    if (is.null(allhprobes)) allhprobes <- rbind(NULL, c(0, 0, 0))
    HomoStretches[[sam]] <- allhprobes
  }
  HomoStretches
}

## ---------------------------------------------------------------------------
## 4. fastAspcf: identical R wrapper; the O(N^2) core aspcfpart is compiled.
##    (aspcfpart_cpp is defined in aspcf_cpp.cpp; if it is not available the
##    pure-R upstream aspcfpart is used instead.)
## ---------------------------------------------------------------------------
fastAspcf_fast <- function(logR, allB, kmin, gamma, isTargetedSeq) {
  if (is.null(isTargetedSeq)) isTargetedSeq <- FALSE
  core <- if (exists("aspcfpart_cpp", inherits = TRUE)) {
    function(...) aspcfpart_cpp(...)
  } else {
    function(logRpart, allBflip, a, b, d, sd1, sd2, N, kmin, gamma)
      aspcfpart(logRpart, allBflip, a, b, d, sd1, sd2, N, kmin, gamma)$breakpts
  }
  N <- length(logR)
  w <- 1000; d <- 100
  startw <- -d; stopw <- w - d
  nseg <- 0; var2 <- 0; var3 <- 0; breakpts <- 0
  repeat {
    from <- max(c(1, startw)); to <- min(c(stopw, N))
    logRpart <- logR[from:to]
    allBpart <- allB[from:to]
    allBflip <- allBpart
    allBflip[allBpart > 0.5] <- 1 - allBpart[allBpart > 0.5]
    sd1 <- getMad(logRpart); sd2 <- getMad(allBflip); sd3 <- getMad(allBpart)
    sd.valid <- c(!is.na(sd1), !is.na(sd2), sd1 != 0, sd2 != 0)
    if (all(sd.valid)) {
      breakptspart <- core(logRpart, allBflip, startw, stopw, d, sd1, sd2, N, kmin, gamma)
      larger <- breakptspart > breakpts[length(breakpts)]
      breakpts <- c(breakpts, breakptspart[larger])
      var2 <- var2 + sd2^2
      var3 <- var3 + sd3^2
      nseg <- nseg + 1
    }
    if (stopw < N + d) {
      startw <- min(stopw - 2 * d + 1, N - 2 * d)
      stopw <- startw + w
    } else break
  }
  breakpts <- unique(c(breakpts, N))
  if (nseg == 0) nseg <- 1
  sd2 <- sqrt(var2 / nseg)
  sd3 <- sqrt(var3 / nseg)
  frst <- breakpts[1:length(breakpts) - 1] + 1
  last <- breakpts[2:length(breakpts)]
  nseg <- length(frst)
  yhat1 <- rep(NA, N); yhat2 <- rep(NA, N)
  for (i in 1:nseg) {
    yhat1[frst[i]:last[i]] <- rep(mean(logR[frst[i]:last[i]]), last[i] - frst[i] + 1)
    yi2 <- allB[frst[i]:last[i]]
    mu <- if (length(yi2) == 0) 0 else mean(abs(yi2 - 0.5))
    if (sqrt(sd2^2 + mu^2) < 2 * sd2) mu <- 0
    if (isTargetedSeq && mu <= 0.05) mu <- 0
    yhat2[frst[i]:last[i]] <- rep(mu + 0.5, last[i] - frst[i] + 1)
  }
  list(yhat1 = yhat1, yhat2 = yhat2)
}

## ---------------------------------------------------------------------------
## 5. ascat.aspcf_fast: identical algorithm; the three hot loops are:
##      (a) the per-probe c() accumulation building logRc   -> preallocated
##      (b) the per-segment c() accumulation building logRd -> preallocated
##      (c) the O(totall^2) breakpoint refinement inner loop -> preallocated d
##    All mean()/which.min() calls are kept exactly as upstream so the output is
##    numerically identical, not merely similar.
## ---------------------------------------------------------------------------
ascat.aspcf_fast <- function(ASCATobj, selectsamples = 1:length(ASCATobj$samples), ascat.gg = NULL,
                             penalty = 70, out.dir = ".", out.prefix = "", seed = as.integer(Sys.time())) {
  set.seed(seed)
  gg <- if (!is.null(ascat.gg)) ascat.gg$germlinegenotypes else
    ASCATobj$Germline_BAF < 0.3 | ASCATobj$Germline_BAF > 0.7
  ghs <- predictGermlineHomozygousStretches_fast(ASCATobj$chr, gg)

  segmentlengths <- unique(c(penalty, 35, 50, 70, 100, 140))
  segmentlengths <- segmentlengths[segmentlengths >= penalty]

  Tumor_LogR_segmented <- matrix(nrow = dim(ASCATobj$Tumor_LogR)[1], ncol = dim(ASCATobj$Tumor_LogR)[2])
  rownames(Tumor_LogR_segmented) <- rownames(ASCATobj$Tumor_LogR)
  colnames(Tumor_LogR_segmented) <- colnames(ASCATobj$Tumor_LogR)
  Tumor_BAF_segmented <- list()

  for (sample in selectsamples) {
    print.noquote(paste("Sample ", ASCATobj$samples[sample], " (", sample, "/", length(ASCATobj$samples), ")", sep = ""))
    logrfilename <- file.path(out.dir, paste(out.prefix, ASCATobj$samples[sample], ".LogR.PCFed.txt", sep = ""))
    baffilename <- file.path(out.dir, paste(out.prefix, ASCATobj$samples[sample], ".BAF.PCFed.txt", sep = ""))
    logRPCFed <- numeric(0); bafPCFed <- numeric(0)

    if (!is.null(ASCATobj$X_nonPAR) && ASCATobj$gender[sample] == "XY") {
      nonPAR_index <- which(ASCATobj$SNPpos[, 1] %in% c("chrX", "X") & ASCATobj$SNPpos[, 2] >= ASCATobj$X_nonPAR[1] &
                              ASCATobj$SNPpos[, 2] <= ASCATobj$X_nonPAR[2] & !is.na(gg[, sample]))
      autosomes_info <- table(gg[which(ASCATobj$SNPpos[, 1] %in% setdiff(ASCATobj$chrs, ASCATobj$sexchromosomes)), sample])
      if (length(nonPAR_index) > 5) {
        gg[nonPAR_index, sample] <- TRUE
        if (!is.null(ASCATobj$Germline_BAF)) {
          DIST <- 1 - sapply(ASCATobj$Germline_BAF[nonPAR_index, sample], function(x) if (x > 0.5) return(x) else return(1 - x))
          gg[nonPAR_index[which(rank(DIST, ties.method = "random") <= round(length(DIST) * (autosomes_info["FALSE"] / sum(autosomes_info))))], sample] <- FALSE
          rm(DIST)
        } else {
          gg[sample(nonPAR_index, round(length(nonPAR_index) * (autosomes_info["FALSE"] / sum(autosomes_info)))), sample] <- FALSE
        }
      }
      rm(nonPAR_index, autosomes_info)
    }

    for (segmentlength in segmentlengths) {
      logRPCFed <- numeric(0); bafPCFed <- numeric(0)
      tbsam <- ASCATobj$Tumor_BAF[, sample]
      names(tbsam) <- rownames(ASCATobj$Tumor_BAF)
      homosam <- gg[, sample]
      for (chrke in 1:length(ASCATobj$chr)) {
        lr <- ASCATobj$Tumor_LogR[ASCATobj$chr[[chrke]], sample]
        lrwins <- vector(mode = "numeric", length = length(lr))
        lrwins[is.na(lr)] <- NA
        lrwins[!is.na(lr)] <- madWins(lr[!is.na(lr)], 2.5, 25)$ywin
        baf <- tbsam[ASCATobj$chr[[chrke]]]
        homo <- homosam[ASCATobj$chr[[chrke]]]
        Select_het <- !homo & !is.na(homo) & !is.na(baf) & !is.na(lr)
        bafsel <- baf[Select_het]
        bafselwinsmirrored <- madWins(ifelse(bafsel > 0.5, bafsel, 1 - bafsel), 2.5, 25)$ywin
        bafselwins <- ifelse(bafsel > 0.5, bafselwinsmirrored, 1 - bafselwinsmirrored)
        indices <- which(Select_het)
        logRaveraged <- NULL
        if (length(indices) != 0) {
          averageIndices <- c(1, (indices[1:(length(indices) - 1)] + indices[2:length(indices)]) / 2, length(lr) + 0.01)
          startindices <- ceiling(averageIndices[1:(length(averageIndices) - 1)])
          endindices <- floor(averageIndices[2:length(averageIndices)] - 0.01)
          if (length(indices) == 1) { startindices <- 1; endindices <- length(lr) }
          logRaveraged <- vector(mode = "numeric", length = length(indices))
          for (i in 1:length(indices)) {
            if (is.na(endindices[i])) endindices[i] <- startindices[i]
            logRaveraged[i] <- mean(lrwins[startindices[i]:endindices[i]], na.rm = TRUE)
          }
        }
        if (length(logRaveraged) > 0) {
          if (length(logRaveraged) < 6) {
            logRASPCF <- rep(mean(logRaveraged), length(logRaveraged))
            bafASPCF <- rep(mean(bafselwinsmirrored), length(logRaveraged))
            if ("isTargetedSeq" %in% names(ASCATobj) && ASCATobj$isTargetedSeq && bafASPCF[1] <= 0.55)
              bafASPCF <- rep(0.5, length(logRaveraged))
          } else {
            PCFed <- fastAspcf_fast(logRaveraged, bafselwins, 6, segmentlength, ASCATobj$isTargetedSeq)
            logRASPCF <- PCFed$yhat1
            bafASPCF <- PCFed$yhat2
          }
          names(bafASPCF) <- names(indices)
          ## ---- map averaged ASPCF levels back onto every probe ---------------
          L <- length(lr)
          np <- length(logRASPCF)
          logRc <- rep(NA_real_, L)
          for (probe in 1:np) {
            i0 <- indices[probe]
            if (probe == 1) logRc[1:i0] <- logRASPCF[probe]
            if (probe == np) {
              if (i0 < L) logRc[(i0 + 1):L] <- logRASPCF[probe]
            } else {
              i1 <- indices[probe + 1]
              totall <- i1 - i0
              if (logRASPCF[probe] == logRASPCF[probe + 1]) {
                logRc[(i0 + 1):i1] <- logRASPCF[probe]
              } else {
                a <- logRASPCF[probe]; b2 <- logRASPCF[probe + 1]
                d <- numeric(totall)          # preallocated (was c() in a loop)
                for (bp in 0:(totall - 1)) {
                  dis <- sum(abs(lr[(1:bp) + i0] - a), na.rm = TRUE)
                  dis <- sum(dis, sum(abs(lr[((bp + 1):totall) + i0] - b2), na.rm = TRUE), na.rm = TRUE)
                  d[bp + 1] <- dis
                }
                bp <- which.min(d) - 1
                if (bp > 0) logRc[(i0 + 1):(i0 + bp)] <- a
                if (bp < totall) logRc[(i0 + bp + 1):i1] <- b2
              }
            }
          }
          if (anyNA(logRc)) stop("internal error: unmapped probes in ASPCF level mapping")
          ## ---- 2nd step: adapt levels (preallocated; same mean()) ----------
          seg <- rle(logRc)$lengths
          ends <- cumsum(seg); starts <- ends - seg + 1L
          logRd <- numeric(L)
          for (i in seq_along(seg)) {
            logRd[starts[i]:ends[i]] <- mean(lr[starts[i]:ends[i]], na.rm = TRUE)
          }
          logRPCFed <- c(logRPCFed, logRd)
          bafPCFed <- c(bafPCFed, bafASPCF)
        } else {
          logRPCFed <- c(logRPCFed, rep(mean(lr, na.rm = TRUE), length(lr)))
        }
        ## ---- correct wrong segments in germline homozygous stretches ------
        homsegs <- ghs[[sample]][ghs[[sample]][, 1] == chrke, ]
        startchr <- min(ASCATobj$chr[[chrke]]); endchr <- max(ASCATobj$chr[[chrke]])
        if (length(homsegs) == 3) homsegs <- t(as.matrix(homsegs))
        if (!is.null(homsegs) && dim(homsegs)[1] != 0) {
          for (i in 1:dim(homsegs)[1]) {
            if (anyNA(homsegs[i, ])) next
            startpos2 <- max(homsegs[i, 2] - 100, startchr); endpos2 <- min(homsegs[i, 3] + 100, endchr)
            startpos3 <- max(homsegs[i, 2] - 5, startchr);   endpos3 <- min(homsegs[i, 3] + 5, endchr)
            towins <- ASCATobj$Tumor_LogR[startpos2:endpos2, sample]
            winsed <- madWins(towins[!is.na(towins)], 2.5, 25)$ywin
            pcfed <- vector(mode = "numeric", length = length(towins))
            pcfed[!is.na(towins)] <- exactPcf(winsed, 6, floor(segmentlength / 4))
            pcfed2 <- pcfed[(startpos3 - startpos2 + 1):(endpos3 - startpos2 + 1)]
            dif <- abs(pcfed2 - logRPCFed[startpos3:endpos3])
            if (!anyNA(dif) && sum(dif > 0.3) > 5) {
              logRPCFed[startpos3:endpos3] <- ifelse(dif > 0.3, pcfed2, logRPCFed[startpos3:endpos3])
            }
          }
        }
      }
      logRPCFed <- fillNA(logRPCFed, zeroIsNA = TRUE)
      ## ---- adapt levels again (preallocated; same mean()/prevlevel logic) --
      seg <- rle(logRPCFed)$lengths
      ends <- cumsum(seg); starts <- ends - seg + 1L
      out <- numeric(length(logRPCFed)); prevlevel <- 0
      for (i in seq_along(seg)) {
        level <- mean(ASCATobj$Tumor_LogR[starts[i]:ends[i], sample], na.rm = TRUE)
        if (is.nan(level)) level <- prevlevel else prevlevel <- level
        out[starts[i]:ends[i]] <- level
      }
      logRPCFed <- out
      names(logRPCFed) <- rownames(ASCATobj$Tumor_LogR)

      if (length(unique(logRPCFed)) < 800) break
    }

    if (!is.na(out.dir)) write.table(logRPCFed, logrfilename, sep = "\t", col.names = FALSE)
    if (!is.na(out.dir)) write.table(bafPCFed, baffilename, sep = "\t", col.names = FALSE)
    bafPCFed <- as.matrix(bafPCFed)
    Tumor_LogR_segmented[, sample] <- logRPCFed
    Tumor_BAF_segmented[[sample]] <- 1 - bafPCFed
  }

  ASCATobj$Tumor_LogR_segmented <- Tumor_LogR_segmented
  ASCATobj$Tumor_BAF_segmented <- Tumor_BAF_segmented
  ASCATobj$failedarrays <- ascat.gg$failedarrays
  ASCATobj
}

## ---------------------------------------------------------------------------
## 6. Wire the fast helpers into an environment whose parent is the ASCAT
##    namespace, so internal helpers (getMad, madWins, exactPcf, fillNA, ...)
##    resolve exactly as upstream, and patch_ascat() can swap them in.
## ---------------------------------------------------------------------------
.ascat_fast_srcenv <- environment()
.ascat_fast_env <- new.env(parent = asNamespace("ASCAT"))
if (exists("aspcfpart_cpp", envir = .ascat_fast_srcenv, inherits = TRUE)) {
  .ascat_fast_env$aspcfpart_cpp <- get("aspcfpart_cpp", envir = .ascat_fast_srcenv, inherits = TRUE)
}
for (.nm in c("make_segments_fast", "create_distance_matrix_fast",
              "predictGermlineHomozygousStretches_fast", "fastAspcf_fast", "ascat.aspcf_fast")) {
  .f <- get(.nm, envir = .ascat_fast_srcenv)
  environment(.f) <- .ascat_fast_env
  assign(.nm, .f, envir = .ascat_fast_env)
  ## also re-bind the user-visible copy, so that calling e.g. ascat.aspcf_fast()
  ## directly resolves ASCAT internals (getMad, madWins, exactPcf, fillNA, ...)
  assign(.nm, .f, envir = .ascat_fast_srcenv)
}
rm(.nm, .f)

## ---------------------------------------------------------------------------
## 7. patch_ascat(): swap the hot upstream functions for the fast ones *inside
##    the installed ASCAT namespace*.  Your existing script then needs no
##    changes at all: ascat.aspcf() and ascat.runAscat() simply get faster.
## ---------------------------------------------------------------------------
## assignInNamespace() only updates the *namespace* binding; the copy that
## library(ASCAT) put on the search path is a separate binding and would keep
## shadowing it.  Patch both.
.patch_binding <- function(name, value) {
  envs <- list(asNamespace("ASCAT"))
  if ("package:ASCAT" %in% search()) envs <- c(envs, list(as.environment("package:ASCAT")))
  for (env in envs) {
    if (!exists(name, envir = env, inherits = FALSE)) next
    locked <- bindingIsLocked(name, env)
    if (locked) unlockBinding(name, env)
    assign(name, value, envir = env)
    if (locked) lockBinding(name, env)
  }
  invisible(TRUE)
}

patch_ascat <- function(verbose = TRUE) {
  e <- .ascat_fast_env
  .patch_binding("make_segments", e$make_segments_fast)
  .patch_binding("create_distance_matrix", e$create_distance_matrix_fast)
  .patch_binding("predictGermlineHomozygousStretches", e$predictGermlineHomozygousStretches_fast)
  .patch_binding("fastAspcf", e$fastAspcf_fast)
  .patch_binding("aspcfpart",
                 function(logRpart, allBflip, a, b, d, sd1, sd2, N, kmin, gamma)
                   list(breakpts = e$aspcfpart_cpp(logRpart, allBflip, a, b, d, sd1, sd2, N, kmin, gamma)))
  ## by far the biggest single win: the c()-growth loops live inside ascat.aspcf
  .patch_binding("ascat.aspcf", e$ascat.aspcf_fast)
  if (verbose) cat("ASCAT patched: ascat.aspcf, make_segments, create_distance_matrix,",
                   "predictGermlineHomozygousStretches, fastAspcf, aspcfpart\n")
  invisible(TRUE)
}

## ---------------------------------------------------------------------------
## 8. Explicit alternatives that do not touch the namespace.
## ---------------------------------------------------------------------------
ascat.runAscat_fast <- function(ASCATobj, ...) {
  ## temporarily patch the two runASCAT helpers, run upstream ascat.runAscat,
  ## then restore.  Guarantees the rest of the code path is byte-for-byte
  ## upstream ASCAT.
  ns <- asNamespace("ASCAT")
  old_ms <- get("make_segments", ns); old_dm <- get("create_distance_matrix", ns)
  on.exit({ assignInNamespace("make_segments", old_ms, ns = "ASCAT")
            assignInNamespace("create_distance_matrix", old_dm, ns = "ASCAT") }, add = TRUE)
  assignInNamespace("make_segments", .ascat_fast_env$make_segments_fast, ns = "ASCAT")
  assignInNamespace("create_distance_matrix", .ascat_fast_env$create_distance_matrix_fast, ns = "ASCAT")
  ASCAT::ascat.runAscat(ASCATobj, ...)
}

runASCAT_fast <- function(lrr, baf, lrrsegmented, bafsegmented, gender, SNPpos, chromosomes, chrnames,
                          sexchromosomes, ...) {
  ns <- asNamespace("ASCAT")
  old_ms <- get("make_segments", ns); old_dm <- get("create_distance_matrix", ns)
  on.exit({ assignInNamespace("make_segments", old_ms, ns = "ASCAT")
            assignInNamespace("create_distance_matrix", old_dm, ns = "ASCAT") }, add = TRUE)
  assignInNamespace("make_segments", .ascat_fast_env$make_segments_fast, ns = "ASCAT")
  assignInNamespace("create_distance_matrix", .ascat_fast_env$create_distance_matrix_fast, ns = "ASCAT")
  ASCAT::runASCAT(lrr, baf, lrrsegmented, bafsegmented, gender, SNPpos, chromosomes, chrnames,
                  sexchromosomes, ...)
}
