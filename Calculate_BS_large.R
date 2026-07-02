#!/usr/bin/env Rscript
#==============================================================
# BorderStrength : LARGE domains (cohesin / condensin / TAD)
#
# Faithful, vectorised re-implementation of Draw_borderStrength.R
# (sum-based border strength) for full ICE2-normalised matrices.
#   BS(i) = (A + B) / C
#     A : sum of intra contacts in the upstream   window
#     B : sum of intra contacts in the downstream window
#     C : sum of inter contacts between the two windows
# Then: cap at mean+2SD, mean-centre, call boundaries (local max)
# and group bins into domains (TADs).
#
# Input : dense contact matrix (.rds or .matrix(.gz)) whose
#         row/col names are "chr:start:end".
#==============================================================

suppressWarnings(suppressMessages(library(data.table)))
suppressPackageStartupMessages(library(optparse))
options(scipen = 10)

option_list <- list(
  make_option(c("-i", "--in"),   help = "dense contact matrix (.rds or .matrix.gz), dimnames chr:start:end"),
  make_option(c("-o", "--out"),  default = "NA", help = "output score file"),
  make_option(c("--domain"),     default = "NA", help = "output domain-interval file [default: <out>_domains.txt]"),
  make_option(c("--window"),     default = "100kb", help = "window size in bp [default %default]"),
  make_option(c("--chr"),        default = "NA", help = "chromosome (default: all chromosomes in the matrix)"),
  make_option(c("--start"),      default = "-1", help = "region start (bp, requires --chr)"),
  make_option(c("--end"),        default = "-1", help = "region end (bp, requires --chr)"),
  make_option(c("--threshold"),  default = "0", help = "min centred score to call a boundary [default %default]"),
  make_option(c("--min_close"),  default = "2", help = "half-width (bins) of the local-maximum test [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

parse_bp <- function(x) {
  x <- as.character(x)
  x <- sub("Mb", "000000", x, ignore.case = TRUE)
  x <- sub("kb", "000", x, ignore.case = TRUE)
  as.numeric(x)
}

FILE_in   <- as.character(opt[["in"]])
FILE_out  <- as.character(opt$out)
FILE_dom  <- as.character(opt$domain)
WINDOW    <- parse_bp(opt$window)
CHR       <- as.character(opt$chr)
START     <- as.numeric(opt$start)
END       <- as.numeric(opt$end)
THRESHOLD <- as.numeric(opt$threshold)
MIN_CLOSE <- as.integer(opt$min_close)

if (is.na(FILE_in) || FILE_in == "NA") stop("--in is required")
if (FILE_out == "NA")                  stop("--out is required")
if (FILE_dom == "NA") {
  FILE_dom <- if (grepl("\\.txt$", FILE_out)) sub("\\.txt$", "_domains.txt", FILE_out) else paste0(FILE_out, "_domains.txt")
}

.this_dir <- (function() {
  a <- commandArgs(FALSE); f <- grep("^--file=", a, value = TRUE)
  if (length(f)) dirname(sub("^--file=", "", f[1])) else "."
})()
source(file.path(.this_dir, "BS_common.R"))

#--------------------------------------------------------------
# core : vectorised sum-based border strength for one matrix
#--------------------------------------------------------------
fast_bs_sum <- function(m, w) {
  N <- nrow(m)
  Mz <- m; Mz[is.na(Mz)] <- 0
  Dmax <- min(2L * w, N - 1L)
  cd <- vector("list", Dmax)                      # cd[[d]][k+1] = sum_{j=1..k} m[j, j+d]
  for (d in 1:Dmax) {
    len <- N - d
    diag_d <- Mz[cbind(1:len, (1 + d):N)]
    cd[[d]] <- c(0, cumsum(diag_d))
  }
  Sv <- function(d, lo, hi) {                     # vectorised sum of diagonal d over j in [lo,hi]
    len <- N - d
    lo <- pmax(1L, lo); hi <- pmin(len, hi)
    val <- cd[[d]][hi + 1L] - cd[[d]][lo]
    val[hi < lo] <- 0
    val
  }
  idx <- (w + 1L):(N - w)
  A <- numeric(length(idx)); B <- numeric(length(idx)); C <- numeric(length(idx))
  for (d in 1:w)          A <- A + Sv(d, idx - w, idx - d)
  for (d in 1:w)          B <- B + Sv(d, idx, idx + w - d)
  for (d in 2:Dmax)       C <- C + Sv(d, pmax(idx - w, idx + 1L - d), pmin(idx - 1L, idx + w - d))
  ifelse(C == 0, 1, (A + B) / C)
}

#--------------------------------------------------------------
# per-chromosome pipeline : DI -> cap -> centre -> boundary
#--------------------------------------------------------------
compute_large_chr <- function(m, loc, w, threshold, min_close) {
  N <- nrow(m)
  if (N < 2 * w + 1) return(NULL)                 # not enough bins for the window
  DI <- fast_bs_sum(m, w)                          # length N-2w

  bc_detect <- c(rep(-99, w), DI, rep(-99, w))     # raw BS, sentinel pads for local-max test
  bc_out    <- c(rep(NA_real_, w), DI, rep(NA_real_, w))

  maxV <- mean(DI) + 2 * sd(DI)
  DIc  <- pmin(DI, maxV)
  avg  <- mean(DIc)
  DIfull <- c(rep(avg, w), DIc, rep(avg, w))
  DInorm <- as.numeric(scale(DIfull, scale = FALSE))

  # boundary : centred score > threshold AND raw BS is a local max over +/-min_close
  mc  <- min_close
  sh  <- lapply(-mc:mc, function(s) shift(bc_detect, n = s))
  rmx <- do.call(pmax, c(sh, list(na.rm = TRUE)))
  boundary <- integer(N)
  valid <- (w + 1L):(N - w)
  boundary[valid] <- as.integer(DInorm[valid] > threshold & bc_detect[valid] >= rmx[valid])

  data.table(chr = loc$chr, start = loc$start, end = loc$end,
             BS = bc_out, BS.norm = DInorm, boundary = boundary)
}

#--------------------------------------------------------------
# load matrix, loop chromosomes
#--------------------------------------------------------------
map <- read_dense_matrix(FILE_in)
loc <- parse_locs(rownames(map))
resolution <- loc$end[1] - loc$start[1] + 1
w <- as.integer(round(WINDOW / resolution))
if (w < 1) stop("window (", WINDOW, " bp) < resolution (", resolution, " bp)")

chroms <- unique(loc$chr)
if (CHR != "NA") {
  if (!CHR %in% chroms) stop("chromosome ", CHR, " not found in matrix")
  chroms <- CHR
}

score_all <- list(); dom_all <- list()
for (cc in chroms) {
  sel <- which(loc$chr == cc)
  sel <- sel[order(loc$start[sel])]
  res <- compute_large_chr(map[sel, sel, drop = FALSE], loc[sel, ], w, THRESHOLD, MIN_CLOSE)
  if (is.null(res)) { message("skip ", cc, " (fewer than ", 2 * w + 1, " bins)"); next }
  if (CHR != "NA") {
    if (START != -1) res <- res[start >= START]
    if (END   != -1) res <- res[end   <= END]
  }
  dm <- call_domains(res$chr, res$start, res$end, res$BS.norm, res$boundary,
                     min_size = 2, depth_frac = 0.3)
  res[, TADid := dm$perbin$TADid]
  res[, TAD   := dm$perbin$TAD]
  score_all[[cc]] <- res
  dom_all[[cc]]   <- dm$domains
}

score <- rbindlist(score_all)
domains <- rbindlist(dom_all)
score[, BS      := formatC(BS,      digits = 4, format = "g")]
score[, BS.norm := formatC(BS.norm, digits = 4, format = "g")]

fwrite(score, FILE_out, sep = "\t", quote = FALSE)
write_domains(domains, FILE_dom)
cat(sprintf("BorderStrength(large): %d bins, %d domains (%d TADs) over %s; window=%d bp (w=%d bins @ %d bp)\n",
            nrow(score), nrow(domains), sum(domains$is_TAD), paste(chroms, collapse = ","),
            as.integer(WINDOW), w, resolution))
