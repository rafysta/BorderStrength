#!/usr/bin/env Rscript
#==============================================================
# BorderStrength : border strength (BS) / boundary score
#   for fine-scale domain boundary definition from Hi-C data
#
# BS(i) = mean( meanA , meanB ) / meanC
#   meanA : mean contact within the upstream   window [i-w, i)
#   meanB : mean contact within the downstream window [i, i+w)
#   meanC : mean contact between the two windows (insulation)
# High BS => strong insulation => domain boundary.
#
# This is a fast re-implementation of
#   HiC/Draw/Draw_borderStrength_high_resolution.R
# The per-position table scans of the original are replaced by
# near-diagonal shifted-vector summation, giving the same BS
# definition orders of magnitude faster.
#==============================================================

suppressWarnings(suppressMessages(library(data.table)))
suppressPackageStartupMessages(library(optparse))
options(scipen = 10)

option_list <- list(
  make_option(c("-i", "--in"),
              help = "input contact map. text map/matrix (loc1<TAB>loc2<TAB>score, loc=chr:start:end), .rds matrix, or hic200-cpp output (bin1<TAB>bin2<TAB>score) when --bin is given"),
  make_option(c("-o", "--out"), default = "NA", help = "output file"),
  make_option(c("--bin"),       default = "NA",
              help = "bin definition file (bin<TAB>chr<TAB>start<TAB>end) produced by hic200-cpp make_bin_def2. Supplying this switches --in to hic200-cpp bin-id format"),
  make_option(c("--window"),    default = "3000", help = "window size in bp [default %default]"),
  make_option(c("--chr"),       default = "NA",  help = "chromosome to analyse (required if input holds >1 chromosome)"),
  make_option(c("--start"),     default = "-1",  help = "region start (bp)"),
  make_option(c("--end"),       default = "-1",  help = "region end (bp)")
)
opt <- parse_args(OptionParser(option_list = option_list))

windowSize <- as.numeric(opt$window)
FILE_in    <- as.character(opt[["in"]])
FILE_out   <- as.character(opt$out)
FILE_bin   <- as.character(opt$bin)
CHR        <- as.character(opt$chr)
START      <- as.numeric(opt$start)
END        <- as.numeric(opt$end)

if (is.na(FILE_in) || FILE_in == "NA") stop("--in is required")
if (FILE_out == "NA")                  stop("--out is required")

#--------------------------------------------------------------
# 1. Load contacts as data.table(chr1, s1, chr2, s2, score)
#    s1, s2 = bin start coordinates (bp)
#--------------------------------------------------------------
reduce_fun <- "max"          # how to collapse duplicated bin pairs

if (FILE_bin != "NA") {
  ## hic200-cpp format : bin ids + separate bin definition table
  reduce_fun <- "sum"        # hic200 scores are already aggregated counts
  bindef <- fread(FILE_bin)
  setnames(bindef, 1:4, c("bin", "chr", "start", "end"))
  binSize <- bindef$end[1] - bindef$start[1] + 1
  D <- fread(FILE_in)
  setnames(D, 1:3, c("bin1", "bin2", "score"))
  D <- merge(D, bindef[, .(bin, chr1 = chr, s1 = start)], by.x = "bin1", by.y = "bin")
  D <- merge(D, bindef[, .(bin, chr2 = chr, s2 = start)], by.x = "bin2", by.y = "bin")
  D <- D[, .(chr1, s1, chr2, s2, score)]
} else if (grepl("\\.rds$", FILE_in, ignore.case = TRUE)) {
  ## dense matrix saved as .rds, rownames "chr:start:end"
  map <- readRDS(FILE_in)
  rn  <- rownames(map)
  lab <- tstrsplit(rn, ":", fixed = TRUE)
  rchr <- as.character(lab[[1]]); rstart <- as.integer(lab[[2]]); rend <- as.integer(lab[[3]])
  binSize <- rend[1] - rstart[1] + 1
  if (CHR != "NA") {
    keep <- which(rchr == CHR)
    map <- map[keep, keep, drop = FALSE]
    rchr <- rchr[keep]; rstart <- rstart[keep]
  }
  idx <- which(!is.na(map) & map != 0, arr.ind = TRUE)
  D <- data.table(chr1 = rchr[idx[, 1]], s1 = rstart[idx[, 1]],
                  chr2 = rchr[idx[, 2]], s2 = rstart[idx[, 2]],
                  score = map[idx])
} else {
  ## text map / matrix : loc1 loc2 score, loc = chr:start:end
  D <- fread(FILE_in)
  setnames(D, 1:3, c("loc1", "loc2", "score"))
  l1 <- tstrsplit(D$loc1, ":", fixed = TRUE)
  l2 <- tstrsplit(D$loc2, ":", fixed = TRUE)
  D <- data.table(chr1 = as.character(l1[[1]]), s1 = as.integer(l1[[2]]),
                  chr2 = as.character(l2[[1]]), s2 = as.integer(l2[[2]]),
                  score = as.numeric(D$score))
  binSize <- as.integer(l1[[3]][1]) - as.integer(l1[[2]][1]) + 1
}

#--------------------------------------------------------------
# 2. Restrict to a single intra-chromosomal region
#--------------------------------------------------------------
D <- D[chr1 == chr2]
if (CHR != "NA") D <- D[chr1 == CHR]
chroms <- unique(D$chr1)
if (length(chroms) != 1) {
  stop("input spans ", length(chroms),
       " chromosomes; specify one with --chr (intra-chromosomal data required)")
}
chromosome <- chroms[1]
if (START != -1) D <- D[s1 >= START & s2 >= START]
if (END   != -1) D <- D[(s1 + binSize - 1) <= END & (s2 + binSize - 1) <= END]
if (nrow(D) == 0) stop("no contacts left after filtering")

w <- as.integer(round(windowSize / binSize))
if (w < 1) stop("window (", windowSize, " bp) is smaller than the bin size (", binSize, " bp)")

#--------------------------------------------------------------
# 3. Map to integer bin indices, keep near-diagonal band, aggregate
#--------------------------------------------------------------
origin <- min(D$s1, D$s2)
D[, a := pmin(s1, s2)]
D[, b := pmax(s1, s2)]
D[, a := as.integer((a - origin) / binSize)]
D[, b := as.integer((b - origin) / binSize)]
D <- D[(b - a) <= 2L * w]                       # only offsets used by A/B/C
D <- D[, .(score = if (reduce_fun == "sum") sum(score) else max(score)), by = .(a, b)]
D[, d := b - a]

Kmax <- max(D$b)
N    <- Kmax + w + 1L                            # extra w so B(t)=A(t+w) is in range
Dmax <- 2L * w

# offset vectors : sumv[[d+1]][a+1] = aggregated score of pair (a, a+d)
sumv <- vector("list", Dmax + 1L)
cntv <- vector("list", Dmax + 1L)
for (d in 0:Dmax) { sumv[[d + 1L]] <- numeric(N); cntv[[d + 1L]] <- numeric(N) }
vd <- D$d; va <- D$a + 1L; vs <- D$score          # plain vectors (no DT scoping)
for (d in 0:Dmax) {
  sel <- which(vd == d)
  if (length(sel)) {
    sumv[[d + 1L]][va[sel]] <- vs[sel]
    cntv[[d + 1L]][va[sel]] <- 1
  }
}

# helper : dst[t] += src[t-e]   (1-based length-N vectors)
add_shift <- function(dst, src, e) {
  if (e >= N) return(dst)
  dst[(e + 1L):N] <- dst[(e + 1L):N] + src[1L:(N - e)]
  dst
}

#--------------------------------------------------------------
# 4. A (upstream intra), C (cross) by shifted summation; B = A shifted by w
#--------------------------------------------------------------
A_sum <- numeric(N); A_cnt <- numeric(N)
C_sum <- numeric(N); C_cnt <- numeric(N)
for (e in 1:w) {
  for (f in 1:e) {                # upstream block, a=t-e b=t-f (f<=e), offset e-f
    d <- e - f
    A_sum <- add_shift(A_sum, sumv[[d + 1L]], e)
    A_cnt <- add_shift(A_cnt, cntv[[d + 1L]], e)
  }
  for (j in 0:(w - 1)) {          # cross block, a=t-e b=t+j, offset e+j
    d <- e + j
    C_sum <- add_shift(C_sum, sumv[[d + 1L]], e)
    C_cnt <- add_shift(C_cnt, cntv[[d + 1L]], e)
  }
}
# B(t) = A(t+w)
B_sum <- numeric(N); B_cnt <- numeric(N)
if (w < N) {
  B_sum[1L:(N - w)] <- A_sum[(w + 1L):N]
  B_cnt[1L:(N - w)] <- A_cnt[(w + 1L):N]
}

A_mean <- ifelse(A_cnt > 0, A_sum / A_cnt, NA_real_)
B_mean <- ifelse(B_cnt > 0, B_sum / B_cnt, NA_real_)
C_mean <- ifelse(C_cnt > 0, C_sum / C_cnt, NA_real_)

num <- rowMeans(cbind(A_mean, B_mean), na.rm = TRUE)   # mean(c(A,B), na.rm=TRUE)
num[is.nan(num)] <- NA_real_
BS  <- ifelse(is.na(C_mean) | C_mean == 0, NA_real_, num / C_mean)

#--------------------------------------------------------------
# 5. Assemble output for split positions t = 1 .. Kmax
#--------------------------------------------------------------
t   <- 1:Kmax
res <- data.table(
  chr   = chromosome,
  start = origin + t * binSize,
  end   = origin + t * binSize + binSize - 1L,
  BS    = BS[t + 1L]
)
res[, BS.norm := BS - mean(BS, na.rm = TRUE)]

# boundary = local maximum of BS.norm over +/-2 bins and >= 0
x  <- res$BS.norm
s1 <- shift(x, -2); s2 <- shift(x, -1); s3 <- shift(x, 1); s4 <- shift(x, 2)
nbr <- pmax(0, s1, s2, s3, s4, na.rm = TRUE)
res[, boundary := as.integer(!is.na(x) & x >= nbr)]
res[is.na(boundary), boundary := 0L]

res[, BS      := formatC(BS,      digits = 4, format = "g")]
res[, BS.norm := formatC(BS.norm, digits = 4, format = "g")]

fwrite(res, FILE_out, sep = "\t", quote = FALSE)
cat(sprintf("BorderStrength: wrote %d rows for %s (binSize=%d bp, window=%d bp, w=%d bins)\n",
            nrow(res), chromosome, binSize, as.integer(windowSize), w))
