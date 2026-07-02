#!/usr/bin/env Rscript
#==============================================================
# BorderStrength : MICRO domains (fine-scale, diagonal-only data)
#
# Fast mean-based border strength (re-implementation of
# Draw_borderStrength_high_resolution.R):
#   BS(i) = mean( meanUpstream , meanDownstream ) / meanInter
# computed by near-diagonal shifted-vector summation.
#
# Input : diagonal-only contact data
#   * text map/matrix : loc1<TAB>loc2<TAB>score  (loc = chr:start:end)
#   * hic200-cpp out  : bin1<TAB>bin2<TAB>score   (+ --bin definition)
#   * .rds matrix     : rownames "chr:start:end"
#==============================================================

suppressWarnings(suppressMessages(library(data.table)))
suppressPackageStartupMessages(library(optparse))
options(scipen = 10)

option_list <- list(
  make_option(c("-i", "--in"),  help = "diagonal-only contact data (map/matrix, hic200-cpp with --bin, or .rds)"),
  make_option(c("-o", "--out"), default = "NA", help = "output score file"),
  make_option(c("--bin"),       default = "NA", help = "bin definition file (bin chr start end) -> hic200-cpp input"),
  make_option(c("--domain"),    default = "NA", help = "output domain-interval file [default: <out>_domains.txt]"),
  make_option(c("--window"),    default = "3000", help = "window size in bp [default %default]"),
  make_option(c("--chr"),       default = "NA", help = "chromosome (default: all chromosomes present)"),
  make_option(c("--start"),     default = "-1", help = "region start (bp, requires --chr)"),
  make_option(c("--end"),       default = "-1", help = "region end (bp, requires --chr)")
)
opt <- parse_args(OptionParser(option_list = option_list))

parse_bp <- function(x) { x <- as.character(x); x <- sub("Mb","000000",x,ignore.case=TRUE); x <- sub("kb","000",x,ignore.case=TRUE); as.numeric(x) }

FILE_in   <- as.character(opt[["in"]])
FILE_out  <- as.character(opt$out)
FILE_bin  <- as.character(opt$bin)
FILE_dom  <- as.character(opt$domain)
windowSize<- parse_bp(opt$window)
CHR       <- as.character(opt$chr)
START     <- as.numeric(opt$start)
END       <- as.numeric(opt$end)

if (is.na(FILE_in) || FILE_in == "NA") stop("--in is required")
if (FILE_out == "NA")                  stop("--out is required")
if (FILE_dom == "NA") {
  FILE_dom <- if (grepl("\\.txt$", FILE_out)) sub("\\.txt$", "_domains.txt", FILE_out) else paste0(FILE_out, "_domains.txt")
}

.this_dir <- (function() { a <- commandArgs(FALSE); f <- grep("^--file=", a, value = TRUE); if (length(f)) dirname(sub("^--file=","",f[1])) else "." })()
source(file.path(.this_dir, "BS_common.R"))

#--------------------------------------------------------------
# 1. Load contacts -> data.table(chr1, s1, chr2, s2, score) + binSize
#--------------------------------------------------------------
reduce_fun <- "max"
if (FILE_bin != "NA") {
  reduce_fun <- "sum"
  bindef <- fread(FILE_bin); setnames(bindef, 1:4, c("bin", "chr", "start", "end"))
  binSize <- bindef$end[1] - bindef$start[1] + 1
  D <- fread(FILE_in); setnames(D, 1:3, c("bin1", "bin2", "score"))
  D <- merge(D, bindef[, .(bin, chr1 = chr, s1 = start)], by.x = "bin1", by.y = "bin")
  D <- merge(D, bindef[, .(bin, chr2 = chr, s2 = start)], by.x = "bin2", by.y = "bin")
  D <- D[, .(chr1, s1, chr2, s2, score)]
} else if (grepl("\\.rds$", FILE_in, ignore.case = TRUE)) {
  map <- readRDS(FILE_in)
  lab <- tstrsplit(rownames(map), ":", fixed = TRUE)
  rchr <- as.character(lab[[1]]); rstart <- as.integer(lab[[2]]); rend <- as.integer(lab[[3]])
  binSize <- rend[1] - rstart[1] + 1
  idx <- which(!is.na(map) & map != 0, arr.ind = TRUE)
  D <- data.table(chr1 = rchr[idx[, 1]], s1 = rstart[idx[, 1]],
                  chr2 = rchr[idx[, 2]], s2 = rstart[idx[, 2]], score = map[idx])
} else {
  D <- fread(FILE_in); setnames(D, 1:3, c("loc1", "loc2", "score"))
  l1 <- tstrsplit(D$loc1, ":", fixed = TRUE); l2 <- tstrsplit(D$loc2, ":", fixed = TRUE)
  binSize <- as.integer(l1[[3]][1]) - as.integer(l1[[2]][1]) + 1
  D <- data.table(chr1 = as.character(l1[[1]]), s1 = as.integer(l1[[2]]),
                  chr2 = as.character(l2[[1]]), s2 = as.integer(l2[[2]]), score = as.numeric(D$score))
}

D <- D[chr1 == chr2]
w <- as.integer(round(windowSize / binSize))
if (w < 1) stop("window (", windowSize, " bp) < bin size (", binSize, " bp)")

chroms <- unique(D$chr1)
if (CHR != "NA") { if (!CHR %in% chroms) stop("chromosome ", CHR, " not found"); chroms <- CHR }

#--------------------------------------------------------------
# core : mean-based BS for one chromosome (shifted-vector method)
#--------------------------------------------------------------
compute_micro_chr <- function(Dc, binSize, w, chromosome, reduce_fun) {
  origin <- min(Dc$s1, Dc$s2)
  a <- as.integer((pmin(Dc$s1, Dc$s2) - origin) / binSize)
  b <- as.integer((pmax(Dc$s1, Dc$s2) - origin) / binSize)
  DT <- data.table(a = a, b = b, score = Dc$score)[(b - a) <= 2L * w]
  DT <- DT[, .(score = if (reduce_fun == "sum") sum(score) else max(score)), by = .(a, b)]
  DT[, d := b - a]

  Kmax <- max(DT$b); N <- Kmax + w + 1L; Dmax <- 2L * w
  sumv <- replicate(Dmax + 1L, numeric(N), simplify = FALSE)
  cntv <- replicate(Dmax + 1L, numeric(N), simplify = FALSE)
  vd <- DT$d; va <- DT$a + 1L; vs <- DT$score
  for (dd in 0:Dmax) {
    sel <- which(vd == dd)
    if (length(sel)) { sumv[[dd + 1L]][va[sel]] <- vs[sel]; cntv[[dd + 1L]][va[sel]] <- 1 }
  }
  add_shift <- function(dst, src, e) { if (e >= N) return(dst); dst[(e + 1L):N] <- dst[(e + 1L):N] + src[1L:(N - e)]; dst }

  A_sum <- numeric(N); A_cnt <- numeric(N); C_sum <- numeric(N); C_cnt <- numeric(N)
  for (e in 1:w) {
    for (f in 1:e) { d <- e - f; A_sum <- add_shift(A_sum, sumv[[d + 1L]], e); A_cnt <- add_shift(A_cnt, cntv[[d + 1L]], e) }
    for (j in 0:(w - 1)) { d <- e + j; C_sum <- add_shift(C_sum, sumv[[d + 1L]], e); C_cnt <- add_shift(C_cnt, cntv[[d + 1L]], e) }
  }
  B_sum <- numeric(N); B_cnt <- numeric(N)
  if (w < N) { B_sum[1L:(N - w)] <- A_sum[(w + 1L):N]; B_cnt[1L:(N - w)] <- A_cnt[(w + 1L):N] }

  A_mean <- ifelse(A_cnt > 0, A_sum / A_cnt, NA_real_)
  B_mean <- ifelse(B_cnt > 0, B_sum / B_cnt, NA_real_)
  C_mean <- ifelse(C_cnt > 0, C_sum / C_cnt, NA_real_)
  num <- rowMeans(cbind(A_mean, B_mean), na.rm = TRUE); num[is.nan(num)] <- NA_real_
  BS  <- ifelse(is.na(C_mean) | C_mean == 0, NA_real_, num / C_mean)

  t <- 1:Kmax
  res <- data.table(chr = chromosome, start = origin + t * binSize,
                    end = origin + t * binSize + binSize - 1L, BS = BS[t + 1L])
  res[, BS.norm := BS - mean(BS, na.rm = TRUE)]
  x <- res$BS.norm
  nbr <- pmax(0, shift(x, -2), shift(x, -1), shift(x, 1), shift(x, 2), na.rm = TRUE)
  res[, boundary := as.integer(!is.na(x) & x >= nbr)]
  res[is.na(boundary), boundary := 0L]
  res
}

#--------------------------------------------------------------
# loop chromosomes
#--------------------------------------------------------------
score_all <- list(); dom_all <- list()
for (cc in chroms) {
  Dc <- D[chr1 == cc]
  res <- compute_micro_chr(Dc, binSize, w, cc, reduce_fun)
  if (CHR != "NA") {
    if (START != -1) res <- res[start >= START]
    if (END   != -1) res <- res[end   <= END]
  }
  dm <- call_domains(res$chr, res$start, res$end, res$BS.norm, res$boundary, min_size = 2, depth_frac = 0.3)
  res[, TADid := dm$perbin$TADid]; res[, TAD := dm$perbin$TAD]
  score_all[[cc]] <- res; dom_all[[cc]] <- dm$domains
}
score <- rbindlist(score_all); domains <- rbindlist(dom_all)
score[, BS := formatC(BS, digits = 4, format = "g")]
score[, BS.norm := formatC(BS.norm, digits = 4, format = "g")]

fwrite(score, FILE_out, sep = "\t", quote = FALSE)
write_domains(domains, FILE_dom)
cat(sprintf("BorderStrength(micro): %d bins, %d domains (%d TADs) over %s; window=%d bp (w=%d bins @ %d bp)\n",
            nrow(score), nrow(domains), sum(domains$is_TAD), paste(chroms, collapse = ","),
            as.integer(windowSize), w, binSize))
