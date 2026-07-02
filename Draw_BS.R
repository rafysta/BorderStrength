#!/usr/bin/env Rscript
#==============================================================
# BorderStrength : plotting utility
#
# Takes a BorderStrength score file (output of Calculate_BS_*.R)
# and draws the centred score (BS.norm) over a chosen region,
# with boundaries marked. Output png / pdf / eps by extension.
#==============================================================

suppressWarnings(suppressMessages(library(data.table)))
suppressPackageStartupMessages(library(optparse))
options(scipen = 10)

option_list <- list(
  make_option(c("-i", "--in"),  help = "BorderStrength score file (chr start end BS BS.norm boundary ...)"),
  make_option(c("-o", "--out"), help = "output image (.png / .pdf / .eps)"),
  make_option(c("--chr"),       default = "NA", help = "chromosome to draw"),
  make_option(c("--start"),     default = "-1", help = "region start (bp)"),
  make_option(c("--end"),       default = "-1", help = "region end (bp)"),
  make_option(c("--domain"),    default = "NA", help = "optional domain-interval file to shade TADs"),
  make_option(c("--column"),    default = "BS.norm", help = "score column to plot [default %default]"),
  make_option(c("--width"),     default = "1200", help = "image width (px) [default %default]"),
  make_option(c("--height"),    default = "160",  help = "image height (px) [default %default]"),
  make_option(c("--ymin"),      default = "NA", help = "y axis minimum"),
  make_option(c("--ymax"),      default = "NA", help = "y axis maximum")
)
opt <- parse_args(OptionParser(option_list = option_list))

parse_bp <- function(x) { x <- as.character(x); x <- sub("Mb","000000",x,ignore.case=TRUE); x <- sub("kb","000",x,ignore.case=TRUE); as.numeric(x) }

FILE_in  <- as.character(opt[["in"]])
FILE_out <- as.character(opt$out)
CHR      <- as.character(opt$chr)
START    <- parse_bp(opt$start)
END      <- parse_bp(opt$end)
COL      <- as.character(opt$column)
W        <- as.numeric(opt$width)
H        <- as.numeric(opt$height)

if (is.na(FILE_in) || FILE_in == "NA") stop("--in is required")
if (is.na(FILE_out))                   stop("--out is required")

D <- fread(FILE_in)
if (!COL %in% names(D)) stop("column '", COL, "' not found in ", FILE_in)
if (CHR != "NA") D <- D[chr == CHR]
if (START != -1) D <- D[start >= START]
if (END   != -1) D <- D[end   <= END]
if (nrow(D) == 0) stop("no rows to plot after filtering")

x1 <- D$start; x2 <- D$end; xc <- (x1 + x2) / 2
y  <- suppressWarnings(as.numeric(D[[COL]])); y[is.na(y)] <- 0
xlo <- if (START != -1) START else min(x1)
xhi <- if (END   != -1) END   else max(x2)
ylo <- if (opt$ymin != "NA") as.numeric(opt$ymin) else min(0, min(y))
yhi <- if (opt$ymax != "NA") as.numeric(opt$ymax) else max(y)

open_dev <- function(f, w, h) {
  if (grepl("\\.png$", f, ignore.case = TRUE)) { png(f, width = w, height = h); return(invisible()) }
  wi <- if (w > 100) w / 72 else w; hi <- if (h > 100) h / 72 else h
  if (grepl("\\.pdf$", f, ignore.case = TRUE)) {
    pdf(f, width = wi, height = hi, useDingbats = FALSE)
  } else if (grepl("\\.eps$", f, ignore.case = TRUE)) {
    postscript(f, horizontal = FALSE, onefile = FALSE, paper = "special", width = wi, height = hi, family = "Helvetica")
  } else {
    stop("output must end in .png, .pdf or .eps")
  }
}

open_dev(FILE_out, W, H)
par(oma = c(0, 0, 0, 0), mar = c(2, 3, 1, 1))
plot(xc, y, type = "n", xaxs = "i", yaxs = "i", xlab = "", ylab = COL,
     ylim = c(ylo, yhi), xlim = c(xlo, xhi))

# shade TAD domains if provided
if (opt$domain != "NA" && file.exists(opt$domain)) {
  dom <- fread(opt$domain)
  if (CHR != "NA") dom <- dom[chr == CHR]
  if ("is_TAD" %in% names(dom)) dom <- dom[is_TAD == 1]
  if (nrow(dom)) rect(dom$start, ylo, dom$end, yhi, col = "#f0f0f0", border = NA)
}

# score polygons: orange positive, blue negative
for (i in seq_along(y)) {
  polygon(c(x1[i], x2[i], x2[i], x1[i]), c(y[i], y[i], 0, 0),
          col = ifelse(y[i] > 0, "orange", "steelblue"), border = NA)
}
abline(h = 0, col = "grey60")

# boundaries
if ("boundary" %in% names(D)) abline(v = xc[D$boundary == 1], col = "red")
invisible(dev.off())
cat(sprintf("BorderStrength: drew %d bins (%s:%s-%s) to %s\n", nrow(D), CHR,
            format(xlo, scientific = FALSE), format(xhi, scientific = FALSE), FILE_out))
