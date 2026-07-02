#==============================================================
# BorderStrength : shared helpers
#   - script directory resolution (for source())
#   - dense matrix reader (.rds / .matrix(.gz))
#   - domain calling (TADid + TAD flag + interval table)
# Sourced by Calculate_BS_large.R and Calculate_BS_micro.R
#==============================================================

suppressWarnings(suppressMessages(library(data.table)))

# directory of the running Rscript, so we can source siblings and
# so callers can locate companion files
get_script_dir <- function() {
  args <- commandArgs(FALSE)
  f <- grep("^--file=", args, value = TRUE)
  if (length(f)) return(dirname(sub("^--file=", "", f[1])))
  "."
}

# parse "chr:start:end" row/col labels into a data.frame
parse_locs <- function(labels) {
  p <- tstrsplit(labels, ":", fixed = TRUE)
  data.frame(chr = as.character(p[[1]]),
             start = as.integer(p[[2]]),
             end = as.integer(p[[3]]),
             stringsAsFactors = FALSE)
}

# read a dense contact matrix whose dimnames are "chr:start:end"
#   .rds                -> readRDS
#   .matrix / .matrix.gz-> read.table (gz auto-detected by file())
read_dense_matrix <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    map <- readRDS(path)
  } else {
    map <- as.matrix(read.table(path, header = TRUE, check.names = FALSE))
  }
  map
}

# convert BS.norm score + boundary flags into domains (TADs)
#   faithful port of the domain logic in Draw_borderStrength.R:
#   TADid = cumsum(boundary), boundary bins get id-0.5;
#   a domain is a TAD when it has >= min_size bins AND the mean of its
#   lowest depth_frac fraction of scores is < 0 (a real depleted interior)
# returns list(perbin = data.table(... TADid, TAD),
#              domains = data.table(chr,start,end,TADid,n_bins,is_TAD))
call_domains <- function(chr, start, end, score, boundary,
                         min_size = 2, depth_frac = 0.3) {
  n  <- length(score)
  bd <- ifelse(is.na(boundary), 0L, as.integer(boundary))
  TADid <- cumsum(bd)
  TADid[bd == 1] <- TADid[bd == 1] - 0.5

  ids <- unique(TADid)
  flag <- numeric(0)
  for (id in ids) {
    sel  <- which(TADid == id)
    size <- length(sel)
    sc   <- sort(score[sel])                 # ascending; NAs dropped by sort
    num  <- as.integer(length(sc) * depth_frac)
    idx  <- if (num >= 1) seq_len(num) else 1L
    bottom <- if (length(sc) >= 1) mean(sc[idx], na.rm = TRUE) else NA_real_
    flag[as.character(id)] <-
      if (size >= min_size && !is.na(bottom) && bottom < 0) 1L else 0L
  }
  TAD <- as.integer(flag[as.character(TADid)])

  perbin <- data.table(chr = chr, start = start, end = end,
                       TADid = TADid, TAD = TAD)
  domains <- perbin[, .(start = min(start), end = max(end),
                        n_bins = .N, is_TAD = TAD[1]),
                    by = .(chr, TADid)]
  setorder(domains, chr, start)
  list(perbin = perbin, domains = domains)
}

# write a domains interval table
write_domains <- function(domains, path) {
  out <- domains[, .(chr, start, end, domain_id = TADid, n_bins, is_TAD)]
  fwrite(out, path, sep = "\t", quote = FALSE)
}
