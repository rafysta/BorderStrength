# BorderStrength

BorderStrength calculates a **border strength (BS) score** (boundary score) from
Hi-C data to define domain boundaries, and groups the resulting bins into domains
(TADs). It covers the whole range of domain scales through two methods that share
the same underlying idea — comparing contact frequency *within* two flanking
windows to contact frequency *between* them — but differ in the input data and in
how the score is aggregated.

| method | target domains | biology (example) | input data | typical resolution | window |
| ------ | -------------- | ----------------- | ---------- | ------------------ | ------ |
| `micro` | micro-domains | S. pombe | diagonal-only matrix (e.g. hic200-cpp) | 200–250 bp | 3 kb |
| `large` | cohesin domains | S. pombe | full ICE2 matrix | 5 kb | ~100 kb |
| `large` | condensin domains | S. pombe | full ICE2 matrix | 10 kb | ~100 kb |
| `large` | TADs | human | full ICE2 matrix (per chr) | 40 kb | ~500 kb |

- **`micro`** (mean-based) is a fast re-implementation of
  `Draw_borderStrength_high_resolution.R`: `BS = mean(meanUp, meanDown) / meanInter`,
  computed by near-diagonal shifted-vector summation. Use it for fine-scale data
  that only stores the diagonal band.
- **`large`** (sum-based) is a faithful, vectorised re-implementation of
  `Draw_borderStrength.R`: `BS = (sumUp + sumDown) / sumInter`, then cap at
  mean+2SD, mean-centre, call boundaries, and group bins into domains. Use it for
  full ICE2-normalised matrices.

Both methods write a per-bin score file and a domain-interval file.

# Software Requirement

+ Git, Bash
+ R (>= 4.0). R packages: `data.table`, `optparse`.

# Install guide

```
git clone https://github.com/rafysta/BorderStrength.git
Rscript install_libraries.R
```

# Input formats

`large` (dense, full matrix):

+ **`.rds`** — a dense R matrix whose row/col names are `chr:start:end`.
+ **`.matrix` / `.matrix.gz`** — the text equivalent (same content as the `.rds`).

`micro` (diagonal-only):

+ **text map** — `loc1  loc2  score`, `loc = chr:start:end` (rfy_hic2 map format).
+ **hic200-cpp** — `bin1  bin2  score` (integer bin ids); pass the bin definition
  file (`bin  chr  start  end`, from `make_bin_def2`) with `--bin`.
+ **`.rds`** — dense matrix as above.

A matrix containing several chromosomes (e.g. `ALL.rds`) is processed one
intra-chromosome at a time and the results are concatenated. Restrict to a single
chromosome with `--chr` (required for `--start/--end` region mode).

# Instructions for use

```
sh BS.sh --method <large|micro> -i <input> -o <score file> [options]
```

Examples:

```
# condensin domains (pombe, 10kb full matrix)
sh BS.sh --method large -i ALL.rds        -o condensin_BS.txt --window 100kb

# TADs (human, 40kb per-chromosome matrix)
sh BS.sh --method large -i chr22.matrix.gz -o tad_BS.txt      --window 500kb

# micro-domains (pombe, diagonal-only hic200-cpp output)
sh BS.sh --method micro -i contacts.txt.gz --bin bin_def.txt -o micro_BS.txt --window 3kb
```

Options (see `sh BS.sh --help`):

```
-m, --method   large (default) | micro
-i, --in       input matrix / data
-o, --out      output score file
    --domain   domain-interval output       [default: <out>_domains.txt]
    --bin      (micro) hic200-cpp bin definition file
-w, --window   window size in bp; accepts kb/Mb  [large=100kb, micro=3kb]
-c, --chr      restrict to one chromosome    [default: all]
    --start/--end   region (bp; requires --chr)
    --threshold     (large) min centred score to call a boundary [0]
    --min_close     (large) half-width (bins) of the local-maximum test [2]
```

# Output

Score file (one row per bin):

| chr | start | end | BS | BS.norm | boundary | TADid | TAD |
| --- | ----- | --- | -- | ------- | -------- | ----- | --- |

`BS` = raw border strength, `BS.norm` = mean-centred, `boundary` = 1 at a called
boundary, `TADid` = domain id (boundary bins get id − 0.5), `TAD` = 1 if the domain
qualifies (≥ 2 bins and a depleted interior).

Domain file (one row per domain):

| chr | start | end | domain_id | n_bins | is_TAD |
| --- | ----- | --- | --------- | ------ | ------ |

Filter `is_TAD == 1` for the called domains.

# Plotting

`Draw_BS.R` renders a score file over a region (png / pdf / eps):

```
Rscript Draw_BS.R -i score.txt -o region.png --chr II --start 800000 --end 1200000 \
        --domain score_domains.txt
```

Positive `BS.norm` is drawn in orange, negative in blue, boundaries as red lines,
and called TADs shaded (when `--domain` is given).

# Demo

`demo_data/` holds small synthetic demos for both methods with expected outputs.
Run everything with `sh demo_data/run_demo.sh`. See `demo_data/About_demo_data.txt`.

# License

GPL-3.0
