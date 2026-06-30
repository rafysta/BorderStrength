# BorderStrength

BorderStrength calculates a **border strength (BS) score**, also called a boundary
score, from Hi-C data to define fine-scale domain boundaries. For every bin position
along a chromosome it compares the average contact frequency *within* the flanking
upstream and downstream windows to the average contact frequency *between* them:

```
BS(i) = mean( meanUpstream , meanDownstream ) / meanInter
```

A high BS marks a strong insulation point, i.e. a domain boundary. The script also
reports a mean-centred score (`BS.norm`) and a 0/1 `boundary` flag set where `BS.norm`
is a local maximum over a +/-2 bin neighbourhood.

This is a standalone, faster re-implementation of `Draw_borderStrength_high_resolution.R`
from the HiC package: the per-position table scans are replaced by near-diagonal
shifted-vector summation, giving the same BS definition with a large speed-up.

# Hardware Requirement

BorderStrength is designed for a Linux environment but runs on any OS with the
appropriate setup. A standard computer with a few GB of RAM is sufficient; memory use
scales with the number of bins on the largest chromosome, not with the square of it.

# Software Requirement

The versions in parentheses were used during development:

+ Git (2.34.1) -- https://git-scm.com/
+ Bash (5.1.16) -- https://www.gnu.org/software/bash/
+ R (4.1.2) -- https://www.r-project.org/

R package dependencies: `data.table`, `optparse`.

# Install guide

```
git clone https://github.com/rafysta/BorderStrength.git
```

Install the required R packages:

```
Rscript install_libraries.R
```

# Input formats

BorderStrength accepts three input layouts and auto-detects them:

1. **Text map / matrix** (default) -- tab-separated `loc1  loc2  score`, where each
   `loc` is `chr:start:end`. This is the map format produced by the rfy_hic2 package
   (e.g. `I_500000.txt.gz`).

2. **hic200-cpp output** -- the `bin1  bin2  score` contact table produced by
   [hic200-cpp](https://github.com/mbyamaguchi/hic200-cpp). Bin ids are integers, so
   you must also pass the bin definition file (`bin  chr  start  end`) made by
   `make_bin_def2` via `--bin`. BorderStrength joins them to recover coordinates.

3. **.rds matrix** -- a dense R matrix whose row/column names are `chr:start:end`.

# Instructions for use

```
sh BS.sh -i <contact map> -o <output file> [--bin <bin def>] [--window 3000] [--chr II]
```

Examples:

```
# rfy_hic2 map format, whole chromosome
sh BS.sh -i I_500000.txt.gz -o I_BS.txt --window 3000

# hic200-cpp format
sh BS.sh -i contacts.txt.gz --bin bin_def.txt -o BS.txt --window 3000 --chr II

# restrict to a region
sh BS.sh -i I_500000.txt.gz -o region_BS.txt --window 3000 --start 800000 --end 900000
```

Options (run `sh BS.sh --help` for the full list):

```
-i, --in       input contact map
-o, --out      output file
    --bin      bin definition file (switches --in to hic200-cpp format)
-w, --window   window size in bp (default 3000; accepts e.g. 3kb)
-c, --chr      chromosome to analyse (required if input has >1 chromosome)
    --start    region start (bp, optional)
    --end      region end (bp, optional)
```

The script may also be called directly:

```
Rscript Calculate_BS.R -i <in> -o <out> [--bin <bin def>] --window 3000
```

# Output

A tab-separated file with one row per bin:

| column   | description                                              |
| -------- | -------------------------------------------------------- |
| chr      | chromosome                                               |
| start    | bin start (bp)                                           |
| end      | bin end (bp)                                             |
| BS       | border strength score                                    |
| BS.norm  | mean-centred BS                                          |
| boundary | 1 if BS.norm is a local maximum (+/-2 bins) and >= 0     |

# Demo

A small synthetic demo (chr II, 250 bp, three engineered boundaries) is in
`demo_data/`, provided in both input formats with the expected output. See
`demo_data/About_demo_data.txt`. Each demo command runs in a few seconds.

# License

GPL-3.0
