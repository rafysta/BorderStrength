#!/bin/bash
# BorderStrength : border strength (BS) / boundary score from Hi-C data
#   --method large  (default) : cohesin / condensin / TAD from full ICE2 matrices (sum-based)
#   --method micro            : fine-scale domains from diagonal-only data (mean-based)

get_usage(){
	cat <<EOF

Usage : $0 [OPTION]

Description
	-h, --help               show help
	-v, --version            show version

	-m, --method [large|micro]
		large (default) : large/mid domains (cohesin, condensin, TAD) from a
		                  full ICE2-normalised matrix (.rds or .matrix.gz)
		micro           : fine-scale (micro) domains from diagonal-only data
		                  (map/matrix, hic200-cpp with --bin, or .rds)

	-i, --in [file]          input contact matrix / data
	-o, --out [file]         output score file (chr start end BS BS.norm boundary TADid TAD)
	--domain [file]          output domain intervals [default: <out>_domains.txt]

	--bin [file]             (micro only) hic200-cpp bin definition file
	-w, --window [bp]        window size in bp; accepts kb/Mb
	                         (default: large=100kb, micro=3kb)
	-c, --chr [chromosome]   restrict to one chromosome (default: all in the matrix)
	--start [bp]             region start (requires --chr)
	--end [bp]               region end   (requires --chr)

	(large only)
	--threshold [0]          min centred score to call a boundary
	--min_close [2]          half-width (bins) of the local-maximum test
EOF
}

get_version(){ echo "${0} version 2.0"; }

SHORT=hvm:i:o:w:c:
LONG=help,version,method:,in:,out:,domain:,bin:,window:,chr:,start:,end:,threshold:,min_close:
PARSED=`getopt --options $SHORT --longoptions $LONG --name "$0" -- "$@"`
if [ $? -ne 0 ]; then exit 2; fi
eval set -- "$PARSED"

METHOD=large
FILE_BIN=NA
FILE_DOM=NA
CHR=NA
START=-1
END=-1
WINDOW=NA
THRESHOLD=0
MIN_CLOSE=2

while true; do
	case "$1" in
		-h|--help)    get_usage; exit 1 ;;
		-v|--version) get_version; exit 1 ;;
		-m|--method)  METHOD="$2"; shift 2 ;;
		-i|--in)      FILE_IN="$2"; shift 2 ;;
		-o|--out)     FILE_OUT="$2"; shift 2 ;;
		--domain)     FILE_DOM="$2"; shift 2 ;;
		--bin)        FILE_BIN="$2"; shift 2 ;;
		-w|--window)  WINDOW="$2"; shift 2 ;;
		-c|--chr)     CHR="$2"; shift 2 ;;
		--start)      START="$2"; shift 2 ;;
		--end)        END="$2"; shift 2 ;;
		--threshold)  THRESHOLD="$2"; shift 2 ;;
		--min_close)  MIN_CLOSE="$2"; shift 2 ;;
		--) shift; break ;;
		*) echo "Programming error"; exit 3 ;;
	esac
done

DIR_LIB=$(dirname $0)
[ ! -n "${FILE_IN}" ]  && echo "Please specify input (--in)"  && exit 1
[ ! -n "${FILE_OUT}" ] && echo "Please specify output (--out)" && exit 1

case "${METHOD}" in
	large)
		[ "${WINDOW}" = "NA" ] && WINDOW=100kb
		Rscript --vanilla "${DIR_LIB}/Calculate_BS_large.R" \
			-i "${FILE_IN}" -o "${FILE_OUT}" --domain="${FILE_DOM}" \
			--window="${WINDOW}" --chr="${CHR}" --start="${START}" --end="${END}" \
			--threshold="${THRESHOLD}" --min_close="${MIN_CLOSE}"
		;;
	micro)
		[ "${WINDOW}" = "NA" ] && WINDOW=3000
		Rscript --vanilla "${DIR_LIB}/Calculate_BS_micro.R" \
			-i "${FILE_IN}" -o "${FILE_OUT}" --bin="${FILE_BIN}" --domain="${FILE_DOM}" \
			--window="${WINDOW}" --chr="${CHR}" --start="${START}" --end="${END}"
		;;
	*)
		echo "Unknown --method '${METHOD}' (use large or micro)"; exit 1 ;;
esac
