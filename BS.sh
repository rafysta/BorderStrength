#!/bin/bash
# BorderStrength : border strength (BS) / boundary score from Hi-C data
# BS = (mean intra-upstream + mean intra-downstream) / mean inter(upstream vs downstream)
# A high BS marks a strong insulation point, i.e. a domain boundary.

get_usage(){
	cat <<EOF

Usage : $0 [OPTION]

Description
	-h, --help
		show help

	-v, --version
		show version

	-i, --in [contact map file]
		input contact map. One of:
		  * text map/matrix : loc1<TAB>loc2<TAB>score  (loc = chr:start:end)
		  * .rds matrix     : rownames "chr:start:end"
		  * hic200-cpp out  : bin1<TAB>bin2<TAB>score   (give --bin as well)

	-o, --out [output file]
		output file (chr start end BS BS.norm boundary)

	--bin [bin definition file]
		bin definition file (bin<TAB>chr<TAB>start<TAB>end) from hic200-cpp
		make_bin_def2. Supplying this switches --in to hic200-cpp bin-id format.

	-w, --window [bp]
		window size in bp (default 3000)

	-c, --chr [chromosome]
		chromosome to analyse (required if input holds >1 chromosome)

	--start [bp]
		region start (optional)

	--end [bp]
		region end (optional)
EOF
}

get_version(){
	echo "${0} version 1.0"
}

SHORT=hvi:o:w:c:
LONG=help,version,in:,out:,bin:,window:,chr:,start:,end:
PARSED=`getopt --options $SHORT --longoptions $LONG --name "$0" -- "$@"`
if [ $? -ne 0 ]; then
	exit 2
fi
eval set -- "$PARSED"

FILE_BIN=NA
CHR=NA
START=-1
END=-1
WINDOW=3000

while true; do
	case "$1" in
		-h|--help)
			get_usage
			exit 1
			;;
		-v|--version)
			get_version
			exit 1
			;;
		-i|--in)
			FILE_IN="$2"
			shift 2
			;;
		-o|--out)
			FILE_OUT="$2"
			shift 2
			;;
		--bin)
			FILE_BIN="$2"
			shift 2
			;;
		-w|--window)
			WINDOW="$2"
			shift 2
			;;
		-c|--chr)
			CHR="$2"
			shift 2
			;;
		--start)
			START="$2"
			shift 2
			;;
		--end)
			END="$2"
			shift 2
			;;
		--)
			shift
			break
			;;
		*)
			echo "Programming error"
			exit 3
			;;
	esac
done

DIR_LIB=$(dirname $0)

[ ! -n "${FILE_IN}" ]  && echo "Please specify input contact map (--in)" && exit 1
[ ! -n "${FILE_OUT}" ] && echo "Please specify output file (--out)"      && exit 1

WINDOW=$(echo "$WINDOW" | sed -e 's/Mb/000kb/' -e 's/kb/000/')

Rscript --vanilla "${DIR_LIB}/Calculate_BS.R" \
	-i "${FILE_IN}" \
	-o "${FILE_OUT}" \
	--bin "${FILE_BIN}" \
	--window "${WINDOW}" \
	--chr "${CHR}" \
	--start "${START}" \
	--end "${END}"
