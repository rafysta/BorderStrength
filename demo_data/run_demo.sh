#!/bin/bash
# Smoke test: run BorderStrength on the demo data (both methods, both input
# formats) and compare against the expected outputs. Run from the repo root:
#   sh demo_data/run_demo.sh
set -e
DIR=$(cd "$(dirname "$0")/.." && pwd)
cd "$DIR"
TMP=$(mktemp -d)

echo "[micro] map / matrix format ..."
sh BS.sh --method micro -i demo_data/demo_map.txt.gz -o "$TMP/map_BS.txt" --window 3000
diff "$TMP/map_BS.txt" demo_data/expected_demo_BS.txt && echo "  micro map score:   OK"

echo "[micro] hic200-cpp format ..."
sh BS.sh --method micro -i demo_data/demo_hic200.txt.gz --bin demo_data/demo_bindef.txt -o "$TMP/h_BS.txt" --window 3000
diff "$TMP/h_BS.txt" demo_data/expected_demo_BS.txt && echo "  micro hic200 score: OK"
diff "$TMP/h_BS_domains.txt" demo_data/expected_demo_domains.txt && echo "  micro domains:     OK"

echo "[large] full matrix ..."
sh BS.sh --method large -i demo_data/demo_large.matrix -o "$TMP/large_BS.txt" --window 100kb
diff "$TMP/large_BS.txt" demo_data/expected_large_BS.txt && echo "  large score:       OK"
diff "$TMP/large_BS_domains.txt" demo_data/expected_large_domains.txt && echo "  large domains:     OK"

echo "[draw] render a region ..."
Rscript --vanilla Draw_BS.R -i "$TMP/large_BS.txt" -o "$TMP/large.png" --chr II --domain "$TMP/large_BS_domains.txt"
[ -s "$TMP/large.png" ] && echo "  draw png:          OK"

rm -rf "$TMP"
echo "Demo finished."
