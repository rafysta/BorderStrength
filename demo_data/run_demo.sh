#!/bin/bash
# Smoke test: run BorderStrength on the demo data in both input formats
# and compare against expected_demo_BS.txt. Run from the repository root:
#   sh demo_data/run_demo.sh
set -e
DIR=$(cd "$(dirname "$0")/.." && pwd)
cd "$DIR"
TMP=$(mktemp -d)

echo "[1/2] map / matrix format ..."
sh BS.sh -i demo_data/demo_map.txt.gz -o "$TMP/map_BS.txt" --window 3000

echo "[2/2] hic200-cpp format ..."
sh BS.sh -i demo_data/demo_hic200.txt.gz --bin demo_data/demo_bindef.txt -o "$TMP/hic200_BS.txt" --window 3000

echo "--- diff map vs expected ---"
diff "$TMP/map_BS.txt" demo_data/expected_demo_BS.txt && echo "map: OK (identical)"
echo "--- diff hic200 vs expected ---"
diff "$TMP/hic200_BS.txt" demo_data/expected_demo_BS.txt && echo "hic200: OK (identical)"

rm -rf "$TMP"
echo "Demo finished."
