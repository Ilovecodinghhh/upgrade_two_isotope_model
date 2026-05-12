#!/bin/bash
# Reassemble the Luo 2024 C4 distribution NetCDF from split parts.
# Source: Zenodo record 10516423 (Luo et al., 2024, Nature Communications 15:1219)
# Original file: C4_distribution_NUS_v2.2.nc (249 MB, md5:1716c1c19b3b4071f44acf9675fa3257)
# Split into 50 MB parts for GitHub upload.

set -e
cd "$(dirname "$0")"

OUT="C4_distribution_NUS_v2.2.nc"
if [ -f "$OUT" ]; then
    echo "$OUT already exists, skipping reassembly."
    exit 0
fi

cat C4_distribution_NUS_v2.2.nc.part_a{a,b,c,d,e} > "$OUT"
echo "Reassembled $OUT ($(du -h "$OUT" | cut -f1))"

# Verify md5 if md5sum is available
if command -v md5sum &>/dev/null; then
    HASH=$(md5sum "$OUT" | cut -d' ' -f1)
    EXPECTED="1716c1c19b3b4071f44acf9675fa3257"
    if [ "$HASH" = "$EXPECTED" ]; then
        echo "MD5 verified: $HASH ✓"
    else
        echo "WARNING: MD5 mismatch! Got $HASH, expected $EXPECTED"
    fi
fi
