#!/bin/bash
# usage: het.sh "prefix"
INPUT=$1
# Check if input file exists
if [[ ! -f "$INPUT.het" ]]; then
    echo "Error: Input file $INPUT.het not found"
    exit 1
fi

# Temporary file for heterozygosity rates
TMP="${INPUT}_het.tmp"

# Format: IID heterozygosity_rate
tail -n +2 "$INPUT.het" | awk '{if( $5 + 0 != 0) print $2, ($5 - $3)/$5}' > "$TMP"

# Check if temp file has data
if [[ ! -s "$TMP" ]]; then
    echo "Error: No valid heterozygosity data found"
    rm -f "$TMP"
    exit 1
fi

# Compute mean ± 3SD
STAT=($(awk '{sum += $2; sumsq += ($2)^2}
    END {
        mean = sum/NR;
       	stddev = sqrt(sumsq/NR - mean^2);
       	printf "%f %f\n", mean - 3*stddev, mean + 3*stddev
    }' "$TMP"))
lo=${STAT[0]}
hi=${STAT[1]}

echo "Heterozygosity thresholds: $lo to $hi"

# Output file: fail list in PLINK format (FID IID) - outputs to current directory
awk -v low="$lo" -v high="$hi" '$2 < low || $2 > high {print 0, $1}' "$TMP" > "fail-het.txt"

# Report result
count=$(wc -l < "fail-het.txt")
echo "$count individuals failed heterozygosity filter"

# Cleanup
rm -f "$TMP"
