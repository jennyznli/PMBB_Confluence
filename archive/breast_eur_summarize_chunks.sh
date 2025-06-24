#!/bin/bash

# Script to summarize SNP filtering across all chunks
# Run this after all extract jobs complete

WDIR="/project/knathans_tecac/jenny/breast"
OUTPUT_DIR="$WDIR/data/imputed/EUR/step2_chunked_pgen"
LOG_DIR="${OUTPUT_DIR}/logs"

echo "====== SNP FILTERING SUMMARY ACROSS ALL CHUNKS ======"
echo "Date: $(date)"
echo ""

# Check if log directory exists
if [[ ! -d "$LOG_DIR" ]]; then
    echo "ERROR: Log directory not found: $LOG_DIR"
    exit 1
fi

# Count total chunks processed
chunk_count=$(ls ${LOG_DIR}/chunk_*_summary.tsv 2>/dev/null | wc -l)
echo "Total chunks processed: $chunk_count"

if [[ $chunk_count -eq 0 ]]; then
    echo "ERROR: No chunk summary files found"
    exit 1
fi

# Create output files
all_chunks_table="${LOG_DIR}/all_chunks_detailed.tsv"
aggregate_summary="${LOG_DIR}/aggregate_filtering_summary.tsv"
temp_file="${LOG_DIR}/temp_combined.tsv"

echo "====== CREATING DETAILED CHUNK TABLE ======"

# Combine all chunk summaries into one large table
echo -e "Chunk\tChr\tChunk_Num\tInitial\tAfter_Extract\tAfter_MAC\tAfter_R2\tAfter_Geno\tSamples" > "$all_chunks_table"

# Sort chunk files by chromosome and chunk number for organized output
for chr in {1..22}; do
    for file in ${LOG_DIR}/chunk_${chr}_*_summary.tsv; do
        if [[ -f "$file" ]]; then
            tail -n +2 "$file" >> "$all_chunks_table"
        fi
    done
done

# Add any remaining files (in case there are other chromosome formats)
for file in ${LOG_DIR}/chunk_*_summary.tsv; do
    if [[ -f "$file" ]]; then
        # Check if this file was already processed above
        chunk_name=$(basename "$file" _summary.tsv | sed 's/chunk_//')
        if ! grep -q "^${chunk_name}" "$all_chunks_table"; then
            tail -n +2 "$file" >> "$all_chunks_table"
        fi
    fi
done

echo "Detailed chunk table created: $all_chunks_table"
echo "Total chunks in table: $(tail -n +2 "$all_chunks_table" | wc -l)"

# Create temporary file for aggregate calculations
cp "$all_chunks_table" "$temp_file"

echo ""
echo "====== CREATING AGGREGATE SUMMARY ======"

# Calculate totals and create aggregate summary
{
    echo "====== AGGREGATE SNP FILTERING SUMMARY ======"
    echo "Date: $(date)"
    echo ""
    
    awk -F'\t' '
    BEGIN {
        print "====== TOTAL COUNTS ACROSS ALL CHUNKS ======"
    }
    NR > 1 {
        initial += $4
        after_extract += $5
        after_mac += $6
        after_r2 += $7
        final += $8
        samples = $9  # Should be same for all chunks
        chunks++
    }
    END {
        printf "Chunks processed: %d\n", chunks
        printf "Samples per chunk: %d\n", samples
        printf "\nVariant counts:\n"
        printf "%-20s: %12s %12s %8s\n", "Step", "Total", "Removed", "% Removed"
        printf "%-20s: %12s %12s %8s\n", "----", "-----", "-------", "---------"
        
        # Initial
        printf "%-20s: %12d %12s %8s\n", "Initial", initial, "-", "-"
        
        # After sample extraction
        removed = initial - after_extract
        pct = (removed * 100.0) / initial
        printf "%-20s: %12d %12d %7.2f%%\n", "After extraction", after_extract, removed, pct
        
        # After MAC
        removed = after_extract - after_mac
        pct = (removed * 100.0) / after_extract
        printf "%-20s: %12d %12d %7.2f%%\n", "After MAC ≥30", after_mac, removed, pct
        
        # After R2
        removed = after_mac - after_r2
        pct = (removed * 100.0) / after_mac
        printf "%-20s: %12d %12d %7.2f%%\n", "After R2 ≥0.2", after_r2, removed, pct
        
        # After geno
        removed = after_r2 - final
        pct = (removed * 100.0) / after_r2
        printf "%-20s: %12d %12d %7.2f%%\n", "After geno <5%%", final, removed, pct
        
        # Overall
        total_removed = initial - final
        total_pct = (total_removed * 100.0) / initial
        printf "\n%-20s: %12d %12d %7.2f%%\n", "OVERALL", final, total_removed, total_pct
    }' "$temp_file"
    
    echo ""
    echo "====== SUMMARY BY CHROMOSOME ======"
    
    awk -F'\t' '
    NR > 1 {
        chr = $2
        chr_initial[chr] += $4
        chr_final[chr] += $8
        chr_chunks[chr]++
    }
    END {
        printf "%-3s %8s %12s %12s %12s %8s\n", "Chr", "Chunks", "Initial", "Final", "Removed", "% Removed"
        printf "%-3s %8s %12s %12s %12s %8s\n", "---", "------", "-------", "-----", "-------", "---------"
        
        for (chr in chr_initial) {
            removed = chr_initial[chr] - chr_final[chr]
            pct = (removed * 100.0) / chr_initial[chr]
            printf "%-3s %8d %12d %12d %12d %7.2f%%\n", chr, chr_chunks[chr], chr_initial[chr], chr_final[chr], removed, pct
        }
    }' "$temp_file" | sort -k1,1n
    
    echo ""
    echo "====== MACHINE-READABLE CHROMOSOME SUMMARY ======"
    echo -e "Chr\tChunks\tInitial\tFinal\tRemoved\tPct_Removed"
    
    awk -F'\t' '
    NR > 1 {
        chr = $2
        chr_initial[chr] += $4
        chr_final[chr] += $8
        chr_chunks[chr]++
    }
    END {
        for (chr in chr_initial) {
            removed = chr_initial[chr] - chr_final[chr]
            pct = (removed * 100.0) / chr_initial[chr]
            print chr "\t" chr_chunks[chr] "\t" chr_initial[chr] "\t" chr_final[chr] "\t" removed "\t" pct
        }
    }' "$temp_file" | sort -k1,1n
    
} > "$aggregate_summary"

# Also display the aggregate summary to stdout
cat "$aggregate_summary"

echo ""
echo "====== FILES CREATED ======"
echo "1. Detailed chunk table: $all_chunks_table"
echo "   - Contains all individual chunk data"
echo "   - $(tail -n +2 "$all_chunks_table" | wc -l) total chunks"
echo ""
echo "2. Aggregate summary: $aggregate_summary" 
echo "   - Contains overall statistics and per-chromosome summaries"
echo ""
echo "3. Individual chunk logs available in: $LOG_DIR"
echo "   - Pattern: chunk_[CHR]_[CHUNK]_snp_counts.log"

# Clean up
rm -f "$temp_file"

echo ""
echo "====== SUMMARY COMPLETE ======"

