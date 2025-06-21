#!/bin/bash

# Script to count total variants extracted from all grafanc_extract log files

LOG_DIR="/home/jennyzli/breast/log"
cd "$LOG_DIR"

echo "====== COUNTING TOTAL EXTRACTED VARIANTS ======"
echo "Searching in directory: $LOG_DIR"
echo "Start time: $(date)"
echo ""

# Find all grafanc_extract .out files
out_files=(grafanc_extract.*.out)
total_files=${#out_files[@]}

echo "Found $total_files log files to process"
echo ""

# Initialize counters
total_variants=0
files_with_variants=0
files_without_variants=0
files_with_errors=0

# Create temporary files for detailed output
variants_summary="variants_summary.txt"
error_summary="error_summary.txt"

echo "Processing files..."
echo "Job_ID,Chunk,Variants_Extracted" > "$variants_summary"

# Process each file
for file in "${out_files[@]}"; do
    if [[ -f "$file" ]]; then
        # Extract variants_extracted line
        variant_line=$(grep "variants_extracted" "$file" 2>/dev/null)
        
        if [[ -n "$variant_line" ]]; then
            # Extract the number (handles various formats)
            variants=$(echo "$variant_line" | grep -o '[0-9]\+' | tail -1)
            
            if [[ -n "$variants" && "$variants" -gt 0 ]]; then
                total_variants=$((total_variants + variants))
                files_with_variants=$((files_with_variants + 1))
                
                # Extract job info for detailed log
                job_id=$(echo "$file" | sed 's/grafanc_extract\.\([0-9_]*\)\.out/\1/')
                echo "$job_id,unknown,$variants" >> "$variants_summary"
                
                # Print progress every 100 files
                if (( files_with_variants % 100 == 0 )); then
                    echo "Processed $files_with_variants files with variants..."
                fi
            else
                files_without_variants=$((files_without_variants + 1))
            fi
        else
            # Check if file has errors
            if grep -q "ERROR\|FAILED\|No ancestry SNPs found" "$file" 2>/dev/null; then
                echo "$file" >> "$error_summary"
                files_with_errors=$((files_with_errors + 1))
            else
                files_without_variants=$((files_without_variants + 1))
            fi
        fi
    fi
done

echo ""
echo "====== SUMMARY RESULTS ======"
printf "%-30s: %d\n" "Total log files processed" "$total_files"
printf "%-30s: %d\n" "Files with variants" "$files_with_variants"
printf "%-30s: %d\n" "Files without variants" "$files_without_variants"
printf "%-30s: %d\n" "Files with errors" "$files_with_errors"
echo ""
printf "%-30s: %'d\n" "TOTAL VARIANTS EXTRACTED" "$total_variants"
echo ""

# Show top 10 chunks with most variants
echo "====== TOP 10 CHUNKS BY VARIANT COUNT ======"
if [[ -f "$variants_summary" ]]; then
    tail -n +2 "$variants_summary" | sort -t',' -k3 -nr | head -10 | \
    while IFS=',' read job chunk variants; do
        printf "Job %-15s: %'8d variants\n" "$job" "$variants"
    done
fi

echo ""
echo "====== FILES CREATED ======"
echo "Detailed summary: $variants_summary"
if [[ -f "$error_summary" ]]; then
    echo "Error files: $error_summary"
fi

echo ""
echo "End time: $(date)"
echo "====== ANALYSIS COMPLETE ======"

# Quick verification command
echo ""
echo "====== VERIFICATION ======"
echo "You can verify this count with:"
echo "grep 'variants_extracted' grafanc_extract.*.out | grep -o '[0-9]\+' | paste -sd+ | bc"
