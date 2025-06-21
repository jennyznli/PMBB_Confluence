#!/bin/bash

# Script to analyze GRAF-ANC extraction results
# Counts SNPs in each chunk and provides summary statistics

# Set up paths
OUTPUT_DIR="/home/jennyzli/breast/imputed/chunked_vcf"
LOG_DIR="/home/jennyzli/breast/log"

echo "====== GRAF-ANC EXTRACTION RESULTS ANALYSIS ======"

module load bcftools/1.20 2>/dev/null

output_files=($(ls ${OUTPUT_DIR}/ancestry_snps_*.vcf.gz 2>/dev/null | sort -V))

if [[ ${#output_files[@]} -eq 0 ]]; then
    echo "ERROR: No output VCF files found in ${OUTPUT_DIR}"
    echo "Looking for pattern: ancestry_snps_*.vcf.gz"
    exit 1
fi

echo "Found ${#output_files[@]} output files to analyze"
echo ""

# Initialize counters
total_variants=0
total_chunks_with_data=0
total_chunks_processed=0

# Create summary files
summary_file="${OUTPUT_DIR}/grafanc_extraction_summary.txt"
detailed_file="${OUTPUT_DIR}/grafanc_extraction_position_summary.txt"

echo "# GRAF-ANC Extraction Summary - $(date)" > "$summary_file"
echo "# Chromosome Chunk StartPos StopPos Variants FileSize" > "$detailed_file"

echo "====== DETAILED RESULTS BY CHUNK ======"
printf "%-8s %-6s %-10s %-10s %-8s %-10s\n" "Chr" "Chunk" "StartPos" "StopPos" "SNPs" "FileSize"
echo "------------------------------------------------------------"

for vcf_file in "${output_files[@]}"; do
    # Extract chunk info from filename
    # Format: ancestry_snps_chr1_chunk1.vcf.gz
    basename_file=$(basename "$vcf_file")
    
    if [[ $basename_file =~ ancestry_snps_(chr[^_]+)_chunk([0-9]+)\.vcf\.gz ]]; then
        chrid="${BASH_REMATCH[1]}"
        chunk="${BASH_REMATCH[2]}"
        
        # Count variants in this file
        variant_count=$(bcftools view -H "$vcf_file" 2>/dev/null | wc -l)
        file_size=$(ls -lh "$vcf_file" | awk '{print $5}')
        
        # Try to get position info from log files if available
        start_pos="N/A"
        stop_pos="N/A"
        
        # Look for corresponding log file
        log_pattern="${LOG_DIR}/grafanc_extract.*_${chunk}.out"
        matching_logs=($(ls $log_pattern 2>/dev/null))
        
        if [[ ${#matching_logs[@]} -gt 0 ]]; then
            log_file="${matching_logs[0]}"
            if [[ -f "$log_file" ]]; then
                start_pos=$(grep "start_pos" "$log_file" 2>/dev/null | awk '{print $2}' | head -1)
                stop_pos=$(grep "stop_pos" "$log_file" 2>/dev/null | awk '{print $2}' | head -1)
                [[ -z "$start_pos" ]] && start_pos="N/A"
                [[ -z "$stop_pos" ]] && stop_pos="N/A"
            fi
        fi
        
        # Update counters
        total_variants=$((total_variants + variant_count))
        total_chunks_processed=$((total_chunks_processed + 1))
        
        if [[ $variant_count -gt 0 ]]; then
            total_chunks_with_data=$((total_chunks_with_data + 1))
        fi
        
        # Print results
        printf "%-8s %-6s %-10s %-10s %-8d %-10s\n" "$chrid" "$chunk" "$start_pos" "$stop_pos" "$variant_count" "$file_size"
        
        # Save to detailed file
        echo "$chrid $chunk $start_pos $stop_pos $variant_count $file_size" >> "$detailed_file"
    else
        echo "WARNING: Could not parse filename: $basename_file"
    fi
done

echo ""
echo "====== SUMMARY STATISTICS ======"
printf "%-30s: %d\n" "Total chunks processed" "$total_chunks_processed"
printf "%-30s: %d\n" "Chunks with ancestry SNPs" "$total_chunks_with_data"
printf "%-30s: %d\n" "Chunks with no SNPs" "$((total_chunks_processed - total_chunks_with_data))"
printf "%-30s: %d\n" "Total ancestry SNPs extracted" "$total_variants"

if [[ $total_chunks_with_data -gt 0 ]]; then
    avg_snps=$((total_variants / total_chunks_with_data))
    printf "%-30s: %d\n" "Average SNPs per chunk (with data)" "$avg_snps"
fi

# Calculate total file size
total_size=$(du -sh "${OUTPUT_DIR}/ancestry_snps_"*.vcf.gz 2>/dev/null | tail -1 | awk '{print $1}')
printf "%-30s: %s\n" "Total output size" "${total_size:-N/A}"

# Get sample count from first file
if [[ ${#output_files[@]} -gt 0 ]]; then
    sample_count=$(bcftools query -l "${output_files[0]}" 2>/dev/null | wc -l)
    printf "%-30s: %d\n" "Samples per file" "$sample_count"
fi

echo ""

# Save summary to file
{
    echo "Total chunks processed: $total_chunks_processed"
    echo "Chunks with ancestry SNPs: $total_chunks_with_data"
    echo "Total ancestry SNPs extracted: $total_variants"
    echo "Analysis completed: $(date)"
} >> "$summary_file"

echo "====== CHROMOSOME BREAKDOWN ======"
# Group by chromosome
declare -A chr_counts
declare -A chr_chunks

for vcf_file in "${output_files[@]}"; do
    basename_file=$(basename "$vcf_file")
    if [[ $basename_file =~ ancestry_snps_(chr[^_]+)_chunk([0-9]+)\.vcf\.gz ]]; then
        chrid="${BASH_REMATCH[1]}"
        variant_count=$(bcftools view -H "$vcf_file" 2>/dev/null | wc -l)
        
        chr_counts["$chrid"]=$((${chr_counts["$chrid"]:-0} + variant_count))
        chr_chunks["$chrid"]=$((${chr_chunks["$chrid"]:-0} + 1))
    fi
done

printf "%-12s %-8s %-10s\n" "Chromosome" "Chunks" "Total_SNPs"
echo "--------------------------------"
for chrid in $(printf '%s\n' "${!chr_counts[@]}" | sort -V); do
    printf "%-12s %-8d %-10d\n" "$chrid" "${chr_chunks[$chrid]}" "${chr_counts[$chrid]}"
done

echo ""
echo "====== OUTPUT FILES ======"
printf "%-20s: %s\n" "Summary file" "$summary_file"
printf "%-20s: %s\n" "Detailed results" "$detailed_file"
printf "%-20s: %s\n" "VCF files location" "$OUTPUT_DIR"

echo ""
echo "====== ANALYSIS COMPLETE ======"

# Check if ready for merge
if [[ $total_variants -gt 0 ]]; then
    echo ""
    echo "✓ Ready for merge step!"
    echo "To merge all chunks:"
    echo "  bsub < grafanc_merge.bsub"
else
    echo ""
    echo "⚠ WARNING: No ancestry SNPs found in any chunks"
    echo "Check extraction logs for issues"
fi
