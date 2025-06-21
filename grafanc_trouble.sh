#!/bin/bash

# Quick script to check sample order across all chromosomes
OUTPUT_DIR="/home/jennyzli/breast/imputed/chunked_vcf"

echo "====== CHECKING SAMPLE ORDER ACROSS ALL CHROMOSOMES ======"
echo "Taking first chunk from each chromosome..."
echo ""

# Load bcftools
module load bcftools/1.20 2>/dev/null

# Check chromosomes 1-22 and X
for chrom in {1..22} X; do
    chrid="chr${chrom}"
    
    # Find first chunk for this chromosome
    first_chunk=$(ls ${OUTPUT_DIR}/ancestry_snps_${chrid}_chunk*.vcf.gz 2>/dev/null | sort -V | head -1)
    
    if [[ -f "$first_chunk" ]]; then
        echo "=== $chrid ($(basename $first_chunk)) ==="
        echo "First 10 samples:"
        bcftools query -l "$first_chunk" | head -10
        echo ""
    else
        echo "=== $chrid ==="
        echo "No chunks found for $chrid"
        echo ""
    fi
done

echo "====== COMPARISON SUMMARY ======"

# Save reference (chr1) sample order
ref_file=$(ls ${OUTPUT_DIR}/ancestry_snps_chr1_chunk*.vcf.gz 2>/dev/null | sort -V | head -1)

if [[ -f "$ref_file" ]]; then
    echo "Using chr1 as reference: $(basename $ref_file)"
    bcftools query -l "$ref_file" > /tmp/ref_samples.txt
    ref_count=$(wc -l < /tmp/ref_samples.txt)
    echo "Reference has $ref_count samples"
    echo ""
    
    # Compare each chromosome to reference
    for chrom in {2..22} X; do
        chrid="chr${chrom}"
        test_file=$(ls ${OUTPUT_DIR}/ancestry_snps_${chrid}_chunk*.vcf.gz 2>/dev/null | sort -V | head -1)
        
        if [[ -f "$test_file" ]]; then
            bcftools query -l "$test_file" > /tmp/test_samples.txt
            test_count=$(wc -l < /tmp/test_samples.txt)
            
            printf "%-5s: " "$chrid"
            
            if [[ $test_count -ne $ref_count ]]; then
                printf "DIFFERENT COUNT (%d vs %d)\n" "$test_count" "$ref_count"
            elif diff /tmp/ref_samples.txt /tmp/test_samples.txt > /dev/null; then
                printf "✓ IDENTICAL ORDER\n"
            else
                printf "✗ DIFFERENT ORDER (same count)\n"
            fi
        else
            printf "%-5s: NO FILES FOUND\n" "$chrid"
        fi
    done
    
    # Clean up temp files
    rm -f /tmp/ref_samples.txt /tmp/test_samples.txt
    
else
    echo "ERROR: No chr1 files found for reference"
fi

echo ""
echo "====== CHECK COMPLETE ======"
