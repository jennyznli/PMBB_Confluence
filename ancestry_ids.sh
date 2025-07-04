#!/bin/bash

# Script to create ancestry-specific ID files from grafanc_extract.txt
# Creates files for all PMBB samples (not just breast cancer)

HOME_DIR="/project/knathans_tecac/jenny/breast"
GRAF_DIR="/project/knathans_tecac/jenny/breast/analysis/grafanc_results"
INPUT_FILE="$GRAF_DIR/grafanc_extract.txt"
OUTPUT_DIR="$HOME_DIR/analysis/ancestry"
BREAST_CASE_CONTROL_IDS="$HOME_DIR/data/case_control_all_ids.txt"

echo "=== CREATING ANCESTRY-SPECIFIC ID FILES ==="
echo "Input file: $INPUT_FILE"
echo "Output directory: $OUTPUT_DIR"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: $INPUT_FILE not found!"
    echo "Please make sure grafanc_extract.txt exists in $GRAF_DIR"
    exit 1
fi

# Get unique ancestry codes from the file
echo "Finding unique ancestry codes..."
unique_ancestries=$(cut -f3 "$INPUT_FILE" | tail -n +2 | sort | uniq | grep -v "^$")

echo "Found ancestry codes: $unique_ancestries"
echo ""

# Initialize counters
total_samples=0
total_breast_ancestry_samples=0

# Sort the breast cancer file once for efficiency (if it exists)
temp_breast_sorted=""
if [ -f "$BREAST_CASE_CONTROL_IDS" ]; then
    echo "Preparing breast cancer overlap analysis..."
    temp_breast_sorted="$OUTPUT_DIR/temp_breast_sorted.txt"
    sort "$BREAST_CASE_CONTROL_IDS" > "$temp_breast_sorted"
    echo "Breast cancer case-control file: $BREAST_CASE_CONTROL_IDS"
    echo ""
fi

# Process each ancestry group (combined loop for efficiency)
for ancestry in $unique_ancestries; do
    echo "Processing $ancestry ancestry..."
    
    # Create simple ID file (just sample IDs)
    ids_file="$OUTPUT_DIR/${ancestry}_ids.txt"
    awk -F'\t' -v anc="$ancestry" 'NR>1 && $3==anc {print $1}' "$INPUT_FILE" > "$ids_file"
    
    # Create PLINK format file (0 in first column, ID in second)
    pids_file="$OUTPUT_DIR/${ancestry}_pids.txt"
    awk -F'\t' -v anc="$ancestry" 'NR>1 && $3==anc {print "0\t" $1}' "$INPUT_FILE" > "$pids_file"
    
    # Count samples for this ancestry
    count=$(wc -l < "$ids_file")
    total_samples=$((total_samples + count))
    
    echo "  - Created $ids_file ($count samples)"
    echo "  - Created $pids_file ($count samples)"
    
    # Create breast cancer overlap files if breast cancer data exists
    if [ -f "$temp_breast_sorted" ]; then
        breast_ancestry_ids_file="$OUTPUT_DIR/${ancestry}_breast_ids.txt"
        breast_ancestry_pids_file="$OUTPUT_DIR/${ancestry}_breast_pids.txt"
        
        # Find intersection with breast cancer samples (using ID files)
        sort "$ids_file" | comm -12 - "$temp_breast_sorted" > "$breast_ancestry_ids_file"
        
        # Create PID version (add "0" prefix)
        awk '{print "0\t" $1}' "$breast_ancestry_ids_file" > "$breast_ancestry_pids_file"
        
        # Count overlapping samples
        overlap_count=$(wc -l < "$breast_ancestry_ids_file")
        total_breast_ancestry_samples=$((total_breast_ancestry_samples + overlap_count))
        
        echo "  - Created $breast_ancestry_ids_file ($overlap_count breast cancer samples)"
        echo "  - Created $breast_ancestry_pids_file ($overlap_count breast cancer samples)"
    fi
    echo ""
done

# Clean up temp file
if [ -f "$temp_breast_sorted" ]; then
    rm "$temp_breast_sorted"
fi

# Create summary file
summary_file="$OUTPUT_DIR/ancestry_summary.txt"
echo "=== ANCESTRY SAMPLE SUMMARY ===" > "$summary_file"
echo "Generated on: $(date)" >> "$summary_file"
echo "Total samples processed: $total_samples" >> "$summary_file"
if [ -f "$BREAST_CASE_CONTROL_IDS" ]; then
    echo "Total breast cancer samples with ancestry data: $total_breast_ancestry_samples" >> "$summary_file"
fi
echo "" >> "$summary_file"
echo "Sample counts by ancestry:" >> "$summary_file"
for ancestry in $unique_ancestries; do
    count=$(wc -l < "$OUTPUT_DIR/${ancestry}_ids.txt")
    breast_count=0
    breast_file="$OUTPUT_DIR/${ancestry}_breast_ids.txt"
    if [ -f "$breast_file" ]; then
        breast_count=$(wc -l < "$breast_file")
    fi
    
    if [ -f "$BREAST_CASE_CONTROL_IDS" ]; then
        printf "%-5s: %d total samples, %d breast cancer samples\n" "$ancestry" "$count" "$breast_count" >> "$summary_file"
    else
        printf "%-5s: %d samples\n" "$ancestry" "$count" >> "$summary_file"
    fi
done
echo "" >> "$summary_file"

echo "=== SUMMARY ==="
echo "Total samples processed: $total_samples"
if [ -f "$BREAST_CASE_CONTROL_IDS" ]; then
    echo "Total breast cancer samples with ancestry data: $total_breast_ancestry_samples"
fi
echo "Ancestry summary saved to: $summary_file"

# Add breast cancer case-control ancestry analysis
echo "=== BREAST CANCER CASE-CONTROL ANCESTRY ANALYSIS ==="

CASE_CONTROL_FILE="$HOME_DIR/data/case_control.txt"
if [ -f "$CASE_CONTROL_FILE" ]; then
    echo "Performing detailed breast cancer case-control ancestry analysis..."
    
    # Create breast cancer ancestry summary
    breast_summary_file="$OUTPUT_DIR/breast_cancer_ancestry_summary.txt"
    echo "=== BREAST CANCER CASE-CONTROL ANCESTRY DEMOGRAPHICS ===" > "$breast_summary_file"
    echo "Generated on: $(date)" >> "$breast_summary_file"
    echo "Case-control file: $CASE_CONTROL_FILE" >> "$breast_summary_file"
    echo "Ancestry file: $INPUT_FILE" >> "$breast_summary_file"
    echo "" >> "$breast_summary_file"
    
    # Join case-control status with ancestry data
    temp_joined="$OUTPUT_DIR/temp_breast_ancestry.txt"
    awk -F'\t' 'BEGIN {OFS="\t"}
        # Read case-control file
        NR==FNR && NR>1 {
            status[$1] = $2
            next
        }
        # Process ancestry file
        NR>FNR && NR>1 {
            if ($1 in status) {
                print $1, $3, status[$1]
            }
        }
    ' "$CASE_CONTROL_FILE" "$INPUT_FILE" > "$temp_joined"
    
    # Count samples
    total_breast_samples=$(wc -l < "$temp_joined")
    cases_total=$(awk '$3=="Case"' "$temp_joined" | wc -l)
    controls_total=$(awk '$3=="Control"' "$temp_joined" | wc -l)
    
    echo "SAMPLE COUNTS:" >> "$breast_summary_file"
    echo "Total breast cancer samples with ancestry data: $total_breast_samples" >> "$breast_summary_file"
    echo "Cases: $cases_total" >> "$breast_summary_file"
    echo "Controls: $controls_total" >> "$breast_summary_file"
    echo "" >> "$breast_summary_file"
    
    echo "ANCESTRY DISTRIBUTION BY STATUS:" >> "$breast_summary_file"
    echo "Format: Count Ancestry Status" >> "$breast_summary_file"
    cut -f2,3 "$temp_joined" | sort | uniq -c | sort -nr >> "$breast_summary_file"
    echo "" >> "$breast_summary_file"
    
    echo "CASES BY ANCESTRY:" >> "$breast_summary_file"
    awk '$3=="Case" {print $2}' "$temp_joined" | sort | uniq -c | sort -nr >> "$breast_summary_file"
    echo "" >> "$breast_summary_file"
    
    echo "CONTROLS BY ANCESTRY:" >> "$breast_summary_file"
    awk '$3=="Control" {print $2}' "$temp_joined" | sort | uniq -c | sort -nr >> "$breast_summary_file"
    echo "" >> "$breast_summary_file"
    
    echo "CASE-CONTROL RATIO BY ANCESTRY:" >> "$breast_summary_file"
    echo "Ancestry  Cases  Controls  Total  Case%" >> "$breast_summary_file"
    for ancestry in $unique_ancestries; do
        cases_count=$(awk -v anc="$ancestry" '$2==anc && $3=="Case"' "$temp_joined" | wc -l)
        controls_count=$(awk -v anc="$ancestry" '$2==anc && $3=="Control"' "$temp_joined" | wc -l)
        total_count=$((cases_count + controls_count))
        if [ $total_count -gt 0 ]; then
            case_percent=$(echo "scale=1; $cases_count * 100 / $total_count" | bc -l 2>/dev/null || echo "0")
            printf "%-8s  %5d  %8d  %5d  %5s%%\n" "$ancestry" "$cases_count" "$controls_count" "$total_count" "$case_percent" >> "$breast_summary_file"
        fi
    done
    
    # Clean up temp file
    rm "$temp_joined"
    
    echo "Breast cancer case-control analysis complete."
    echo "Results saved to: $breast_summary_file"
    else
    echo "Case-control file not found: $CASE_CONTROL_FILE"
    echo "Skipping detailed breast cancer ancestry analysis."
fi
