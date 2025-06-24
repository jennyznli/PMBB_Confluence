#!/bin/bash

# Script to create ancestry-specific ID files from grafanc_extract.txt
# Creates files for all PMBB samples (not just breast cancer)

HOME_DIR="/project/knathans_tecac/jenny/breast"
GRAF_DIR="/project/knathans_tecac/jenny/breast/analysis/grafanc_results"
INPUT_FILE="$GRAF_DIR/grafanc_extract.txt"
OUTPUT_DIR="$HOME_DIR/analysis/ancestry"

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

# Process each ancestry group
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
    echo ""
done

# Display results
echo "=== SUMMARY ==="
echo "Total samples processed: $total_samples"
echo ""
echo "Sample counts by ancestry:"
for ancestry in $unique_ancestries; do
    count=$(wc -l < "$OUTPUT_DIR/${ancestry}_ids.txt")
    printf "%-5s: %d samples\n" "$ancestry" "$count"
done

echo ""
echo "Files created in $OUTPUT_DIR:"

echo "Summary saved to: $summary_file"
echo ""

# Add breast cancer case-control ancestry analysis
echo ""
echo "=== BREAST CANCER CASE-CONTROL ANCESTRY ANALYSIS ==="

CASE_CONTROL_FILE="$HOME_DIR/data/case_control.txt"
if [ -f "$CASE_CONTROL_FILE" ]; then
    echo "Analyzing ancestry distribution in breast cancer cases and controls..."
    
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
    
    # Display on screen
    echo ""
    echo "Breast cancer sample counts:"
    echo "Total: $total_breast_samples (Cases: $cases_total, Controls: $controls_total)"
    echo ""
    echo "Cases by ancestry:"
    awk '$3=="Case" {print $2}' "$temp_joined" | sort | uniq -c | sort -nr
    echo ""
    echo "Controls by ancestry:"
    awk '$3=="Control" {print $2}' "$temp_joined" | sort | uniq -c | sort -nr
    echo ""
    echo "Case-control ratios:"
    echo "Ancestry  Cases  Controls  Total  Case%"
    for ancestry in $unique_ancestries; do
        cases_count=$(awk -v anc="$ancestry" '$2==anc && $3=="Case"' "$temp_joined" | wc -l)
        controls_count=$(awk -v anc="$ancestry" '$2==anc && $3=="Control"' "$temp_joined" | wc -l)
        total_count=$((cases_count + controls_count))
        if [ $total_count -gt 0 ]; then
            case_percent=$(echo "scale=1; $cases_count * 100 / $total_count" | bc -l 2>/dev/null || echo "0")
            printf "%-8s  %5d  %8d  %5d  %5s%%\n" "$ancestry" "$cases_count" "$controls_count" "$total_count" "$case_percent"
        fi
    done
    
    # Clean up temp file
    rm "$temp_joined"
    
    echo ""
    echo "Breast cancer ancestry analysis saved to: $breast_summary_file"
    else
    echo "Case-control file not found: $CASE_CONTROL_FILE"
    echo "Skipping breast cancer ancestry analysis."
    echo "Run the breast cancer ancestry script first to create case_control.txt"
fi
