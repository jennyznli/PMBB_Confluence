#!/bin/bash
module load plink/1.9-20210416
module load htslib/1.9

HOME_DIR="/project/knathans_tecac/jenny/breast"
PMBB_DIR="/static/PMBB/PMBB-Release-2024-3.0"
PHE_DIR="/static/PMBB/PMBB-Release-2024-3.0/Phenotypes/3.0"
PER="$PHE_DIR/PMBB-Release-2024-3.0_phenotype_person.txt"
COV="$PHE_DIR/PMBB-Release-2024-3.0_covariates.txt"
GRAF_DIR="/project/knathans_tecac/jenny/breast/grafanc_results"
GRAF="/project/knathans_tecac/jenny/breast/grafanc_results/all_pmbb_ancestry_pops.txt"

cd $GRAF_DIR

# Define case/control file paths
CASES_FILE="/project/knathans_tecac/jenny/breast/cases.txt"
CONTROLS_FILE="/project/knathans_tecac/jenny/breast/controls.txt"

echo "=== BREAST CANCER CASE-CONTROL ANCESTRY ANALYSIS ==="

# Create combined case-control sample list
echo "Creating breast cancer case-control sample list..."
echo -e "Sample\tStatus" > $HOME_DIR/case_control.txt
awk '{print $1"\tCase"}' "$CASES_FILE" >> $HOME_DIR/case_control.txt
awk '{print $1"\tControl"}' "$CONTROLS_FILE" >> $HOME_DIR/case_control.txt

# Get just sample IDs for filtering
cut -f1 $HOME_DIR/case_control.txt | tail -n +2 > breast_samples.txt

echo "Sample counts:"
echo "Cases: $(wc -l < "$CASES_FILE")"
echo "Controls: $(wc -l < "$CONTROLS_FILE")"
echo "Total breast cancer samples: $(wc -l < breast_samples.txt)"

echo "Extracting ancestry data from each file (breast cancer samples only)..."

# Extract from GRAF file and subset to breast cancer samples
echo -e "Sample\tAncGroupID" > grafanc_raw_breast.txt
cut -f1,27 "$GRAF" | tail -n +2 | grep -F -f breast_samples.txt >> grafanc_raw_breast.txt

echo "GRAF breast cancer samples found: $(tail -n +2 grafanc_raw_breast.txt | wc -l)"

# Map AncGroupID to continental and subcontinental ancestry using the mapping file
echo "Mapping GRAF ancestry codes..."
echo -e "Sample\tGrafancAncGroupID\tGrafancCont\tGrafancSub" > grafanc_extract_breast.txt

awk -F"[ \t]+" 'BEGIN {OFS="\t"}
    NR==FNR && FNR>1 {
        gsub(/^[ \t]+|[ \t]+$/, "", $1)
        map[$1] = $2 "\t" $3
        next
    }
    NR>FNR && FNR>1 {
        gsub(/^[ \t]+|[ \t]+$/, "", $2)
        if ($2 in map) {
            print $1, $2, map[$2]
        } else {
            print $1, $2, "UNKNOWN", "UNKNOWN"
        }
    }
' grafanc_mapping.txt grafanc_raw_breast.txt >> grafanc_extract_breast.txt

# Extract from covariates file and subset to breast cancer samples
echo -e "person_id\tCovClass" > cov_extract_breast.txt
cut -f1,5 "$COV" | tail -n +2 | grep -F -f breast_samples.txt >> cov_extract_breast.txt

echo "Covariate breast cancer samples found: $(tail -n +2 cov_extract_breast.txt | wc -l)"

# Extract from person file and subset to breast cancer samples
echo -e "person_id\trace_concept_id\trace_source_value\tethnicity_concept_id\tethnicity_source_value" > person_extract_breast.txt
cut -f1,7,15,8,17 "$PER" | tail -n +2 | grep -F -f breast_samples.txt >> person_extract_breast.txt

echo "Person breast cancer samples found: $(tail -n +2 person_extract_breast.txt | wc -l)"

# Get unique sample IDs from grafanc_breast, which has less samples than PMBB files
cut -f1 grafanc_extract_breast.txt | tail -n +2 | sort > samples_graf_breast.txt

# Use grafanc_extract_breast as the base and add covariates data
echo "Adding covariates data..."
awk -F'\t' 'BEGIN {OFS="\t"}
    # Read covariates file into array
    NR==FNR && NR>1 {
        cov[$1] = $2
        next
    }
    # Process grafanc_extract_breast
    NR>FNR {
        if (FNR == 1) {
            print $0, "CovAnc"
        } else {
            if ($1 in cov) {
                print $0, cov[$1]
            } else {
                print $0, "NA"
            }
        }
    }' cov_extract_breast.txt grafanc_extract_breast.txt > temp_with_cov_breast.txt

echo "Adding person data (race and ethnicity)..."

awk -F'\t' 'BEGIN {OFS="\t"}
    # Read person file into array
    FNR > 1 && NR == FNR {
        person[$1] = $2 "\t" $3 "\t" $4 "\t" $5
        next
    }
    # Process temp_with_cov_breast.txt
    NR > FNR {
        if (FNR == 1) {
            # Handle header line
            print $0, "RaceID", "RaceValue", "EthnicityID", "EthnicityValue"
        } else {
            if ($1 in person) {
                print $0, person[$1]
            } else {
                print $0, "NA", "NA", "NA", "NA"
            }
        }
    }
' person_extract_breast.txt temp_with_cov_breast.txt > temp_ancestry_breast.txt

# Add case-control status to final file
echo "Adding case-control status..."
awk -F'\t' 'BEGIN {OFS="\t"}
    # Read case-control status
    NR==FNR && NR>1 {
        status[$1] = $2
        next
    }
    # Process ancestry file and add status
    NR>FNR {
        if (FNR == 1) {
            print $0, "Status"
        } else {
            if ($1 in status) {
                print $0, status[$1]
            } else {
                print $0, "Unknown"
            }
        }
    }
' $HOME_DIR/case_control.txt temp_ancestry_breast.txt > ancestry_breast_final.txt

# Clean up temp files
rm temp_with_cov_breast.txt temp_ancestry_breast.txt

INPUT_FILE="ancestry_breast_final.txt"
OUTPUT_FILE="ancestry_breast_analysis_summary.txt"

echo "=== BREAST CANCER ANCESTRY COMPARISON ANALYSIS ===" > $OUTPUT_FILE
echo "Analysis of: $INPUT_FILE" >> $OUTPUT_FILE
echo "Generated on: $(date)" >> $OUTPUT_FILE
echo "" >> $OUTPUT_FILE

# Basic counts by case-control status
echo "=== SAMPLE COUNTS ===" >> $OUTPUT_FILE
total_samples=$(tail -n +2 $INPUT_FILE | wc -l)
cases_count=$(awk -F'\t' 'NR>1 && $10=="Case"' $INPUT_FILE | wc -l)
controls_count=$(awk -F'\t' 'NR>1 && $10=="Control"' $INPUT_FILE | wc -l)

echo "Total breast cancer samples: $total_samples" >> $OUTPUT_FILE
echo "Cases: $cases_count" >> $OUTPUT_FILE
echo "Controls: $controls_count" >> $OUTPUT_FILE
echo "" >> $OUTPUT_FILE

# 1. GRAF Continental vs Covariates Ancestry
echo "=== GRAF CONTINENTAL vs COVARIATES ANCESTRY ===" >> $OUTPUT_FILE
echo "GrafancCont -> CovAnc comparison:" >> $OUTPUT_FILE
cut -f3,5 $INPUT_FILE | tail -n +2 | sort | uniq -c | sort -nr >> $OUTPUT_FILE
echo "" >> $OUTPUT_FILE

# Agreement rate between GRAF and Covariates
graf_cov_agree=$(cut -f3,5 $INPUT_FILE | tail -n +2 | awk '$1==$2' | wc -l)
graf_cov_percent=$(echo "scale=1; $graf_cov_agree * 100 / $total_samples" | bc)
echo "GRAF-Covariate Agreement: $graf_cov_agree/$total_samples (${graf_cov_percent}%)" >> $OUTPUT_FILE
echo "" >> $OUTPUT_FILE

# 2. GRAF Continental vs Self-Reported Race
echo "=== GRAF CONTINENTAL vs SELF-REPORTED RACE ===" >> $OUTPUT_FILE
echo "GrafancCont -> RaceValue comparison:" >> $OUTPUT_FILE
cut -f3,8 $INPUT_FILE | tail -n +2 | sort | uniq -c | sort -nr >> $OUTPUT_FILE
echo "" >> $OUTPUT_FILE

# 3. GRAF Subcontinental vs Self-Reported Race
echo "=== GRAF SUBCONTINENTAL vs SELF-REPORTED RACE ===" >> $OUTPUT_FILE
echo "GrafancSub -> RaceValue comparison:" >> $OUTPUT_FILE
cut -f4,8 $INPUT_FILE | tail -n +2 | sort | uniq -c | sort -nr >> $OUTPUT_FILE
echo "" >> $OUTPUT_FILE

# 4. Covariates vs Self-Reported Race
echo "=== COVARIATES vs SELF-REPORTED RACE ===" >> $OUTPUT_FILE
echo "CovAnc -> RaceValue comparison:" >> $OUTPUT_FILE
cut -f5,8 $INPUT_FILE | tail -n +2 | sort | uniq -c | sort -nr >> $OUTPUT_FILE
echo "" >> $OUTPUT_FILE

# 5. Ancestry by Case-Control Status
echo "=== ANCESTRY BY CASE-CONTROL STATUS ===" >> $OUTPUT_FILE
echo "GRAF Continental by Status:" >> $OUTPUT_FILE
cut -f3,10 $INPUT_FILE | tail -n +2 | sort | uniq -c | sort -nr >> $OUTPUT_FILE
echo "" >> $OUTPUT_FILE

echo "Self-Reported Race by Status:" >> $OUTPUT_FILE
cut -f8,10 $INPUT_FILE | tail -n +2 | sort | uniq -c | sort -nr >> $OUTPUT_FILE
echo "" >> $OUTPUT_FILE

# 6. Individual ancestry category breakdowns
echo "=== INDIVIDUAL CATEGORY BREAKDOWNS ===" >> $OUTPUT_FILE

echo "GRAF Continental Categories:" >> $OUTPUT_FILE
cut -f3 $INPUT_FILE | tail -n +2 | sort | uniq -c | sort -nr >> $OUTPUT_FILE
echo "" >> $OUTPUT_FILE

echo "GRAF Subcontinental Categories:" >> $OUTPUT_FILE
cut -f4 $INPUT_FILE | tail -n +2 | sort | uniq -c | sort -nr >> $OUTPUT_FILE
echo "" >> $OUTPUT_FILE

echo "Covariates Ancestry Categories:" >> $OUTPUT_FILE
cut -f5 $INPUT_FILE | tail -n +2 | sort | uniq -c | sort -nr >> $OUTPUT_FILE
echo "" >> $OUTPUT_FILE

echo "Self-Reported Race Categories:" >> $OUTPUT_FILE
cut -f8 $INPUT_FILE | tail -n +2 | sort | uniq -c | sort -nr >> $OUTPUT_FILE
echo "" >> $OUTPUT_FILE

echo "Self-Reported Ethnicity Categories:" >> $OUTPUT_FILE
cut -f9 $INPUT_FILE | tail -n +2 | sort | uniq -c | sort -nr >> $OUTPUT_FILE
echo "" >> $OUTPUT_FILE

# 7. Discordant cases analysis
echo "=== DISCORDANT CASES ANALYSIS ===" >> $OUTPUT_FILE

echo "Breast cancer cases where GRAF != Covariates:" >> $OUTPUT_FILE
awk -F'\t' 'NR>1 && $3!=$5 {print $1"\t"$3"\t"$5"\t"$8"\t"$10}' $INPUT_FILE > discord_graf_cov_breast.txt
awk -F'\t' 'NR>1 && $3!=$5 {print $1"\t"$3"\t"$5"\t"$8"\t"$10}' $INPUT_FILE | head -20 >> $OUTPUT_FILE
echo "" >> $OUTPUT_FILE

echo "Cases where genetic ancestry != self-reported race:" >> $OUTPUT_FILE
echo "Sample\tGRAF\tCov\tSelf-Reported\tStatus" >> $OUTPUT_FILE
awk -F'\t' 'NR>1 {
    if (($3=="EUR" && $8!="White" && $7!="White/Caucasian") ||
        ($3=="AFR" && $8!="Black or African American") ||
        ($3=="AMR" && ($8!~"Hispanic" && $8!~"Latino"))) {
        print $1"\t"$3"\t"$5"\t"$8"\t"$10
    }
}' $INPUT_FILE | head -20 >> $OUTPUT_FILE
echo "" >> $OUTPUT_FILE

# 8. Summary statistics
echo "=== SUMMARY STATISTICS ===" >> $OUTPUT_FILE
echo "Total breast cancer samples analyzed: $total_samples" >> $OUTPUT_FILE
echo "Cases: $cases_count" >> $OUTPUT_FILE
echo "Controls: $controls_count" >> $OUTPUT_FILE
echo "GRAF-Covariate agreement rate: ${graf_cov_percent}%" >> $OUTPUT_FILE

# Count major discrepancies
major_discrepancies=$(awk -F'\t' 'NR>1 && $3!=$5' $INPUT_FILE | wc -l)
echo "Major ancestry discrepancies (GRAF≠Cov): $major_discrepancies" >> $OUTPUT_FILE

echo "" >> $OUTPUT_FILE
echo "Analysis complete! Check $OUTPUT_FILE for detailed results." >> $OUTPUT_FILE

# Clean up intermediate files
rm grafanc_raw_breast.txt samples_graf_breast.txt breast_samples.txt 2>/dev/null || true

# Display summary to screen
echo "=== BREAST CANCER ANALYSIS COMPLETE ==="
echo "Results saved to: $OUTPUT_FILE"
echo ""
echo "Quick Summary:"
echo "- Total breast cancer samples: $total_samples"
echo "- Cases: $cases_count, Controls: $controls_count"
echo "- GRAF-Covariate agreement: $graf_cov_agree/$total_samples (${graf_cov_percent}%)"
echo "- Major discrepancies: $major_discrepancies"
echo ""
echo "Output files:"
echo "- Main results: ancestry_breast_final.txt"
echo "- Summary: ancestry_breast_analysis_summary.txt"
echo "- Case-control list: $HOME_DIR/case_control.txt"
echo "- Discordant samples: discord_graf_cov_breast.txt"
echo ""
echo "Final columns:"
echo "1. Sample ID, 2. AncGroupID, 3. GrafancCont, 4. GrafancSub"
echo "5. CovAnc, 6. RaceID, 7. RaceValue, 8. EthnicityID"
echo "9. EthnicityValue, 10. Status (Case/Control)"
