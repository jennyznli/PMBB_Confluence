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

echo "Extracting ancestry data from each file..."

# Extract from GRAF file (columns 1 and 27: Sample ID and AncGroupID)
echo -e "Sample\tAncGroupID" > grafanc_raw.txt
cut -f1,27 "$GRAF" | tail -n +2 >> grafanc_raw.txt

# Map AncGroupID to continental and subcontinental ancestry using the mapping file
echo "Mapping GRAF ancestry codes..."
echo -e "Sample\tGrafancAncGroupID\tGrafancCont\tGrafancSub" > grafanc_extract.txt

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
' grafanc_mapping.txt grafanc_raw.txt >> grafanc_extract.txt

# Extract from covariates file (person_id and Class - column 5)
echo -e "person_id\tCovClass" > cov_extract.txt
cut -f1,5 "$COV" | tail -n +2 >> cov_extract.txt

# Extract from person file (person_id, race_concept_id, race_source_value, ethnicity_concept_id, ethnicity_source_value)
# Columns: 1, 7, 15, 8, 17
echo -e "person_id\trace_concept_id\trace_source_value\tethnicity_concept_id\tethnicity_source_value" > person_extract.txt
cut -f1,7,15,8,17 "$PER" | tail -n +2 >> person_extract.txt

# Get unique sample IDs from grafanc, which has less samples than PMBB files
cut -f1 grafanc_extract.txt | tail -n +2 | sort > samples_graf.txt

# Use grafanc_extract as the base and add covariates data
echo "Adding covariates data..."
awk -F'\t' 'BEGIN {OFS="\t"}
    # Read covariates file into array
    NR==FNR && NR>1 {
        cov[$1] = $2
        next
    }
    # Process grafanc_extract
    NR>FNR {
        if (NR==1) {
            print $0, "CovClass"
        } else {
            if ($1 in cov) {
                print $0, cov[$1]
            } else {
                print $0, "CovAnc"
            }
        }
    }' cov_extract.txt grafanc_extract.txt > temp_with_cov.txt

echo "Adding person data (race and ethnicity)..."

awk -F'\t' 'BEGIN {OFS="\t"}
    # Read person file into array
    FNR > 1 && NR == FNR {
        person[$1] = $2 "\t" $3 "\t" $4 "\t" $5
        next
    }
    # Process temp_with_cov.txt
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
' person_extract.txt temp_with_cov.txt > ancestry_final.txt

# Clean up temp file
rm temp_with_cov.txt
rm 

INPUT_FILE="ancestry_final.txt"
OUTPUT_FILE="ancestry_analysis_summary.txt"

echo "=== ANCESTRY COMPARISON ANALYSIS ===" > $OUTPUT_FILE
echo "Analysis of: $INPUT_FILE" >> $OUTPUT_FILE
echo "Generated on: $(date)" >> $OUTPUT_FILE
echo "" >> $OUTPUT_FILE

# Basic counts
echo "=== SAMPLE COUNTS ===" >> $OUTPUT_FILE
total_samples=$(tail -n +2 $INPUT_FILE | wc -l)
echo "Total samples: $total_samples" >> $OUTPUT_FILE
echo "" >> $OUTPUT_FILE

# 1. GRAF Continental vs Covariates Ancestry
echo "=== GRAF CONTINENTAL vs COVARIATES ANCESTRY ===" >> $OUTPUT_FILE
echo "GrafancCont -> CovAnc comparison:" >> $OUTPUT_FILE
cut -f3,5 $INPUT_FILE | tail -n +2 | sort | uniq -c | sort -nr >> $OUTPUT_FILE
echo "" >> $OUTPUT_FILE

# should do more analysis with this

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

# Agreement rate between Covariates and Self-Reported Race
cov_race_agree=$(awk -F'\t' 'NR>1 {
    if (($5=="EUR" && ($7=="White" || $7=="White/Caucasian")) || 
        ($5=="AFR" && ($7=="Black or African American")) ||
        ($5=="AMR" && ($7=="Hispanic" || $7=="Latino"))) {
        count++
    }
} END {print count+0}' $INPUT_FILE)
cov_race_percent=$(echo "scale=1; $cov_race_agree * 100 / $total_samples" | bc)
echo "Covariate-Race Agreement (estimated): $cov_race_agree/$total_samples (${cov_race_percent}%)" >> $OUTPUT_FILE
echo "" >> $OUTPUT_FILE

# 5. Individual ancestry category breakdowns
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

# 6. Discordant cases analysis
echo "=== DISCORDANT CASES ANALYSIS ===" >> $OUTPUT_FILE

echo "Cases where GRAF != Covariates:" >> $OUTPUT_FILE
awk -F'\t' 'NR>1 && $3!=$5 {print $1"\t"$3"\t"$5"\t"$8}' $INPUT_FILE > discord_graf_cov.txt
awk -F'\t' 'NR>1 && $3!=$5 {print $1"\t"$3"\t"$5"\t"$8}' $INPUT_FILE | head -20 >> $OUTPUT_FILE
echo "" >> $OUTPUT_FILE


# i'm not sure if this section is correct ... 

echo "Cases where ancestry predictions != self-reported race:" >> $OUTPUT_FILE
echo "Sample\tGRAF\tCov\tSelf-Reported" >> $OUTPUT_FILE
awk -F'\t' 'NR>1 {
    if (($3=="EUR" && $7!="White" && $7!="White/Caucasian") ||
        ($3=="AFR" && $7!="Black or African American") ||
        ($3=="AMR" && $7!~"Hispanic" && $7!~"Latino")) {
        print $1"\t"$3"\t"$5"\t"$7
    }
}' $INPUT_FILE | head -20 >> $OUTPUT_FILE
echo "" >> $OUTPUT_FILE

# 7. Create cross-tabulation matrices
echo "=== CROSS-TABULATION MATRICES ===" >> $OUTPUT_FILE

echo "GRAF Continental vs Covariates Matrix:" >> $OUTPUT_FILE
# Create a simple cross-tab
cut -f3,5 $INPUT_FILE | tail -n +2 | sort | uniq -c | \
awk '{print $2,$3,$1}' | sort | \
awk '{
    matrix[$1","$2] = $3
    row[$1]++
    col[$2]++
}
END {
    # Print header
    printf "%-10s", ""
    for (c in col) printf "%-10s", c
    printf "\n"
    
    # Print rows
    for (r in row) {
        printf "%-10s", r
        for (c in col) {
            printf "%-10s", (matrix[r","c] ? matrix[r","c] : 0)
        }
        printf "\n"
    }
}' >> $OUTPUT_FILE
echo "" >> $OUTPUT_FILE

# 8. Summary statistics
echo "=== SUMMARY STATISTICS ===" >> $OUTPUT_FILE
echo "Total samples analyzed: $total_samples" >> $OUTPUT_FILE
echo "GRAF-Covariate agreement rate: ${graf_cov_percent}%" >> $OUTPUT_FILE
echo "Estimated Covariate-Race agreement rate: ${cov_race_percent}%" >> $OUTPUT_FILE

# Count major discrepancies
major_discrepancies=$(awk -F'\t' 'NR>1 && $3!=$5' $INPUT_FILE | wc -l)
echo "Major ancestry discrepancies (GRAF≠Cov): $major_discrepancies" >> $OUTPUT_FILE

echo "" >> $OUTPUT_FILE
echo "Analysis complete! Check $OUTPUT_FILE for detailed results." >> $OUTPUT_FILE

# Display summary to screen
echo "=== ANALYSIS COMPLETE ==="
echo "Results saved to: $OUTPUT_FILE"
echo ""
echo "Quick Summary:"
echo "- Total samples: $total_samples"
echo "- GRAF-Covariate agreement: $graf_cov_agree/$total_samples (${graf_cov_percent}%)"
echo "- Major discrepancies: $major_discrepancies"
echo ""
echo "First 20 lines of analysis:"
head -20 $OUTPUT_FILE
