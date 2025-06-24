#!/bin/bash
#BSUB -J commonsnps_qc
#BSUB -o log/commonsnps_qc.%J.out
#BSUB -e log/commonsnps_qc.%J.err
#BSUB -n 24
#BSUB -M 32GB
#BSUB -R "rusage[mem=32GB]"
#BSUB -W 04:00

module load plink/1.9-20210416
module load htslib/1.9

set -e 
WDIR="/project/knathans_tecac/jenny/breast"
HOME_DIR="/project/knathans_tecac/jenny/breast/data/imputed/EUR/step1"
mkdir -p "$HOME_DIR"
cd "$HOME_DIR"

PMBB_DIR="/static/PMBB/PMBB-Release-2024-3.0"
PHE_DIR="/static/PMBB/PMBB-Release-2024-3.0/Phenotypes/3.0"
IMP_DIR="/static/PMBB/PMBB-Release-2024-3.0/Imputed"
COM_DIR="$IMP_DIR/common_snps_LD_pruned"
PER="$PHE_DIR/PMBB-Release-2024-3.0_phenotype_person.txt"
SNP_FILE="$COM_DIR/PMBB-Release-2024-3.0_genetic_imputed.commonsnps"

EUR_FILE="$WDIR/analysis/ancestry/EUR_breast_pids.txt"

CASES_FILE="$WDIR/data/cases.txt"
CONTROLS_FILE="$WDIR/data/controls.txt"
HET_FILE="$WDIR/scripts/het.sh"
IBD_FILE="$WDIR/scripts/ibd.sh"

# Initialize QC tracking log
QC_LOG="$WDIR/log/eur_qc_summary.tsv"
echo "Starting QC at $(date)" > "$QC_LOG"

# Initialize tracking variables
declare -A counts
declare -A variants

### TRACKING FUNCTION ###
track_counts() {
    local step="$1"
    local bfile="$2"
    
    if [[ -f "${bfile}.fam" && -f "${bfile}.bim" ]]; then
        counts["${step}_samples"]=$(wc -l < "${bfile}.fam")
        variants["${step}_variants"]=$(wc -l < "${bfile}.bim")
    fi
}

### CHECK REQUIRED FILES ###
required_files=("$PER" "$SNP_FILE.bed" "$SNP_FILE.bim" "$SNP_FILE.fam" "$EUR_FILE")
for file in "${required_files[@]}"; do
    [[ ! -f "$file" ]] && { echo "ERROR: Required file not found: $file"; exit 1; }
done

# Check for case/control files
if [[ ! -f "$CASES_FILE" ]] || [[ ! -f "$CONTROLS_FILE" ]]; then
    echo "ERROR: cases.txt and controls.txt files required"
    echo "Expected locations:"
    echo "  Cases: $CASES_FILE"
    echo "  Controls: $CONTROLS_FILE"
    exit 1
fi

### PREP SEX INFO FOR QC ###
echo "=== PREPARING SEX INFORMATION ==="

# Copy original files with raw suffix
cp $SNP_FILE.bim common_snps_raw.bim
cp $SNP_FILE.bed common_snps_raw.bed
cp $SNP_FILE.fam common_snps_raw.fam

# Track initial counts
track_counts "initial" "common_snps_raw"

# Get sex information from phenotypes (ID and sex code)
awk 'NR > 1 {print $1, $2}' "$PER" > all_ids_sex.txt

echo "Sex codes found:"
awk '{print $2}' all_ids_sex.txt | sort | uniq -c

# Update FAM file with correct sex information
# 8532 = Female (2), 8507 = Male (1)
awk '
BEGIN { OFS="\t" }
NR==FNR {
    # Read sex information
    sex_code = ($2 == 8507 ? 1 : ($2 == 8532 ? 2 : 0));
    sex[$1] = sex_code;
    next
}
{
    if ($2 in sex) {
        $5 = sex[$2]
    } else {
        $5 = 0  # Unknown sex
    }
    print $1, $2, $3, $4, $5, $6
}' all_ids_sex.txt $SNP_FILE.fam > common_snps.fam
cp common_snps_raw.bim common_snps.bim
cp common_snps_raw.bed common_snps.bed

### CREATE EUR ANCESTRY DATASET ###
echo "=== CREATING EUR ANCESTRY DATASET ==="

# Subset to only EUR ancestry individuals
plink --bfile common_snps \
      --keep $EUR_FILE \
      --make-bed \
      --threads 16 \
      --out breast_eur

# Track after ancestry filtering
track_counts "eur_ancestry" "breast_eur"

### UPDATE PHENOTYPE INFORMATION ###
echo "=== UPDATING PHENOTYPE INFORMATION ==="

# Update phenotype column (column 6) in FAM file
awk '
BEGIN {
    # Read cases file
    while ((getline id < "'$CASES_FILE'") > 0) {
        cases[id] = 1
    }
    close("'$CASES_FILE'")
    # Read controls file
    while ((getline id < "'$CONTROLS_FILE'") > 0) {
        controls[id] = 1
    }
    close("'$CONTROLS_FILE'")
}
{
    # Update phenotype based on case/control status
    if ($2 in cases) {
        $6 = 2  # Case
    } else if ($2 in controls) {
        $6 = 1  # Control
    } else {
        $6 = 0  # Unknown/Missing
    }
    print $1, $2, $3, $4, $5, $6
}' breast_eur.fam > breast_eur_updated.fam

# Replace the FAM file
mv breast_eur_updated.fam breast_eur.fam

echo "Updated phenotypes in FAM file"

### VARIANT QC ###
echo "=== VARIANT QC CHECKS ==="

# maf > 0.05 - this was already done but just double checking 
echo "Checking minor allele frequency..."
plink --bfile breast_eur \
      --maf 0.05 \
      --make-bed \
      --threads 16 \
      --out breast_eur_maf

# Track after MAF filtering
track_counts "maf_filter" "breast_eur_maf"

echo "Checking variant missingness..."
plink --bfile breast_eur_maf \
      --geno 0.05 \
      --make-bed \
      --threads 16 \
      --out breast_eur_variants_qc

# Track after variant missingness filtering
track_counts "variant_miss" "breast_eur_variants_qc"

# hwe 1e-6 was already done by PMBB
rm breast_eur_maf*

### INDIVIDUAL QC ###
echo "=== INDIVIDUAL QC CHECKS ==="

# Check individual missingness
echo "Checking individual missingness..."
plink --bfile breast_eur_variants_qc \
      --missing \
      --threads 16 \
      --out breast_eur_variants_qc

awk 'NR > 1 && $6 > 0.05 {print $1, $2}' breast_eur_variants_qc.imiss > fail-missing-ind.txt
missing_fail=$(wc -l < fail-missing-ind.txt)
echo "Individuals with >5% missingness: $missing_fail"

# Check heterozygosity
echo "Checking heterozygosity..."
plink --bfile breast_eur_variants_qc \
      --het \
      --threads 16 \
      --out breast_eur_variants_qc

# Run heterozygosity QC script
$HET_FILE breast_eur_variants_qc
het_fail=$(wc -l < fail-het.txt)
echo "Heterozygosity outliers: $het_fail"

### IBD CHECK ###
echo "=== IBD/RELATEDNESS CHECK ==="
echo "Computing pairwise IBD..."

# Increase memory and threads for genome calculation
plink --bfile breast_eur_variants_qc \
      --genome \
      --threads 24 \
      --memory 32000 \
      --out breast_eur_variants_qc

# Run IBD QC script
$IBD_FILE breast_eur_variants_qc
ibd_fail=$(wc -l < fail-ibd.txt)
echo "Related individuals to remove: $ibd_fail"

### COMBINE QC FAILURES ###
echo "=== COMBINING QC FAILURES ==="

# Combine all QC failures
cat fail-missing-ind.txt fail-het.txt fail-ibd.txt | sort | uniq > fail-qc-inds.txt
total_fail=$(wc -l < fail-qc-inds.txt)

echo "Total individuals failing QC: $total_fail"

# Create clean dataset
plink --bfile breast_eur_variants_qc \
      --remove fail-qc-inds.txt \
      --threads 16 \
      --make-bed \
      --out breast_eur_qc_final

# Track final counts
track_counts "final" "breast_eur_qc_final"

# Create list of final QC'd individuals
awk '{print $1, $2}' breast_eur_qc_final.fam > breast_eur_qc_final_pids.txt

### CREATE SUMMARY LOG FILE ###
echo "=== CALCULATING REMOVAL STATISTICS ==="

# Calculate removals at each step
samples_removed_ancestry=$((${counts[initial_samples]} - ${counts[eur_ancestry_samples]}))
variants_removed_maf=$((${variants[eur_ancestry_variants]} - ${variants[maf_filter_variants]}))
variants_removed_geno=$((${variants[maf_filter_variants]} - ${variants[variant_miss_variants]}))
samples_removed_indqc=$((${counts[variant_miss_samples]} - ${counts[final_samples]}))

# Count cases and controls in final dataset
final_cases=$(awk '$6 == 2' breast_eur_qc_final.fam | wc -l)
final_controls=$(awk '$6 == 1' breast_eur_qc_final.fam | wc -l)

# Calculate removal statistics for logging
samples_removed_ancestry=$((${counts[initial_samples]} - ${counts[eur_ancestry_samples]}))
variants_removed_maf=$((${variants[eur_ancestry_variants]} - ${variants[maf_filter_variants]}))
variants_removed_geno=$((${variants[maf_filter_variants]} - ${variants[variant_miss_variants]}))
samples_removed_indqc=$((${counts[variant_miss_samples]} - ${counts[final_samples]}))

# Count cases and controls in final dataset
final_cases=$(awk '$6 == 2' breast_eur_qc_final.fam | wc -l)
final_controls=$(awk '$6 == 1' breast_eur_qc_final.fam | wc -l)

### CREATE SUMMARY LOG FILE ###

# Create header if file doesn't exist or is empty
if [[ ! -s "$QC_LOG" ]] || ! grep -q "^Dataset" "$QC_LOG"; then
    echo -e "Dataset\tDate\tInitial_Samples\tInitial_Variants\tEUR_Samples\tEUR_Removed\tMAF_Variants\tMAF_Removed\tGeno_Variants\tGeno_Removed\tMissing_Fail\tHet_Fail\tIBD_Fail\tTotal_Ind_Fail\tFinal_Samples\tFinal_Variants\tFinal_Cases\tFinal_Controls" > "$QC_LOG"
fi

# Append summary row
echo -e "breast_eur\t$(date '+%Y-%m-%d %H:%M:%S')\t${counts[initial_samples]}\t${variants[initial_variants]}\t${counts[eur_ancestry_samples]}\t$samples_removed_ancestry\t${variants[maf_filter_variants]}\t$variants_removed_maf\t${variants[variant_miss_variants]}\t$variants_removed_geno\t$missing_fail\t$het_fail\t$ibd_fail\t$total_fail\t${counts[final_samples]}\t${variants[final_variants]}\t$final_cases\t$final_controls" >> "$QC_LOG"

echo "QC complete. Summary log: $QC_LOG"
