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
HOME_DIR="/project/knathans_tecac/jenny/breast/common_snps"
cd "$HOME_DIR"

PMBB_DIR="/static/PMBB/PMBB-Release-2024-3.0"
PHE_DIR="/static/PMBB/PMBB-Release-2024-3.0/Phenotypes/3.0"
IMP_DIR="/static/PMBB/PMBB-Release-2024-3.0/Imputed"
COM_DIR="$IMP_DIR/common_snps_LD_pruned"

PER="$PHE_DIR/PMBB-Release-2024-3.0_phenotype_person.txt"
SNP_FILE="$COM_DIR/PMBB-Release-2024-3.0_genetic_imputed.commonsnps"

# Define case/control file paths
CASES_FILE="/project/knathans_tecac/jenny/breast/cases.txt"
CONTROLS_FILE="/project/knathans_tecac/jenny/breast/controls.txt"
HET_FILE="/project/knathans_tecac/jenny/breast/het.sh"
IBD_FILE="/project/knathans_tecac/jenny/breast/ibd.sh"

### CHECK REQUIRED FILES ###
required_files=("$PER" "$SNP_FILE.bed" "$SNP_FILE.bim" "$SNP_FILE.fam")
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
    # Update FAM file sex column (column 5)
    if ($2 in sex) {
        $5 = sex[$2]
    } else {
        $5 = 0  # Unknown sex
    }
    print $1, $2, $3, $4, $5, $6
}' all_ids_sex.txt $SNP_FILE.fam > common_snps.fam
cp common_snps_raw.bim common_snps.bim
cp common_snps_raw.bed common_snps.bed

### CREATE CASE/CONTROL DATASET ###
echo "=== CREATING CASE/CONTROL DATASET ==="

# Combine cases and controls into one file
cat "$CASES_FILE" "$CONTROLS_FILE" > case_control_all_ids.txt

# Create keep file for PLINK (FID IID format)
awk '{print 0, $1}' case_control_all.txt > case_control_all_pids.txt

echo "Cases: $(wc -l < "$CASES_FILE")"
echo "Controls: $(wc -l < "$CONTROLS_FILE")"
echo "Total case/control samples: $(wc -l < case_control_all_pids.txt)"

# Subset to only cases and controls
plink --bfile common_snps \
      --keep case_control_all_pids.txt \
      --make-bed \
      --threads 16 \
      --out breast_all

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
}' breast_all.fam > breast_all_updated.fam

# Replace the FAM file
mv breast_all_updated.fam breast_all.fam

echo "Updated phenotypes in FAM file"

### INDIVIDUAL QC ###
echo "=== INDIVIDUAL QC CHECKS ==="

# Check individual missingness
echo "Checking individual missingness..."
plink --bfile breast_all \
      --missing \
      --threads 16 \
      --out breast_all

awk 'NR > 1 && $6 > 0.05 {print $1, $2}' breast_all.imiss > fail-missing-ind.txt
echo "Individuals with >5% missingness: $(wc -l < fail-missing-ind.txt)"

# Check heterozygosity
echo "Checking heterozygosity..."
plink --bfile breast_all \
      --het \
      --threads 16 \
      --out breast_all

# Run heterozygosity QC script
$HET_FILE breast_all
echo "Heterozygosity outliers: $(wc -l < fail-het.txt)"

### IBD CHECK ###
echo "=== IBD/RELATEDNESS CHECK ==="
echo "Computing pairwise IBD (this may take a while)..."

# Increase memory and threads for genome calculation
plink --bfile breast_all \
      --genome \
      --threads 24 \
      --memory 32000 \
      --out breast_all

# Run IBD QC script
$IBD_FILE breast_all
echo "Related individuals to remove: $(wc -l < fail-ibd.txt)"

### COMBINE QC FAILURES ###
echo "=== COMBINING QC FAILURES ==="

# Combine all QC failures
cat fail-missing-ind.txt fail-het.txt fail-ibd.txt > fail-qc-inds.txt

echo "Total individuals failing QC: $(wc -l < fail-qc-inds.txt)"

# Create clean dataset
plink --bfile breast_all \
      --remove fail-qc-inds.txt \
      --threads 16 \
      --make-bed \
      --out breast_qc_ind

awk '{print $1, $2}' breast_qc_ind.fam > breast_qc_ind_pids.txt

echo ""
echo "=== QC SUMMARY ==="
echo "Original samples: $(wc -l < breast_all.fam)"
echo "Failed missingness: $(wc -l < fail-missing-ind.txt)"
echo "Failed heterozygosity: $(wc -l < fail-het.txt)"
echo "Failed IBD: $(wc -l < fail-ibd.txt)"
echo "Total failed QC: $(wc -l < fail-qc-inds.txt)"
echo "Clean samples: $(wc -l < breast_qc_ind.fam)"
