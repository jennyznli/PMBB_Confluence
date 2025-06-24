#!/bin/bash
#BSUB -J eur_qc
#BSUB -o ../log/eur_qc.%J.out
#BSUB -e ../log/eur_qc.%J.err
#BSUB -n 24
#BSUB -M 32GB
#BSUB -R "rusage[mem=32GB]"
#BSUB -W 04:00

module load plink/1.9-20210416

set -e 
WDIR="/project/knathans_tecac/jenny/breast"
HOME_DIR="/project/knathans_tecac/jenny/breast/data/EUR/step1"

# Create directory if it doesn't exist
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
QC_LOG="$HOME_DIR/eur_qc_summary.tsv"
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
      --out breast_eur_raw

# Track after ancestry filtering
track_counts "eur_ancestry" "breast_eur_raw"

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
}' breast_eur_raw.fam > breast_eur_phenotypes.fam
cp breast_eur_raw.bim breast_eur_phenotypes.bim
cp breast_eur_raw.bed breast_eur_phenotypes.bed

echo "Updated phenotypes in FAM file"

### QC STEP 1: MISSINGNESS > 0.05 ###
echo "=== QC STEP 1: MISSINGNESS FILTERING ==="

# Variant missingness (geno 0.05)
echo "Filtering variants with >5% missingness..."
plink --bfile breast_eur_phenotypes \
      --geno 0.05 \
      --make-bed \
      --threads 16 \
      --out breast_eur_geno

track_counts "geno_filter" "breast_eur_geno"

# Individual missingness (mind 0.05)
echo "Filtering individuals with >5% missingness..."
plink --bfile breast_eur_geno \
      --mind 0.05 \
      --make-bed \
      --threads 16 \
      --out breast_eur_mind

track_counts "mind_filter" "breast_eur_mind"

### QC STEP 2: MAF CHECK (DOUBLE CHECK) ###
echo "=== QC STEP 2: MAF CHECK ==="

# Check current MAF distribution
plink --bfile breast_eur_mind \
      --freq \
      --threads 16 \
      --out breast_eur_mind

low_maf=$(awk 'NR > 1 && $5 < 0.05' breast_eur_mind.frq | wc -l)
echo "Variants with MAF < 0.05: $low_maf"

# Apply MAF filter if needed
if [[ $low_maf -gt 0 ]]; then
    plink --bfile breast_eur_mind \
          --maf 0.05 \
          --make-bed \
          --threads 16 \
          --out breast_eur_maf
else
    cp breast_eur_mind.bed breast_eur_maf.bed
    cp breast_eur_mind.bim breast_eur_maf.bim
    cp breast_eur_mind.fam breast_eur_maf.fam
fi

track_counts "maf_filter" "breast_eur_maf"

### QC STEP 3: HWE CHECK ###
echo "=== QC STEP 3: HARDY-WEINBERG EQUILIBRIUM CHECK ==="

plink --bfile breast_eur_maf \
      --hardy \
      --threads 16 \
      --out breast_eur_maf

# Count HWE violations
hwe_fail=$(awk 'NR > 1 && $9 < 1e-6' breast_eur_maf.hwe | wc -l)
echo "HWE violations (p < 1e-6): $hwe_fail"

# Apply HWE filter
if [[ $hwe_fail -gt 0 ]]; then
    plink --bfile breast_eur_maf \
          --hwe 1e-6 \
          --make-bed \
          --threads 16 \
          --out breast_eur_hwe
else
    cp breast_eur_maf.bed breast_eur_hwe.bed
    cp breast_eur_maf.bim breast_eur_hwe.bim
    cp breast_eur_maf.fam breast_eur_hwe.fam
fi

track_counts "hwe_filter" "breast_eur_hwe"

### QC STEP 4: HETEROZYGOSITY ###
echo "=== QC STEP 4: HETEROZYGOSITY CHECK ==="

plink --bfile breast_eur_hwe \
      --het \
      --threads 16 \
      --out breast_eur_hwe

# Run heterozygosity QC script
$HET_FILE breast_eur_hwe
het_fail=$(wc -l < fail-het.txt)
echo "Heterozygosity outliers: $het_fail"

# Remove heterozygosity outliers
if [[ $het_fail -gt 0 ]]; then
    plink --bfile breast_eur_hwe \
          --remove fail-het.txt \
          --make-bed \
          --threads 16 \
          --out breast_eur_het
else
    cp breast_eur_hwe.bed breast_eur_het.bed
    cp breast_eur_hwe.bim breast_eur_het.bim
    cp breast_eur_hwe.fam breast_eur_het.fam
fi

track_counts "het_filter" "breast_eur_het"

### QC STEP 5: IBD/RELATEDNESS ###
echo "=== QC STEP 5: IBD/RELATEDNESS CHECK ==="

plink --bfile breast_eur_het \
      --genome \
      --threads 24 \
      --memory 32000 \
      --out breast_eur_het

# run missing to get call rates for IBD removal 
plink --bfile breast_eur_het \
      --missing \
      --threads 24 \
      --memory 32000 \
      --out breast_eur_het

# Run IBD QC script
$IBD_FILE breast_eur_het
ibd_fail=$(wc -l < fail-ibd.txt)
echo "Related individuals to remove: $ibd_fail"

# Remove related individuals
if [[ $ibd_fail -gt 0 ]]; then
    plink --bfile breast_eur_het \
          --remove fail-ibd.txt \
          --make-bed \
          --threads 16 \
          --out breast_eur_ibd
else
    cp breast_eur_het.bed breast_eur_ibd.bed
    cp breast_eur_het.bim breast_eur_ibd.bim
    cp breast_eur_het.fam breast_eur_ibd.fam
fi

track_counts "ibd_filter" "breast_eur_ibd"

### QC STEP 6: DIFFERENTIAL MISSINGNESS ###
echo "=== QC STEP 6: DIFFERENTIAL MISSINGNESS ==="

# Test for differential missingness between cases and controls
plink --bfile breast_eur_ibd \
      --test-missing \
      --threads 16 \
      --out breast_eur_ibd

# Count significant differential missingness
diff_miss=$(awk 'NR > 1 && $5 < 1e-4' breast_eur_ibd.missing | wc -l)
echo "Differential missingness violations (p < 1e-4): $diff_miss"

# Remove variants with differential missingness
if [[ $diff_miss -gt 0 ]]; then
    awk 'NR > 1 && $5 < 1e-4 {print $2}' breast_eur_ibd.missing > fail-diff-miss.txt
    plink --bfile breast_eur_ibd \
          --exclude fail-diff-miss.txt \
          --make-bed \
          --threads 16 \
          --out breast_eur_qc_final
else
    cp breast_eur_ibd.bed breast_eur_qc_final.bed
    cp breast_eur_ibd.bim breast_eur_qc_final.bim
    cp breast_eur_ibd.fam breast_eur_qc_final.fam
    touch fail-diff-miss.txt
fi

track_counts "final" "breast_eur_qc_final"

# Create list of final QC'd individuals
awk '{print $1, $2}' breast_eur_qc_final.fam > breast_eur_qc_final_pids.txt

### CALCULATE REMOVAL STATISTICS ###

# Calculate removals at each step
samples_removed_ancestry=$((${counts[initial_samples]} - ${counts[eur_ancestry_samples]}))
variants_removed_geno=$((${variants[eur_ancestry_variants]} - ${variants[geno_filter_variants]}))
samples_removed_mind=$((${counts[geno_filter_samples]} - ${counts[mind_filter_samples]}))
variants_removed_maf=$((${variants[mind_filter_variants]} - ${variants[maf_filter_variants]}))
variants_removed_hwe=$((${variants[maf_filter_variants]} - ${variants[hwe_filter_variants]}))
samples_removed_het=$((${counts[hwe_filter_samples]} - ${counts[het_filter_samples]}))
samples_removed_ibd=$((${counts[het_filter_samples]} - ${counts[ibd_filter_samples]}))
variants_removed_diffmiss=$((${variants[ibd_filter_variants]} - ${variants[final_variants]}))

# Count cases and controls in final dataset
final_cases=$(awk '$6 == 2' breast_eur_qc_final.fam | wc -l)
final_controls=$(awk '$6 == 1' breast_eur_qc_final.fam | wc -l)

### CREATE SUMMARY LOG FILE ###

# Create header if file doesn't exist or is empty
if [[ ! -s "$QC_LOG" ]] || ! grep -q "^Dataset" "$QC_LOG"; then
    echo -e "Dataset\tDate\tInitial_Samples\tInitial_Variants\tEUR_Samples\tEUR_Removed\tGeno_Variants_Removed\tMind_Samples_Removed\tMAF_Variants_Removed\tHWE_Variants_Removed\tHet_Samples_Removed\tIBD_Samples_Removed\tDiffMiss_Variants_Removed\tFinal_Samples\tFinal_Variants\tFinal_Cases\tFinal_Controls" > "$QC_LOG"
fi

# Append summary row
echo -e "breast_eur\t$(date '+%Y-%m-%d %H:%M:%S')\t${counts[initial_samples]}\t${variants[initial_variants]}\t${counts[eur_ancestry_samples]}\t$samples_removed_ancestry\t$variants_removed_geno\t$samples_removed_mind\t$variants_removed_maf\t$variants_removed_hwe\t$samples_removed_het\t$samples_removed_ibd\t$variants_removed_diffmiss\t${counts[final_samples]}\t${variants[final_variants]}\t$final_cases\t$final_controls" >> "$QC_LOG"

### CLEANUP INTERMEDIATE FILES ###
echo "=== CLEANING UP INTERMEDIATE FILES ==="

# Remove intermediate PLINK files but keep QC failure lists and final files
rm -f breast_eur_raw.* breast_eur_phenotypes.* breast_eur_geno.* breast_eur_mind.* breast_eur_maf.* breast_eur_hwe.* breast_eur_het.* breast_eur_ibd.*
rm -f common_snps_raw.* common_snps.* all_ids_sex.txt
rm -f *.log *.nosex *.frq *.hwe *.het *.genome *.missing

echo "QC complete. Summary log: $QC_LOG"
echo "Final dataset: breast_eur_qc_final.{bed,bim,fam}"
echo "QC failure files preserved: fail-*.txt"
