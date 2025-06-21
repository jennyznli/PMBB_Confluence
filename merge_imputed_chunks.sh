#!/bin/bash
#BSUB -q knmkln_normal
#BSUB -J "QC[1-980]"                 # adjust 980 to the total number of chunks
#BSUB -o log/qc.%J_%I.out
#BSUB -e log/qc.%J_%I.err
#BSUB -n 1
#BSUB -M 32000
#BSUB -R "rusage[mem=16GB]"
#BSUB -W 8:00

##########  modules  ##########
module load plink/2.0
module load htslib/1.9
module load bcftools/1.20

##########  directories  ##########
PMBB_DIR="/static/PMBB/PMBB-Release-2024-3.0/Imputed"
PGEN_DIR="/static/PMBB/PMBB-Release-2024-3.0/Imputed/chunked_pgen_files"
HOME_DIR="/home/jennyzli/breast/step2"
OUTPUT_DIR="/home/jennyzli/breast/step2/pgen"
cd "$HOME_DIR"

##########  read manifest  ##########
input_manifest="${PMBB_DIR}/imputed_variant_chunked_input_manifest.tsv"
read chrid chunk start_pos stop_pos < <(
  awk -v i="$LSB_JOBINDEX" '$1==i {print $2, $3, $4, $5}' "$input_manifest"
)

##########  derive names  ##########
BASE="${PGEN_DIR}/PMBB-Release-2024-3.0_genetic_imputed.${chrid}_chunk${chunk}_${start_pos}_${stop_pos}"

##########  get initial counts  ##########
# n_expected_variants=$(wc -l < ${bed_prefix}.bim)
# echo "InitialVariants: $n_expected_variants" > table_${base}.txt

##########  STEP1: QC filters  ##########
# Step 1: Filter variants by MAC (Minor Allele Count) ≥ 30
plink2 --pfile ${BASE} \
    --mac 30 \
    --make-pgen \
    --out "$OUTPUT_DIR/${BASE}_mac"

# Step 2: Filter by Rsq > 0.2 (already done by PMBB) 

# Step 3: Filter variants with missingness > 5%
plink2 --pfile "$OUTPUT_DIR/${BASE}_mac" \
    --geno 0.05 \
    --make-pgen \
    --out "$OUTPUT_DIR/${BASE}_mac_snp"

# Step 4: Filter samples with > 5% missingness
plink2 --pfile "$OUTPUT_DIR/${BASE}_mac_snp" \
    --mind 0.05 \
    --make-pgen \
    --out "$OUTPUT_DIR/${BASE}_mac_snp_sample"

# HWE?? 

##########  record final count  ##########
echo "FinalVariants: $(wc -l < ${base}_qc.bim)" >> table_${base}.txt



