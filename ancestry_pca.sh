#!/bin/bash
module load plink/2.0

# Performs PCA, makes covariate and phenotype files

txt_files=(breast/ancestry/*.txt)

# Process each ancestry-specific subset
for sample_file in "${txt_files[@]}"; do
    ANC_GROUP=$(basename "$sample_file" .txt)
    INPUT_PREFIX="${ANC_GROUP}/${ANC_GROUP}"

    echo "Processing ${ANC_GROUP}..."

    if plink2 --bfile "$INPUT_PREFIX" \
              --pca 10 approx \
              --out "$INPUT_PREFIX"; then
        echo "Successfully ran PCA for ${ANC_GROUP}"
    else
        echo "Error processing ${ANC_GROUP}" >&2
        continue  # skip to next group on failure
    fi

    # COVARIATE FILE
    if [[ -f "${INPUT_PREFIX}.eigenvec" ]]; then
        (echo "FID IID V1 V2 V3 V4 V5 V6 V7 V8 V9 V10"
         sed '1{/^#/d};' "${INPUT_PREFIX}.eigenvec") > "${ANC_GROUP}/step1_covariate.txt"
    else
        echo "Missing eigenvec file for ${ANC_GROUP}" >&2
        continue
    fi

    # PHENOTYPE FILE (from .fam: FID IID SEX or trait in column 6)
    (echo "FID IID Y"
     awk '{print $1, $2, $6}' "$INPUT_PREFIX.fam") > "${ANC_GROUP}/step1_phenotype.txt"

done

