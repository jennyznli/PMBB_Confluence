#!/bin/bash
module load plink/2.0

# perform regenie step 1 on each ancestry group

txt_files=(breast/ancestry/*.txt)

# Process each ancestry-specific subset
for sample_file in "${txt_files[@]}"; do
    ANC_GROUP=$(basename "$sample_file" .txt)
    INPUT_PREFIX="${ANC_GROUP}/${ANC_GROUP}"

    echo "Processing ${ANC_GROUP}..."

if regenie --step 1 \
  --bed "${ANC_GROUP}/${ANC_GROUP}" \
  --phenoFile "${ANC_GROUP}/step1_phenotype.txt" \
  --covarFile "${ANC_GROUP}/step1_covariates.txt" \
  --covarColList PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
  --strict \
  --bsize 1000 --bt \
  --lowmem --lowmem-prefix tmp_rg \
  --gz --threads 8 \
  --use-relative-path \
  --out regenie_step1_out \
  --loocv
        echo "Successfully ran REGENIE step 1 for ${ANC_GROUP}"
    else
        echo "Error processing ${ANC_GROUP}" >&2
        continue  # skip to next group on failure
    fi
done

