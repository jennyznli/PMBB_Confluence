#!/bin/bash
module load plink/2.0


txt_files=(breast/ancestry/*.txt)

# Process each ancestry-specific subset
for sample_file in "${txt_files[@]}"; do
    ANC_GROUP=$(basename "$sample_file" .txt)
    INPUT_PREFIX="${ANC_GROUP}/${ANC_GROUP}"

    echo "Processing ${ANC_GROUP}..."


        echo "Successfully ran REGENIE step 1 for ${ANC_GROUP}"
    else
        echo "Error processing ${ANC_GROUP}" >&2
        continue  # skip to next group on failure
    fi
done

