
#!/bin/bash

##########  modules  ##########
module load plink/1.9-20210416
module load htslib/1.9
module load bcftools/1.20

##########  directories  ##########
PMBB_DIR="/static/PMBB/PMBB-Release-2024-3.0/Imputed/chunked_bed_files"
HOME_DIR="/home/jennyzli/breast/Imputed"
cd "$HOME_DIR"

FILE="$PMBB_DIR/PMBB-Release-2024-3.0_genetic_imputed.chrX_chunk5_12151110_15441771"

plink --bfile "$FILE" \
      --make-bed \
      --out test_chunk

