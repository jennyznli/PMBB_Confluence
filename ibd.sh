#!/bin/bash
# Usage: ./ibd.sh "prefix"

INPUT=$1

# Check if files exist
if [[ ! -f "$INPUT.genome" ]]; then
    echo "Error: $INPUT.genome not found"
    exit 1
fi
if [[ ! -f "$INPUT.imiss" ]]; then
    echo "Error: $INPUT.imiss not found" 
    exit 1
fi

# Step 1: Filter genome file for related pairs
awk 'NR>1 && $10 > 0.185 {print $1, $2, $3, $4, $10}' "$INPUT.genome" > related-pairs.txt

# Step 2: Compare call rates and remove person with lower call rate 
awk -v INPUT="$INPUT" '
BEGIN {
    # Load call rates from imiss file
    while ((getline line < (INPUT ".imiss")) > 0) {
        split(line, fields)
        if (fields[1] == "FID") continue
        fid = fields[1]; iid = fields[2]; f_miss = fields[6]
        call_rate[fid"_"iid] = 1 - f_miss
    }
    close(INPUT ".imiss")
}
{
    id1 = $1"_"$2
    id2 = $3"_"$4
    if (id1 in call_rate && id2 in call_rate) {
        if (call_rate[id1] < call_rate[id2]) {
            print $1, $2  # Remove person 1
        } else {
            print $3, $4  # Remove person 2
        }
    }
}
' related-pairs.txt | sort | uniq > fail-ibd.txt

