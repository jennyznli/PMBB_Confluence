#!/bin/bash

#BSUB -J pmbb_commonsnps_qc
#BSUB -o pmbb_commonsnps_qc.%J.out
#BSUB -e pmbb_commonsnps_qc.%J.err
#BSUB -n 8
#BSUB -M 20GB
#BSUB -R "rusage[mem=20GB]"
#BSUB -W 02:00

module load plink/1.9-20210416
set -e 

SNP_FILE="/static/PMBB/PMBB-Release-2024-3.0/Imputed/common_snps_LD_pruned/PMBB-Release-2024-3.0_genetic_imputed.commonsnps"
HOME_DIR="/home/jennyzli/breast"
cd "$HOME_DIR"

### INDIVIDUAL QC ###

cp $SNP_FILE common_snps

# sex check
plink --bfile common_snps --check-sex --out common_snps
grep PROBLEM common_snps.sexcheck > common_snps.sexprobs
awk '{print $1"\t"$2}' common_snps.sexprobs > fail-sex.txt

# check ppl missingness
plink --bfile common_snps --missing --out common_snps

# check --het 
plink --bfile common_snps --het --out common_snps
# sort out ppl with heterozygosity  more/less than 3SD away from mean
bash het2.sh
# sort out ppl with genotype failure rate more than 0.03
bash missing-ind.sh

# LD pruning 
plink --bfile common_snps --indep-pairwise 250kb 50 0.2 --out common_snps

# pairwise IBS with --extract 
nohup plink --bfile common_snps --extract common_snps.prune.in --genome --out common_snps > ibs.log 2>&1 & 
tial -f nohup.out

cat fail-* | sort -k1 | uniq > fail-qc-inds.txt
plink --bfile common_snps --remove fail-qc-inds.txt --keep common_snps.prune.in --make-bed --out common_snps_ind

## SNP LEVEL ###

# cehck snps missingness 
0.05 
plink --bfile common_snps_ind --missing --out comon_snps_ind

# skipped diferential call rates for now 
plink --bfile common_snps_ind \
  --maf 0.05 --geno 0.05 --mind 0.05 \
  --make-bed --out common_snps_qc 

echo "Done with common SNP QC."

