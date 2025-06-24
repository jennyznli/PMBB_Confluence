#!/bin/bash
#BSUB -q MY_QUEUE
#BSUB -J "jobname[1-23]"
#BSUB -o log/jobname.%J_%I.out
#BSUB -e log/jobname.%J_%I.err
#BSUB -n 1
#BSUB -M 32000
#BSUB -R "rusage[mem=16GB]"
#BSUB -W 8:00

# optionally set working directory and data directory
WDIR="/path/to/working/directory"
cd $WDIR
# create log directory if it doesn't exist
mkdir -p log
# Get the chromosome number from the job index
CHROMOSOMES=($(seq 1 22) "X") # set up bash array for chromosomes
# select chromosome, adjust for 0-based indexing in bash
CHROMOSOME=${CHROMOSOMES[$((LSB_JOBINDEX-1))]}
# get file with chromosome
input_file="/path/to/input_file/chr${CHROMOSOME}_input_file.txt"


# STEP1: DO THE THING
# ------------------------
echo "STEP1: DO THE THING"
# your code here

# SYNTAX TO SUBMIT JOB
# ---------------------------
# submit chromosome 22
# bsub -hl -J "jobname[22]" < /path/to/this/script.bsub

# submit chromosomes 5-7
# bsub -hl -J "jobname[5-7]" < /path/to/this/script.bsub

# submit all chunks
# bsub -hl -J "jobname[1-23]" < /path/to/this/script.bsub
# bsub -hl < /path/to/this/script.bsub # use default job array in this script