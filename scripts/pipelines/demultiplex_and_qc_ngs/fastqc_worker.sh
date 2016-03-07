#!/usr/bin/env bash
#$ -S /bin/bash
# -V
#$ -pe smp 2

# Parse the file path
OLFIFS="${IFS}"
IFS=","
read -a FILE_PATHS <<< "${FILE_LIST}"
IFS="${OLFIFS}"
FILE_PATH=${FILE_PATHS[${SGE_TASK_ID}]}

# Run FastQC on our specific file
module load bioinformatics/fastqc/0.11.3
fastqc ${FILE_PATH}
