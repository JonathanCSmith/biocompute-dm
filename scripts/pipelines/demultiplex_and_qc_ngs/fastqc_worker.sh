#!/usr/bin/env bash
#$ -S /bin/bash
#$ -V
#$ -pe smp 2

# Debug
echo "Array ID: ${SGE_TASK_ID}"
echo "File List: ${FILE_LIST}"

# Parse the file path
OLFIFS="${IFS}"
IFS=","
read -a FILE_PATHS <<< "${FILE_LIST}"
IFS="${OLFIFS}"
echo "File Paths: ${FILE_PATHS}"

INDEX=$((${SGE_TASK_ID}-1))
echo "File worker index: ${INDEX}"

FILE_PATH=${FILE_PATHS[${INDEX}]}
echo "Working on file: ${FILE_PATH}"

# Run FastQC on our specific file
module load bioinformatics/FastQC/0.11.3
fastqc "${FILE_PATH}"

