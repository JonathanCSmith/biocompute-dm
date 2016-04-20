#!/usr/bin/env bash
#$ -j y
#$ -S /bin/bash
#S -V

OLDIFS="${IFS}"

# ================================================= BUILD OUR IO VALUES ===============================================
echo "Beginning load module"

# Note we are only interested in the first line as we are only expecting 1 folder!
IFS="," read NAME DATA_INPUT_DIRECTORY DATA_OUTPUT_DIRECTORY EXTRA < <(sed -n 1p < "${SAMPLE_CSV}")

echo "Submission input directory: ${DATA_INPUT_DIRECTORY}"
echo "Submission output directory: ${DATA_OUTPUT_DIRECTORY}"
echo "Extra: ${EXTRA}"
echo "Module output directory: ${MODULE_OUTPUT_DIRECTORY}"
# ========================================== FINISHED BUILDING OUR IO VALUES ==========================================

# ============================================ IDENTIFY AND MOVE OUR FILES ============================================
echo "================== Processing read 1s ==================="

for f in ${DATA_INPUT_DIRECTORY}/*_R1.fastq # We are only interested in demuxed files!
do
    echo "Processing File: ${f}"

    NAME="${f##*/}"
    NAME_WITHOUT_EXTENSION="${NAME%_R1.*}"
    SAMPLE_DIR="${DATA_OUTPUT_DIRECTORY}/${NAME_WITHOUT_EXTENSION}"

    # Make the directory if it does not exist yet
    echo "Constructing directory: ${SAMPLE_DIR} for ${NAME}"
    mkdir "${SAMPLE_DIR}"

    # Move the file to the directory
    mv "${f}" "${SAMPLE_DIR}/${NAME}"
done
echo "================== Processing read 2s ==================="

for f in ${DATA_INPUT_DIRECTORY}/*_R2.fastq
do
    echo "Processing File: ${f}"

    NAME="${f##*/}"
    NAME_WITHOUT_EXTENSION="${NAME%_R2.*}"
    SAMPLE_DIR="${DATA_OUTPUT_DIRECTORY}/${NAME_WITHOUT_EXTENSION}"

    # Make the directory if it does not exist yet
    echo "Constructing directory: ${SAMPLE_DIR} for ${NAME}"
    mkdir "${SAMPLE_DIR}"

    # Move the file to the directory
    mv "${f}" "${SAMPLE_DIR}/${NAME}"
done
# ============================================ IDENTIFY AND MOVE OUR FILES ============================================

IFS="${OLDIFS}"