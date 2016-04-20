#!/usr/bin/env bash
#$ -j y
#$ -S /bin/bash
#S -V

OLDIFS="${IFS}"
echo "Beginning load module"
RX_1="*_R1.fastq"
RX_2="*_R2.fastq"

# ============================================ IDENTIFY AND MOVE OUR FILES ============================================
IFS=","
while read NAME DATA_INPUT_DIRECTORY DATA_OUTPUT_DIRECTORY EXTRA; do
    if [[ ${DATA_INPUT_DIRECTORY} =~ ${RX_1} ]]; then
        NAME="${DATA_INPUT_DIRECTORY##*/}"
        NAME_WITHOUT_EXTENSION="${NAME%_R1.*}"

    elif [[ ${DATA_INPUT_DIRECTORY} =~ ${RX_2} ]]; then
        NAME="${DATA_INPUT_DIRECTORY##*/}"
        NAME_WITHOUT_EXTENSION="${NAME%_R2.*}"

    else
        echo "Unknown file: ${DATA_INPUT_DIRECTORY}, this sample will not be processed."
        continue
    fi

    # Make the directory if it does not exist yet
    SAMPLE_DIR="${DATA_OUTPUT_DIRECTORY}/${NAME_WITHOUT_EXTENSION}"
    echo "Constructing directory: ${SAMPLE_DIR} for ${NAME}"
    mkdir "${SAMPLE_DIR}"

    # Move the file to the directory
    mv "${DATA_INPUT_DIRECTORY}" "${SAMPLE_DIR}/${NAME}"
done < "${SAMPLE_CSV}"
# ============================================ IDENTIFY AND MOVE OUR FILES ============================================

IFS="${OLDIFS}"