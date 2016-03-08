#!/usr/bin/env bash
#$ -S /bin/bash
# -V

OLDIFS="${IFS}"

# ================================================= BUILD OUR IO VALUES ===============================================
echo "Beginning demultiplexing module"

# Note we are only interested in the first line as we are only expecting 1 folder!
IFS="," read NAME DATA_INPUT_DIRECTORY DATA_OUTPUT_DIRECTORY EXTRA < <(sed -n 1p < "${SAMPLE_CSV}")

echo "Submission input directory: ${DATA_INPUT_DIRECTORY}"
echo "Submission output directory: ${DATA_OUTPUT_DIRECTORY}"
echo "Extra: ${EXTRA}"
echo "Module output directory: ${MODULE_OUTPUT_DIRECTORY}"
# ========================================== FINISHED BUILDING OUR IO VALUES ==========================================

# ==================================================== SAMPLE DATA ====================================================
# Traverse the output path to generate a csv series of file paths
FILE_LIST=""
FILE_COUNT=0
for f in ${DATA_OUTPUT_DIRECTORY}/*; # We are only interested in demuxed files!
do
    if [[ ${f} == *.fastq.gz ]]; then
        echo "Processing ${f}"
        NAME="${f##*/}"
        NAME_WITHOUT_EXTENSION="${NAME%%.*}"
        NAME_WITHOUT_EXTENSION="${NAME%%_*}"

        # Identify the correct directory name
        IFS="_" read -r ID leftover <<< "${NAME_WITHOUT_EXTENSION}"
        echo "Sample identifier parsed as: ${ID}"

        DIR="${DATA_OUTPUT_DIRECTORY}"
        if [[ ${ID} == Undetermine* ]]; then
            DIR+="/${NAME_WITHOUT_EXTENSION}"
        else
            DIR+="/${ID}"
        fi

        # Make the directory if it does not exist yet
        echo "Constructing directory: ${DIR}"
        mkdir "${DIR}"

        # Move the file to the directory
        mv "${DATA_OUTPUT_DIRECTORY}/${f}" "${DATA_OUTPUT_DIRECTORY}/${DIR}/${f}"

    elif [[ ${f} == *_fastqc.html || ${f} == *_fastqc.zip ]]; then
        echo "Processing ${f}"
        NAME="${f##*/}"
        NAME_WITHOUT_EXTENSION="${NAME%.*}"
        NAME_WITHOUT_EXTENSION="${NAME%_fastqc}"
        NAME_WITHOUT_EXTENSION="${NAME%%_*}"

        # Identify the correct directory name
        IFS="_" read -r ID leftover <<< "${NAME_WITHOUT_EXTENSION}"
        echo "Sample identifier parsed as: ${ID}"

        DIR="${DATA_OUTPUT_DIRECTORY}"
        if [[ ${ID} == Undetermine* ]]; then
            DIR+="/${NAME_WITHOUT_EXTENSION}"
        else
            DIR+="/${ID}"
        fi

        # Make the directory if it does not exist yet
        echo "Constructing directory: ${DIR}"
        mkdir "${DIR}"

        # Move the file to the directory
        mv "${DATA_OUTPUT_DIRECTORY}/${f}" "${DATA_OUTPUT_DIRECTORY}/${DIR}/${f}"
    fi
done
# ==================================================== SAMPLE DATA ====================================================

# ===================================================== REPORTING =====================================================


IFS="${OLDIFS}"
# ===================================================== REPORTING =====================================================