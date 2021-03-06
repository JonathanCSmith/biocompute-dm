#!/usr/bin/env bash
#$ -S /bin/bash
#$ -V

echo "Beginning collation process"

OLDIFS="${IFS}"

FILE_COUNT=$(wc -l < "${SAMPLE_CSV}")

if [ ${FILE_COUNT} -le 0 ]; then

    # There are no files identified
    >&2 echo "The provided sample csv did not contain any files."

# Ping back our info to the webserver
ssh ${USERNAME}@${HPC_IP} << EOF
    curl --form event="module_error" ${SERVER}\'/message/pipelines|${TICKET}\'
EOF

    exit

fi

IFS=','
REGEX="_alignment_metrics.txt"
REGEX_2=".*\\*.*"
while read SAMPLE_NAME SAMPLE_INPUT_PATH SAMPLE_OUTPUT_PATH EXTRA
do
    if [[ "$SAMPLE_INPUT_PATH}" =~ $REGEX ]]; then

        SAMPLE_PATH=$(dirname "${SAMPLE_INPUT_PATH}")
        ALIGNMENT="${SAMPLE_INPUT_PATH}"
        DUPLICATE=""=

        # We are looking for a specific file type
        for f in "${SAMPLE_PATH}"/*_duplicate_metrics.txt; do
            if [ "${DUPLICATE}" ]; then
                echo "More that one duplicate metrics was identified. ChIPSeq post align qc cannot determine which you wish to use. New file is: ${f}. This is a programming error and indicative of a current flaw in Biocompute that will be addressed asap"

# Ping back our info to the webserver
ssh ${USERNAME}@${HPC_IP} << EOF
curl --form event="module_error" ${SERVER}\'/message/pipelines|${TICKET}\'
EOF

                exit
            else
                DUPLICATE="${f}"
                echo "Identified duplicate file: ${DUPLICATE}"
            fi
        done

        # Validate not null
        if [[ -z "${ALIGNMENT}" || -z "${DUPLICATE}" ]]; then
            echo "One of the expected metrics for ${SAMPLE_NAME} was missing, this sample will be skipped"
            continue

        elif [[ "${ALIGNMENT}" =~ $REGEX_2 ]]; then
            echo "One of the expected metrics for ${SAMPLE_NAME} was missing, this sample will be skipped"
            continue

        elif [[ "${DUPLICATE}" =~ $REGEX_2 ]]; then
            echo "One of the expected metrics for ${SAMPLE_NAME} was missing, this sample will be skipped"
            continue

        fi

        # Move to current working directory
        echo "Moving metrics for: ${SAMPLE_NAME} to ./sample_${SAMPLE_NAME}/"
        mkdir "./sample_${SAMPLE_NAME}"
        cp "${ALIGNMENT}" "./sample_${SAMPLE_NAME}/sample_${SAMPLE_NAME}_alignment_metrics.txt"
        cp "${DUPLICATE}" "./sample_${SAMPLE_NAME}/sample_${SAMPLE_NAME}_duplicate_metrics.txt"

    fi
done < "${SAMPLE_CSV}"

# Validate cwd is not empty
FILE_COUNT=$(find ./ -maxdepth 1 -type d -name 'sample_*' | wc -l)
echo "Identified ${FILE_COUNT} viable metric sets to process"

# Execute perl
echo "Executing perl script"
"${PIPELINE_SOURCE}"/collect_ChIP_metrics.pl

# Cleanup cwd
for d in ./sample_*/; do
    echo "Removing directory: ${d}"
    rm -rf "${d}"
done

# Store this in pipeline outputs for later usage
mv "./merged_qc_table" "./pipeline_output/merged_qc_table.txt"

echo "Collation process finished"