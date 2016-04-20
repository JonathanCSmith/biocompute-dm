#!/usr/bin/env bash
#$ -S /bin/bash
#$ -V

echo "Beginning colation process"

OLDIFS="${IFS}"

FILE_COUNT=$(wc -l < "${SAMPLE_CSV}")

if [ ${FILE_COUNT} -le 0 ]; then

    # There are no files identified
    >&2 echo "The provided sample csv did not contain any files."

# Ping back our info to the webserver
ssh ${USERNAME}@${HPC_IP} << EOF
    curl --form event="module_error" ${SERVER}\'/message/pipelines|${TICKET}\'
EOF

else

    # Loop through samples
    echo "Beginning loop for samples"
    IFS=","
    while read SAMPLE_NAME SAMPLE_INPUT_PATH SAMPLE_OUTPUT_PATH EXTRA
    do

        ALIGNMENT=""
        DUPLICATE=""
        COVERAGE=""

        # Expectation of nested pipeline datasets within sample dir (current layout of biocompute - may change, see the two errors below)
        for d in "${SAMPLE_INPUT_PATH}/*/"; do

            # We are looking for a specific file type
            for f in d/*_alignment_metrics.fastq; do
                if [ "${ALIGNMENT}" ]; then
                    echo "More that one alignment metrics was identified. Exome post align qc cannot determine which you wish to use. New file is: ${ALIGNMENT}. This is a programming error and indicative of a current flaw in Biocompute that will be addressed asap"

# Ping back our info to the webserver
ssh ${USERNAME}@${HPC_IP} << EOF
    curl --form event="module_error" ${SERVER}\'/message/pipelines|${TICKET}\'
EOF

                    exit
                else
                    echo "Identified aligment file: ${ALIGNMENT}"
                    ALIGNMENT="${f}"
                fi
            done

            # We are looking for a specific file type
            for f in d/*_duplicate_metrics.fastq; do
                if [ "${DUPLICATE}" ]; then
                    echo "More that one duplicate metrics was identified. Exome post align qc cannot determine which you wish to use. New file is: ${DUPLICATE}. This is a programming error and indicative of a current flaw in Biocompute that will be addressed asap"

# Ping back our info to the webserver
ssh ${USERNAME}@${HPC_IP} << EOF
    curl --form event="module_error" ${SERVER}\'/message/pipelines|${TICKET}\'
EOF

                    exit
                else
                    echo "Identified aligment file: ${DUPLICATE}"
                    DUPLICATE="${f}"
                fi
            done

            # We are looking for a specific file type
            for f in d/*_coverage_metrics.fastq; do
                if [ "${COVERAGE}" ]; then
                    echo "More that one coverage metrics was identified. Exome post align qc cannot determine which you wish to use. New file is: ${COVERAGE}. This is a programming error and indicative of a current flaw in Biocompute that will be addressed asap"

# Ping back our info to the webserver
ssh ${USERNAME}@${HPC_IP} << EOF
    curl --form event="module_error" ${SERVER}\'/message/pipelines|${TICKET}\'
EOF

                    exit
                else
                    echo "Identified aligment file: ${COVERAGE}"
                    COVERAGE="${f}"
                fi
            done

            # Validate not null
            if [[ -z "${ALIGNMENT}" || -z "${DUPLICATE}" || -z "${COVERAGE}" ]]; then
                echo "One of the expected metrics for ${SAMPLE_NAME} was missing, this sample will be skipped"
                continue
            fi

            # Move to current working directory
            echo "Moving metrics for: ${SAMPLE_NAME} to ./sample_${SAMPLE_NAME}/"
            mkdir "./sample_${SAMPLE_NAME}"
            cp "${ALIGNMENT}" "./sample_${SAMPLE_NAME}/${SAMPLE_NAME}_alignment_metrics"
            cp "${DUPLICATE}" "./sample_${SAMPLE_NAME}/${SAMPLE_NAME}_duplicate_metrics"
            cp "${COVERAGE}" "./sample_${SAMPLE_NAME}/${SAMPLE_NAME}_coverage_metrics"

        done
    done < "${SAMPLE_CSV}"

    # Validate cwd is not empty
    FILE_COUNT=$(wc -l < "./sample_*")
    echo "Identified ${FILE_COUNT} viable metric sets to process"

    # Execute perl
    echo "Executing perl script"
    "${PIPELINE_SOURCE}"/get_qc_data.pl

    # Cleanup cwd
    for d in "./sample_*/"; do
        echo "Removing directory: ${d}"
        rm -rf "${d}"
    done
fi

echo "Colation process finished"