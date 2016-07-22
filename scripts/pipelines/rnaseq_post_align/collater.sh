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

else

    # Loop through samples
    echo "Beginning loop for samples"
    IFS=","
    while read SAMPLE_NAME SAMPLE_INPUT_PATH SAMPLE_OUTPUT_PATH EXTRA
    do

        ALIGNMENT=""
        DUPLICATE=""
        METRICS=""
        METRICS0_4=""
        METRICS4_8=""
        METRICS0_8=""

        # We are looking for a specific file type
        for f in ${d}/*_alignment_metrics.txt; do
            if [ "${ALIGNMENT}" ]; then
                echo "More that one alignment metrics was identified. RNAseq post align qc cannot determine which you wish to use. New file is: ${f}."

# Ping back our info to the webserver
ssh ${USERNAME}@${HPC_IP} << EOF
curl --form event="module_error" ${SERVER}\'/message/pipelines|${TICKET}\'
EOF

                exit
            else
                ALIGNMENT="${f}"
                echo "Identified alignment file: ${ALIGNMENT}"
            fi
        done

        # We are looking for a specific file type
        for f in ${d}/*_duplicate_metrics.txt; do
            if [ "${DUPLICATE}" ]; then
                echo "More that one duplicate metrics was identified. RNAseq post align qc cannot determine which you wish to use. New file is: ${f}."

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

        # We are looking for a specific file type
        for f in ${d}/*_RnaSeqMetrics.txt; do
            if [ "${METRICS}" ]; then
                echo "More that one rnaseq metrics was identified. RNAseq post align qc cannot determine which you wish to use. New file is: ${f}."

# Ping back our info to the webserver
ssh ${USERNAME}@${HPC_IP} << EOF
curl --form event="module_error" ${SERVER}\'/message/pipelines|${TICKET}\'
EOF

                exit
            else
                METRICS="${f}"
                echo "Identified metrics file: ${METRICS}"
            fi
        done

        # We are looking for a specific file type
        for f in ${d}/*_RnaSeqMetrics_upto4kb.txt; do
            if [ "${METRICS0_4}" ]; then
                echo "More that one duplicate metrics was identified. RNAseq post align qc cannot determine which you wish to use. New file is: ${f}."

# Ping back our info to the webserver
ssh ${USERNAME}@${HPC_IP} << EOF
curl --form event="module_error" ${SERVER}\'/message/pipelines|${TICKET}\'
EOF

                exit
            else
                METRICS0_4="${f}"
                echo "Identified metrics (up to 4kb) file: ${METRICS0_4}"
            fi
        done

        # We are looking for a specific file type
        for f in ${d}/*_RnaSeqMetrics_4to8kb.txt; do
            if [ "${METRICS4_8}" ]; then
                echo "More that one duplicate metrics was identified. RNAseq post align qc cannot determine which you wish to use. New file is: ${f}."

# Ping back our info to the webserver
ssh ${USERNAME}@${HPC_IP} << EOF
curl --form event="module_error" ${SERVER}\'/message/pipelines|${TICKET}\'
EOF

                exit
            else
                METRICS4_8="${f}"
                echo "Identified metrics (up to 4kb) file: ${METRICS4_8}"
            fi
        done

        # We are looking for a specific file type
        for f in ${d}/*_RnaSeqMetrics_8kb.txt; do
            if [ "${METRICS0_8}" ]; then
                echo "More that one duplicate metrics was identified. RNAseq post align qc cannot determine which you wish to use. New file is: ${f}."

# Ping back our info to the webserver
ssh ${USERNAME}@${HPC_IP} << EOF
curl --form event="module_error" ${SERVER}\'/message/pipelines|${TICKET}\'
EOF

                exit
            else
                METRICS0_8="${f}"
                echo "Identified metrics (8kb) file: ${METRICS0_8}"
            fi
        done

        # Validate not null
        if [[ -z "${ALIGNMENT}" || -z "${DUPLICATE}" || -z "${METRICS}" || -z "${METRICS0_4}" || -z "${METRICS4_8}"|| -z "${METRICS0_8}" ]]; then
            echo "One of the expected metrics for ${SAMPLE_NAME} was missing, this sample will be skipped"
            continue

        elif [[ "${ALIGNMENT}" =~ ".**.*" ]]; then
            echo "One of the expected metrics for ${SAMPLE_NAME} was missing, this sample will be skipped"
            continue

        elif [[ "${DUPLICATE}" =~ ".**.*" ]]; then
            echo "One of the expected metrics for ${SAMPLE_NAME} was missing, this sample will be skipped"
            continue

        elif [[ "${METRICS}" =~ ".**.*" ]]; then
            echo "One of the expected metrics for ${SAMPLE_NAME} was missing, this sample will be skipped"
            continue

        elif [[ "${METRICS0_4}" =~ ".**.*" ]]; then
            echo "One of the expected metrics for ${SAMPLE_NAME} was missing, this sample will be skipped"
            continue

        elif [[ "${METRICS4_8}" =~ ".**.*" ]]; then
            echo "One of the expected metrics for ${SAMPLE_NAME} was missing, this sample will be skipped"
            continue

        elif [[ "${METRICS0_8}" =~ ".**.*" ]]; then
            echo "One of the expected metrics for ${SAMPLE_NAME} was missing, this sample will be skipped"
            continue

        fi

        # Move to current working directory
        echo "Moving metrics for: ${SAMPLE_NAME} to ./sample_${SAMPLE_NAME}/"
        mkdir "./sample_${SAMPLE_NAME}"
        cp "${ALIGNMENT}" "./sample_${SAMPLE_NAME}/sample_${SAMPLE_NAME}_alignment_metrics.txt"
        cp "${DUPLICATE}" "./sample_${SAMPLE_NAME}/sample_${SAMPLE_NAME}_duplicate_metrics.txt"
        cp "${DUPLICATE}" "./sample_${SAMPLE_NAME}/sample_${SAMPLE_NAME}_RnaSeqMetrics.txt"
        cp "${DUPLICATE}" "./sample_${SAMPLE_NAME}/sample_${SAMPLE_NAME}_RnaSeqMetrics_upto4Kb.txt"
        cp "${DUPLICATE}" "./sample_${SAMPLE_NAME}/sample_${SAMPLE_NAME}_RnaSeqMetrics_4To8Kb.txt"
        cp "${DUPLICATE}" "./sample_${SAMPLE_NAME}/sample_${SAMPLE_NAME}_RnaSeqMetrics_8Kb.txt"

    done < "${SAMPLE_CSV}"

    # Validate cwd is not empty
    FILE_COUNT=$(find ./ -maxdepth 1 -type d -name 'sample_*' | wc -l)
    echo "Identified ${FILE_COUNT} viable metric sets to process"

    # Execute perl
    echo "Executing perl script"
    "${PIPELINE_SOURCE}"/collect_RNASeq_metrics.pl

    # Cleanup cwd
    for d in ./sample_*/; do
        echo "Removing directory: ${d}"
        rm -rf "${d}"
    done
fi

echo "Collation process finished"