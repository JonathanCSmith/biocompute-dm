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

# ================================================== SUBMISSION DATA ==================================================
# Traverse the output path to generate a csv series of file paths
FILE_LIST=""
FILE_COUNT=0
for f in ${DATA_OUTPUT_DIRECTORY}/*.fastq.gz # We are only interested in demuxed files!
do
    FILE_LIST+="${f},"
    FILE_COUNT=$((FILE_COUNT+1))
done

echo "File Count ${FILE_COUNT}"
echo "File List ${FILE_LIST}"

# Error catching for no files (logic is inverted here because my IDE was playing havoc with the heredoc below)
if [ ${FILE_COUNT} -ne 0 ]; then
# Get rid of the extra comma
FILE_LIST="${FILE_LIST%?}"

# Pass to an array job to handle
ssh ${USERNAME}@${HPC_IP} << END
    source /etc/profile;
    qsub -V -t 1-${FILE_COUNT}:1 -v "FILE_LIST=${FILE_LIST},LOGGING_DIR=${MODULE_OUTPUT_DIRECTORY}" "${PIPELINE_SOURCE}//fastqc_worker.sh"
END

# Safe exit so we don't trip our error code below
exit
fi
# ================================================== SUBMISSION DATA ==================================================

# =================================================== ERROR HANDLING ==================================================
>&2 echo "The provided directory did not contain any fastq.gz files - this indicates that the demux process did not run properly"
# Ping back our info to the webserver - TODO: Silence it?
ssh ${USERNAME}@${HPC_IP} << EOF
    curl --form event="module_error" ${SERVER}\'/message/pipelines|${TICKET}\'
EOF
IFS="${OLDIFS}"
# =================================================== ERROR HANDLING ==================================================