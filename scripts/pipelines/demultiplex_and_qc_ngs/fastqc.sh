#!/usr/bin/env bash
#$ -S /bin/bash
#$ -V

OLDIFS="${IFS}"

# ================================================= BUILD OUR IO VALUES ===============================================
echo "Beginning fastqc module"

# Note we are only interested in the first line as we are only expecting 1 folder!
IFS="," read NAME DATA_INPUT_DIRECTORY DATA_OUTPUT_DIRECTORY EXTRA < <(sed -n 1p < "${SAMPLE_CSV}")

echo "Submission input directory: ${DATA_INPUT_DIRECTORY}"
echo "Submission output directory: ${DATA_OUTPUT_DIRECTORY}"
echo "Extra: ${EXTRA}"
echo "Module output directory: ${MODULE_OUTPUT_DIRECTORY}"
# ========================================== FINISHED BUILDING OUR IO VALUES ==========================================

# ============================================== SUBMISSION PROPERTIES=================================================
# Traverse the output path to generate a csv series of file paths
FILE_LIST="${MODULE_OUTPUT_DIRECTORY}/samples.csv"
FILE_COUNT=0
for f in ${DATA_OUTPUT_DIRECTORY}/*.fastq.gz # We are only interested in demuxed files!
do
    FILE_COUNT=$((FILE_COUNT+1))
    echo "${f}"
done > "${FILE_LIST}"

echo "File Count ${FILE_COUNT}"
# ============================================== SUBMISSION PROPERTIES=================================================

# ================================================== SUBMISSION DATA ==================================================
if [ ${FILE_COUNT} -le 0 ]; then

    # There are no files identified
    >&2 echo "The provided directory did not contain any fastq.gz files - this indicates that the demux process did not run properly"

# Ping back our info to the webserver
ssh ${USERNAME}@${HPC_IP} << EOF
    curl --form event="module_error" ${SERVER}\'/message/pipelines|${TICKET}\'
EOF

else

    # Only 1 file present - array job is not suitable
    if [ ${FILE_COUNT} -eq 1 ]; then

# Don't call as an array job
JOBID=$(ssh ${USERNAME}@${HPC_IP} << END
    source /etc/profile;
    JOBID=\$(qsub -v "FILE_LIST=${FILE_LIST},SGE_TASK_ID=1" -o "${MODULE_OUTPUT_DIRECTORY}//fastqc_worker_out.txt" -e "${MODULE_OUTPUT_DIRECTORY}//fastqc_worker_error.txt" "${PIPELINE_SOURCE}//fastqc_worker.sh" | cut -d ' ' -f 3);
    echo \$JOBID
END
)

        echo "Retrieved Job Id: ${JOBID}"

        HAS_RUNNING=1
        while [ ${HAS_RUNNING} -eq 1 ]
        do
            # Don't poll too often
            sleep 100

# Poll the job id
RESULT=$(ssh ${USERNAME}@${HPC_IP} << END
    source /etc/profile;
    RESULT=\$(qacct -j ${JOBID} 2>&1);
    echo \$RESULT
END
)

            # If our job was not present in qacct
            REGEX="error*"
            if [[ "${RESULT}" =~ ${REGEX} ]]; then
                continue
            else
                HAS_RUNNING=0
            fi
        done

        echo "Successfully waited for the job to finish."

    # More than one file - use an array job
    else

# Pass to an array job to handle
echo "Confirmed file count: ${FILE_COUNT}"
JOBID=$(ssh ${USERNAME}@${HPC_IP} << END
    source /etc/profile;
    JOBID=\$(qsub -t 1-${FILE_COUNT} -v "FILE_LIST=${FILE_LIST}" -o "${MODULE_OUTPUT_DIRECTORY}" -e "${MODULE_OUTPUT_DIRECTORY}" "${PIPELINE_SOURCE}//fastqc_worker.sh" | cut -d ' ' -f 3);
    echo \$JOBID
END
)

        echo "Retrieved Job Id: ${JOBID}"
        IFS='.' read -r JOBID leftover <<< "${JOBID}"
        echo "Parsed Job Id: ${JOBID}"

        # Poll job ids and wait for completion
        HAS_RUNNING=1
        while [ ${HAS_RUNNING} -eq 1 ]
        do

            # Don't poll too often
            sleep 100

            # Poll each job id
            NOT_PRESENT=0
            for i in $(seq 1 ${FILE_COUNT})
            do

                echo "Searching for job: ${i}."

RESULT=$(ssh ${USERNAME}@${HPC_IP} << END
    source /etc/profile;
    RESULT=\$(qacct -j ${JOBID} -t ${i} 2>&1);
    echo \$RESULT
END
)

                # If our job was not present in qacct
                echo "Job Query returned: ${RESULT}"
                REGEX="error*"
                if [[ "${RESULT}" =~ ${REGEX} ]]; then
                    NOT_PRESENT=1
                    break;
                fi
            done

            # Evaluate the outcome of the for loop
            if [ ${NOT_PRESENT} -eq 0 ]; then
                HAS_RUNNING=0
            fi
        done

        echo "Successfully waited for all array jobs to finish, beginning log move."

        # Rename logs and exit
        for i in $(seq 1 ${FILE_COUNT})
        do
            mv "${MODULE_OUTPUT_DIRECTORY}//fastqc_worker.sh.o${JOBID}.${i}" "${MODULE_OUTPUT_DIRECTORY}//fastqc_worker_arrayid_${i}_out.txt"
            mv "${MODULE_OUTPUT_DIRECTORY}//fastqc_worker.sh.po${JOBID}.${i}" "${MODULE_OUTPUT_DIRECTORY}//fastqc_worker_arrayid_${i}_error.txt"
        done
    fi
fi
# ================================================== SUBMISSION DATA ==================================================
IFS="${OLDIFS}"
