#!/usr/bin/env bash
#$ -S /bin/bash
#$ -V

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

    # Pass to an array job to handle
    echo "Confirmed file count: ${FILE_COUNT}"
JOBID=$(ssh ${USERNAME}@${HPC_IP} << END
    source /etc/profile;
    JOBID=\$(qsub -t 1-${FILE_COUNT} -v "DATA_FILE="${SAMPLE_CSV}",ref="${ref}",ann="${ann}",PIPELINE_SOURCE="${PIPELINE_SOURCE} -o "${MODULE_OUTPUT_DIRECTORY}" -e "${MODULE_OUTPUT_DIRECTORY}" "${PIPELINE_SOURCE}//align_worker.sh" | cut -d ' ' -f 3);
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
            if [[ ${RESULT} =~ ${REGEX} ]]; then
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
        mv "${MODULE_OUTPUT_DIRECTORY}//align_worker.sh.o${JOBID}.${i}" "${MODULE_OUTPUT_DIRECTORY}//align_worker_arrayid_${i}_out.txt"
        mv "${MODULE_OUTPUT_DIRECTORY}//align_worker.sh.po${JOBID}.${i}" "${MODULE_OUTPUT_DIRECTORY}//align_worker_arrayid_${i}_error.txt"
    done
fi

IFS="${OLDIFS}"