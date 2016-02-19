#!/usr/bin/env bash
echo "Beginning job submission"

# Take all script input arguments and build strings for them
for i in "$@"
do
case ${i} in
    # Module name for logging etc
    -m=*)
    MODULE="${i#*=}"
    shift
    ;;

    # Ticket for transactions - also used as the job name for the outermost pipeline specific job
    -t=*|--ticket=*)
    TICKET="${i#*=}"
    shift
    ;;

    # Script execution path - note this must be relative to the OGS
    -e=*)
    SCRIPT_STRING="${i#*=}"
    shift
    ;;

    # Local Output Directory
    -l=*)
    LOCAL_OUTPUT_DIRECTORY="${i#*=}"
    shift
    ;;

    # Pipeline files directory
    -p=*)
    PIPELINE_SOURCE_DIRECTORY="${i#*=}"
    shift
    ;;

    # Working directory argument
    -w=*|--working=*)
    WORKING_DIRECTORY="${i#*=}"
    shift
    ;;

    # Specific argument for the inputs file
    -i=*|--inputs=*)
    INPUTS_STRING="${i#*=}"
    shift
    ;;

    # Generic string containing all csvs that should be passed to qsub - note this MUST conform to qsub expectations
    -v=*|--variables)
    VARIABLES_STRING="${i#*=}"
    shift
    ;;

    # Ping back location
    -s=*)
    SERVER="${i#*=}"
    shift
    ;;

    # Unknown
    *)
    ;;
esac
done

# Combine the strings in a meaningful manner
VARS="ticket=${TICKET},pipeline_source=${PIPELINE_SOURCE_DIRECTORY},samples=${INPUTS_STRING}"
if [ "${VARIABLES_STRING}" ]; then
    VARS="${VARS},${VARIABLES_STRING}"
fi

# Copy over script file
echo "Current: ${PWD}"
echo "Pipeline: ${WORKING_DIRECTORY}"
scp ./cleanup.sh biocis@10.202.64.28:~

echo Module to Submit: ${MODULE}
echo Ticket Id: ${TICKET}
echo Working Directory: ${WORKING_DIRECTORY}
echo Variables: ${VARS}
echo Script: ${SCRIPT_STRING}

# Submit the job and its monitor
OUTPUT_FILE=${LOCAL_OUTPUT_DIRECTORY}
OUTPUT_FILE+="/module_submission.txt"
ssh biocis@10.202.64.28 << EOF > ${OUTPUT_FILE} 2>&1
    echo Beginning submission log for module: ${MODULE}
    JOBID=\$(qsub -o ${MODULE}_output.log -e ${MODULE}_error.log -N job-${TICKET} -wd ${WORKING_DIRECTORY} -v ${VARS} ${SCRIPT_STRING} | cut -d ' ' -f 3);
    echo Job Id: \$JOBID
    qsub -hold_jid \$JOBID -N cleanup-"${TICKET}" -v SERVER="${SERVER}",TICKET="${TICKET}",JOBID=\$JOBID ./cleanup.sh
EOF

echo "Job submission complete"
