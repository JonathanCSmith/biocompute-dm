#!/usr/bin/env bash
echo "Beginning job submission"

# Take all script input arguments and build strings for them
for i in "$@"
do
case ${i} in
    # HPC Username
    -u=*)
    USERNAME="${i#*=}"
    shift
    ;;

    -h=*)
    HPC_IP="${i#*=}"
    shift
    ;;

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

    # Cleanup script path
    -c=*)
    CLEANUP_SCRIPT="${i#*=}"
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

# Good Logging
echo "Current: ${PWD}"
echo "Pipeline: ${WORKING_DIRECTORY}"
echo Module to Submit: ${MODULE}
echo Ticket Id: ${TICKET}
echo Working Directory: ${WORKING_DIRECTORY}
echo Variables: ${VARS}
echo Script: ${SCRIPT_STRING}
echo Beginning SSH

# Submit the job and its monitor
ssh ${USERNAME}@${HPC_IP} << EOF
    echo Beginning submission log for module: ${MODULE}
    JOBID=\$(qsub -N job-${TICKET} -o ${MODULE}_module_output.log -e ${MODULE}_module_error.log -wd ${WORKING_DIRECTORY} -v ${VARS} ${SCRIPT_STRING} | cut -d ' ' -f 3);
    echo Job Id: \$JOBID
    qsub -hold_jid \$JOBID -N cleanup-${TICKET} -o ${MODULE}_module_cleanup.log -e ${MODULE}_module_cleanup.log -v USERNAME=${USERNAME},HPC_IP=${HPC_IP},SERVER=${SERVER},TICKET=${TICKET},JOBID=\$JOBID ${CLEANUP_SCRIPT}
EOF

echo "Job submission complete"
