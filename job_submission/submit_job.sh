#!/usr/bin/env bash
echo "Beginning job submission"

# Take all script input arguments and build strings for them
for i in "$@"
do
case ${i} in
    # Ticket for transactions - also used as the job name for the outermost pipeline specific job
    -t=*|--ticket=*)
    TICKET="${i#*=}"
    shift
    ;;

    # Script execution path - note this must be relative to the OGS
    -s=*|--script=*)
    SCRIPT_STRING="${i#*=}"
    shift
    ;;

    # Local Output Directory
    -l=*)
    LOCAL_OUTPUT_DIRECTORY="${i#*=}"
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

    # Unknown
    *)
    ;;
esac
done

# Combine the strings in a meaningful manner
VARS="ticket=${TICKET},working_directory=${WORKING_DIRECTORY},samples=${INPUTS_STRING}"
if [ "${VARIABLES_STRING}" ]; then
    VARS+=",${VARIABLES_STRING}"
fi

# Copy over script file
echo "Current: ${PWD}"
echo "Pipeline: ${WORKING_DIRECTORY}"
scp ./cleanup.sh biocis@10.202.64.28:~

# Create command string
STRING="JOBID=\$\(qsub -cwd -N ${TICKET} -v ${VARS} ${SCRIPT_STRING}\)"

# Submit the job and its monitor
OUTPUT_FILE=${LOCAL_OUTPUT_DIRECTORY}
OUTPUT_FILE+="/header_node_output.txt"
ssh biocis@10.202.64.28 "
    eval "${STRING}"
    eval "JOBID=\$\( echo \$\{JOBID\} \| grep -o -E '[0-9]+' \)"
    qsub -hold_jid "${JOBID}" -N CLEANUP ./cleanup.sh -v "TICKET=${TICKET},JOBID=${JOBID}" # The monitor
" > ${OUTPUT_FILE} 2>&1

# Move the output into our working directory
#mv ./header_node_output.txt ${WORKING_DIRECTORY}
#rm ./header_node_output.txt

echo "Job submission complete"
