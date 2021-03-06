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
    PIPELINE_SOURCE="${i#*=}"
    shift
    ;;

    # Working directory argument
    -w=*|--working=*)
    WORKING_DIRECTORY="${i#*=}"
    shift
    ;;

    # Specific argument for the inputs file
    -i=*|--inputs=*)
    SAMPLE_CSV="${i#*=}"
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
VARS="USERNAME=${USERNAME},HPC_IP=${HPC_IP},TICKET=${TICKET},PIPELINE_SOURCE=${PIPELINE_SOURCE},SAMPLE_CSV=${SAMPLE_CSV},PIPELINE_OUTPUT_DIRECTORY=${PIPELINE_OUTPUT_DIRECTORY},MODULE_OUTPUT_DIRECTORY=${MODULE_OUTPUT_DIRECTORY}"
if [ "${VARIABLES_STRING}" ]; then
    VARS="${VARS},${VARIABLES_STRING}"
fi

# Good logging
echo "Current: ${PWD}"
echo "Pipeline: ${WORKING_DIRECTORY}"
echo "Module to Submit: ${MODULE}"
echo "Ticket Id: ${TICKET}"
echo "Working Directory: ${WORKING_DIRECTORY}"
echo "Variables: ${VARS}"
echo "Script: ${SCRIPT_STRING}"
echo "Beginning SSH"

# Pingback for status - TODO: Silence it?
source ${SCRIPT_STRING}  # Run locally - use the current env variables
curl --form event="module_end" --form sub="${SUB_TIME}" --form stat="${START_TIME}" --form end="${END_TIME}" 127.0.0.1:${SERVER}/"message/pipelines|${TICKET}"
echo "Job submission complete"
