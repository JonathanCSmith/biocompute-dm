#!/usr/bin/env bash
# =====================================================================================================================
#
# Submits modules for running in pipelines that have already been started. Should be submitted from the HPC
#
# Many of the variables will need to be added manually
#
# =====================================================================================================================

TICKET="NOT_VALID-NO_REPORTING_ALLOWED" # UUID of the module ticket - allows reporting

MODULE="" # Module name
MODULE_OUTPUT_DIRECTORY="" # Module root
WORKING_DIRECTORY="" # Pipeline root
PIPELINE_SOURCE_DIRECTORY="" # Pipeline sh files
INPUTS_STRING="" # Location of the csv file
SCRIPT_STRING="" # Location of the script to execute
USERNAME="" # Username on the hpc
HPC_IP="" # IP of the hpc as accessible from compute nodes
SERVER="" # IP of the webserver as accessible from head nodes
CLEANUP_SCRIPT="" # Location of the cleanup script

VARS="USERNAME=${USERNAME},HPC_IP=${HPC_IP},TICKET=${TICKET},PIPELINE_SOURCE=${PIPELINE_SOURCE_DIRECTORY},samples=${INPUTS_STRING},MODULE_OUTPUT_DIRECTORY=${MODULE_OUTPUT_DIRECTORY}"

# =============================================== ADDITIONAL VARS =====================================================
ADDITIONAL_VARS=""
VARS+="${ADDITIONAL_VARS}"
# =============================================== ADDITIONAL VARS =====================================================

# Submit the job and its monitor
ssh ${USERNAME}@${HPC_IP} << EOF
    source /etc/profile;
    JOBID=\$(qsub -V -N job-${TICKET} -o ${MODULE_OUTPUT_DIRECTORY}//${MODULE}_q_submission_output.log -e ${MODULE_OUTPUT_DIRECTORY}//${MODULE}_q_submission_error.log -wd ${WORKING_DIRECTORY} -v ${VARS} ${SCRIPT_STRING} | cut -d ' ' -f 3);
    echo Job Id: \$JOBID
    qsub -V -hold_jid \$JOBID -N cleanup-${TICKET} -o ${MODULE_OUTPUT_DIRECTORY}//${MODULE}_module_cleanup.log -e ${MODULE_OUTPUT_DIRECTORY}//${MODULE}_module_cleanup.log -v USERNAME=${USERNAME},HPC_IP=${HPC_IP},SERVER=${SERVER},TICKET=${TICKET},JOBID=\$JOBID ${CLEANUP_SCRIPT}
EOF