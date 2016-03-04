#!/usr/bin/env bash
# =====================================================================================================================
#
# Submits modules for running in pipelines that have already been started. Should be submitted from the HPC
#
# Many of the variables will need to be added manually
#
# =====================================================================================================================

TICKET="NOT_VALID-NO_REPORTING_ALLOWED"

MODULE=""
OUTPUT_DIRECTORY=""
WORKING_DIRECTORY=""
PIPELINE_SOURCE_DIRECTORY=""
INPUTS_STRING=""
SCRIPT_STRING=""
USERNAME=""
HPC_IP=""
SERVER=""
CLEANUP_SCRIPT=""

VARS="USERNAME=${USERNAME},HPC_IP=${HPC_IP},TICKET=${TICKET},PIPELINE_SOURCE=${PIPELINE_SOURCE_DIRECTORY},samples=${INPUTS_STRING},MODULE_OUTPUT_DIRECTORY=${OUTPUT_DIRECTORY}"

# =============================================== ADDITIONAL VARS =====================================================
ADDITIONAL_VARS=""
VARS+="${ADDITIONAL_VARS}"
# =============================================== ADDITIONAL VARS =====================================================

echo Beginning submission log for module: ${MODULE}
JOBID=$(qsub -V -N job-${TICKET} -o ${OUTPUT_DIRECTORY}//${MODULE}_q_submission_output.log -e ${OUTPUT_DIRECTORY}//${MODULE}_q_submission_error.log -wd ${WORKING_DIRECTORY} -v ${VARS} ${SCRIPT_STRING} | cut -d ' ' -f 3);
echo Job Id: ${JOBID}
qsub -V -hold_jid ${JOBID} -N cleanup-${TICKET} -o ${OUTPUT_DIRECTORY}//${MODULE}_module_cleanup.log -e ${OUTPUT_DIRECTORY}//${MODULE}_module_cleanup.log -v USERNAME=${USERNAME},HPC_IP=${HPC_IP},SERVER=${SERVER},TICKET=${TICKET},JOBID=${JOBID} ${CLEANUP_SCRIPT}