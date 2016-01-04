#!/usr/bin/env bash

TICKET = ${TICKET}
JOBID = ${JOBID}

# Notify that cleanup should begin on the webserver
ssh biocis@10.202.64.28 "
    curl --form status="cleanup,${JOBID},${JOB_ID}" 10.202.64.27:190008/messages|"${TICKET}"
"

# Begin calculating stats whilst on a node and inform the webserver

# Remove this file?