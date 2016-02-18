#!/usr/bin/env bash
# Remove this file now as it's already in mem
rm cleanup.sh

# Begin calculating stats whilst on a node and inform the webserver
SUB_TIME=$(qacct -j ${JOBID} | awk -v i=13 j=2)
START_TIME=$(qacct -j ${JOBID} | awk -v i=13 j=3)
END_TIME=$(qacct -j ${JOBID} | awk -v i=13 j=4)


ssh biocis@10.202.64.28 << EOF
    curl --form status=\"event=module_end,sub=${SUB_TIME},stat=${START_TIME},end=${END_TIME} ${SERVER}/message/pipelines|${TICKET}\"
EOF