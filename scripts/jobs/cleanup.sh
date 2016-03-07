#!/usr/bin/env bash

echo "Begging job cleanup"

ssh ${USERNAME}@${HPC_IP} << EOF
    # Begin calculating stats whilst on a node and inform the webserver
    SUB_TIME=\$(qacct -j ${JOBID} | awk -v i=13 j=2);
    START_TIME=\$(qacct -j ${JOBID} | awk -v i=13 j=3);
    END_TIME=\$(qacct -j ${JOBID} | awk -v i=13 j=4);

    # TODO: Check on the whether the job was dropped

    # Ping back our info to the webserver - TODO: Silence it?
    curl --silent --form "event=module_end" --form "sub=\$SUB_TIME" --form "stat=\$START_TIME" --form "end=\$END_TIME" \$SERVER\'/message/pipeline\|${TICKET}\'
EOF