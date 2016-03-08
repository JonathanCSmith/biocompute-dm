#!/usr/bin/env bash

echo "Begging job cleanup"

USERNAME=
HPC_IP=
TICKET=
JOBID=
SERVER=

ssh ${USERNAME}@${HPC_IP} << EOF
    # Begin calculating stats whilst on a node and inform the webserver
    SUB_TIME=\$(qacct -j ${JOBID} | awk 'FNR == 13 {print \$2 " " \$3 " "  \$4 " " \$5 " " \$6}');
    START_TIME=\$(qacct -j ${JOBID} | awk 'FNR == 14 {print \$2 " " \$3 " " \$4 " " \$5 " " \$6}');
    END_TIME=\$(qacct -j ${JOBID} | awk 'FNR == 15 {print \$2 " " \$3 " " \$4 " " \$5 " " \$6}');

    # Debug
    echo "Sub Time: \${SUB_TIME}"
    echo "Start Time: \${START_TIME}"
    echo "End Time: \${END_TIME}"

    # TODO: Check on the whether the job was dropped

    # Ping back our info to the webserver
    ADDRESS="${SERVER}/message/pipelines|${TICKET}"
    echo \${ADDRESS}
    curl --form "event=module_end" --form "sub=\$SUB_TIME" --form "start=\$START_TIME" --form "end=\$END_TIME" "\${ADDRESS}"

EOF