#!/usr/bin/python

import os
import sys

try:
    runID = sys.argv[1]
except Exception:
    print
    "You must supply a sequence run ID"
    sys.exit()

os.environ['PATH'] = '/var/www/biocis/cgi-bin/API'

from getRunByID import getRunByID
from runQuery import runQuery

try:
    i = getRunByID(runID)
except Exception:
    print
    "Unable to find seqRun with ID " + runID
    sys.exit()

print
i.seqRunName
for t in range(0, len(i.seqProjs)):
    # print i.seqProjs[t].seqProjectName

    try:
        i.seqProjs[t].deleteSeqProjectTree()
    except Exception:
        print
        "Error! Unable to delete seq project data tree for seqProjID: " + str(i.seqProjs[t].seqProjectID)
        sys.exit()

    DBquery = "delete from seqRun where seqRunID=" + str(i.seqRunID)
    try:
        res = runQuery(DBquery)
    except Exception:
        print
        "unable to remove seq run with ID " + str(i.seqRunID)
