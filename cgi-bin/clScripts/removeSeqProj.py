#!/usr/bin/python

import os
import sys

try:
    projID = sys.argv[1]
except Exception:
    print
    "You must supply a seq project ID"
    sys.exit()

os.environ['PATH'] = '/var/www/cgi-bin/coreInSys/API'

from seqProject import seqProject

try:
    i = seqProject(0)
    i.getSeqProjByID(int(projID))

except Exception:
    print
    "Unable to load seq project with ID " + projID
    sys.exit()

print
i.seqProjectName

try:
    i.deleteSeqProjectTree()
except Exception:
    print
    "Error! Unable to delete seq project data tree for seqProjID: " + projID
    sys.exit()
