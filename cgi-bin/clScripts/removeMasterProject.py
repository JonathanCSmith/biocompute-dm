#!/usr/bin/python

import os
import sys

try:
    projID = sys.argv[1]
except Exception:
    print
    "You must supply a master project ID"
    sys.exit()

os.environ['PATH'] = '/var/www/cgi-bin/coreInSys/API'

from runQuery import runQuery

DBquery = "delete from masterProject where masterProjectID=" + projID
try:
    res = runQuery(DBquery)
except:
    print
    "Unable to remove master project."
    sys.exit()

DBquery = "delete from typeLinker where parentID=" + projID
try:
    res = runQuery(DBquery)
except:
    print
    "Unable to remove master project link."
    sys.exit()
