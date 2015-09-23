#!/usr/bin/python

import os
import sys

try:
    seqProjectID = sys.argv[1]
except Exception:
    print
    "You must supply a seq project ID"
    sys.exit()

os.environ['PATH'] = '/var/www/biocis/cgi-bin/API'

from runQuery import runQuery

DBquery = "delete from demultiplex where seqProjectID=" + seqProjectID
try:
    res = runQuery(DBquery)
except Exception:
    print
    "unable to remove seq project with ID " + seqProjectID
