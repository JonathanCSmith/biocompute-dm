#!/usr/bin/python

import os,sys

try:
	exptID=sys.argv[1]
except Exception:
	print "You must supply a sequence experiment ID"
	sys.exit()

os.environ['PATH']='/var/www/cgi-bin/coreInSys/API'



from getExptByID import getExptByID
from runQuery import runQuery

import os
try:
	i=getExptByID(exptID)
except Exception:
	print "Unable to find seqExperiment with ID "+exptID
	sys.exit()


print i.seqExptName
for t in range(0,len(i.seqProjs)):
	#print i.seqProjs[t].seqProjectName	

	try:
		i.seqProjs[t].deleteSeqProjectTree()
	except Exception:
		print "Error! Unable to delete seq project data tree for seqProjID: "+str(i.seqProjs[t].seqProjectID)
		sys.exit()


	DBquery="delete from seqExperiment where seqExptID="+str(i.seqExptID)
	try:
		res=runQuery(DBquery)
	except Exception:
		print "unable to remove seq experiment with ID "+str(i.seqExptID)



