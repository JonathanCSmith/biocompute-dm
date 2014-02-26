#!/usr/bin/python


#Can we use the setuid bit to allow a root, passwordless ssh command to one of the nsds to perform the move?




import os,sys

try:
	projID=sys.argv[1]
except Exception:
	print "You must supply a seq project ID and an sftp username"
	sys.exit()


try:
	sftpUser=sys.argv[2]
	sftpPath="/gpfs/home/transfer/"+sftpUser+"/data/"
except Exception:
	print "You must supply a seq project ID and an sftp username"
	sys.exit()



os.environ['PATH']='/var/www/cgi-bin/coreInSys/API'


from seqProject import seqProject
from runQuery import runQuery
from getRunByID import getRunByID
#from seqRun import seqRun

import os
try:
	i=seqProject(0)
	i.getSeqProjByID(int(projID))

except Exception:
	print "Unable to load seq project with ID "+projID
	sys.exit()

DBquery="select location from demultiplex where seqProjectID="+projID

res=runQuery(DBquery)
#print res[0][0]
sourcePATH="/gpfs"+res[0][0]
sourcePATH=sourcePATH+"/processed"

#print sourcePATH
#print i.seqProjectName
#print i.seqRunID

j=getRunByID(int(i.seqRunID))

#print j.seqRunName

sftpPath=sftpPath+j.seqRunName

if os.path.islink(sourcePATH):
	print "Source directory is a symbolic link. Aborting"
	print "(This probably means the data has already been copied to the sftp area)"
	sys.exit()

else:

	command="mv "+sourcePATH+" "+sftpPath

	print command

	linkCommand="ln -s "+sftpPath+" "+sourcePATH

	print linkCommand

"""
try:
	i.deleteSeqProjectTree()
except Exception:
	print "Error! Unable to delete seq project data tree for seqProjID: "+projID
	sys.exit()
"""




