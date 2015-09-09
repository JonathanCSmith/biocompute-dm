#!/usr/bin/python


import cgi
link = cgi.FieldStorage()


CH=open("./cgi_header","r")
cgiHeader=CH.readlines()
CH.close()

HF=open("./template_top","r")
topLines=HF.readlines()
HF.close()

FF=open("./template_bottom","r")
bottomLines=FF.readlines()
FF.close()

#from API.makeRunObj import makeRunObj
from API.getSeqProjByID import getSeqProjByID
from API.getRunByID import getRunByID



for ch in cgiHeader:
        print ch[:-1]
for he in topLines:
	print he[:-1]



#This is the seqProjectID
projectID=link["projID"].value


i=getSeqProjByID(0,projectID)
#print i.seqProjs[0].seqRunID
if i.seqRunID=="NULL":
	print "<p>This project isn't associated with any experiment.</p>"	
	print '<p><a href="/cgi-bin/coreInSys/linkProject?projName='+projectName+'&projType='+projectType+'&projID='+projectID+'">Click here</a> to link this project to an experiment</p>'
else:

	#get experimental details

	j=getRunByID(i.seqRunID)
	projectName=j.seqRunName

	#print "<p>",j.flowcellID,"</p>"
	#print "<p>",j.dataLocation,"</p>"
	#print "<p>",j.indexTagCycles,"</p>"
	#print "<p>",j.readCycles,"</p>"

	from PL_CONTROL.postAlignQC import postAlignQC
	pq=postAlignQC(i,j)
	#print "<p>fq.QCreportsDir",fq.QCreportsDir,"</p>"
	#print "<p>fq.fastQCdir",fq.fastQCdir,"</p>"
	#fq.makeDir(fq.QCreportsDir)
        #fq.makeDir(fq.fastQCdir)      
	pq.establishPaths()
	pq.setupAllSamples()
	
print '<script language="javascript" type="text/javascript">'
print 'window.location.href="/cgi-bin/coreInSys/seqProjDetails?projID='+projectID+'";'
print '</script>'


for fo in bottomLines:
	print fo[:-1]




