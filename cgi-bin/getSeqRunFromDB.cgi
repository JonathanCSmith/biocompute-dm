#!/usr/bin/python

import sys
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
from API.runObj2array import runObj2array
from API.runOnly2array import runOnly2array
from API.getRunByID import getRunByID


for ch in cgiHeader:
        print ch[:-1]
for he in topLines:
	print he[:-1]


iD = link["projID"].value
#id="P211"
#print "<p>",iD,"</p>"



i=getRunByID(iD)

#for fo in bottomLines:
#        print fo[:-1]
#sys.exit()


"""
print "<p>",len(i.seqProjs),"</p>"
for sqPro in range(0,len(i.seqProjs)):
	print "<p>",i.seqProjs[sqPro].seqProjectName,"</p>"
	print "<p>i.seqProjs[sqPro].seqProjectID",i.seqProjs[sqPro].seqProjectID,"</p>"
"""
expArr=[]
expArr.append(i)

#for fo in bottomLines:
#	print fo[:-1]
#sys.exit()





#print "<p>len(i.seqProjs)",len(i.seqProjs),"</p>"
#print "<p>i.seqProjs[0].seqProjectID",i.seqProjs[0].seqProjectID,"</p>"

#print "<p>",d,"</p>"

if len(i.seqProjs)==1 and i.seqProjs[0].seqProjectID==0:



	d=runOnly2array(expArr)



	print '<a href="/cgi-bin/coreInSys/linkSeqProj2Run?runID='+iD+'">Link Seq Project</a>'
	#print "<p>",expArr[0].seqRunName,"</p>" 
	print "<script src='/CoreInSys/js/runOnlyTable.js'></script>"
	print "<script>"
        print d
        print "experOnlyTable(runArr);"
        print "</script>"



else:


	d=runObj2array(expArr)

	


	print '<a href="/cgi-bin/coreInSys/linkSeqProj2Run?runID='+iD+'">Link Seq Project</a>'
	print "<script src='/CoreInSys/js/runTable.js'></script>"
	print "<script>"
	print d
	print "experimentTable(runArr);"
	print "</script>"

	print '<a href="/cgi-bin/coreInSys/editSeqRunFromDB.cgi?projID='+iD+'"> EDIT</a>'
	print '<br>'
	print '<a href="/CoreInSys/staging/'+i.dataLocation+'/First_Base_Report.htm" target="_blank">SEQUENCING STATUS</a>'

print '<br>'


for fo in bottomLines:
	print fo[:-1]




