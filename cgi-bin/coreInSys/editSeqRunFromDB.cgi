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


from API.makeRunObj import makeRunObj
from API.editRunObj2array import runObj2array
from API.getRunByID import getRunByID


for ch in cgiHeader:
        print ch[:-1]
for he in topLines:
	print he[:-1]



remLanesString="var remLanes=[];"
remSamplesString="var remSamples=[];"
remSeqProjsString="var remSeqProjs=[];"

id = link["projID"].value
#id="P211"

i=getRunByID(id)
#print "<p>i.seqRunID",i.seqRunID,"</p>"
expArr=[]
expArr.append(i)
d=runObj2array(expArr)


#print "<p>",d,"</p>"
print '<form id="newSeqRun" name="newSeqRun" method="post" action="/cgi-bin/coreInSys/updateSeqRun">'
#print "<script src='/CoreInSys/js/experimentFormTable.js'></script>"
print "<script src='/CoreInSys/js/runUpdateFormTable.js'></script>"
print "<script>"

print d
print remLanesString
print remSamplesString
print remSeqProjsString

print "experimentTable(runArr);"

print "</script>"
print '<input id="Insert" type="submit" value="Update">'
print '</form>'



for fo in bottomLines:
	print fo[:-1]


