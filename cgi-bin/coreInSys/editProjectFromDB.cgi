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


from API.makeExptObj import makeExptObj
from API.masterProjSeq2array import masterProjSeq2array
from API.getProjByID import getProjByID


for ch in cgiHeader:
        print ch[:-1]
for he in topLines:
	print he[:-1]

projectType=link["projType"].value

#print "<p>"+projectType+"</p>"
projectName = link["projName"].value
#print "<p>"+projectName+"</p>"
projectID=link["projID"].value
#print "<p>"+projectID+"</p>"

#id="P211"
i=getProjByID(projectID,projectType)


projArr=[]
projArr.append(i)
#print projArr

d=masterProjSeq2array(projArr)

print '<form id="upDateForm" name="upDate_Form" method="post" action="/cgi-bin/coreInSys/updateProject">'
print "<script src='/CoreInSys/js/projectFormTable.js'></script>"
print "<script>"
print d
print "experimentTable(experArr);"
print "</script>"
print '<input id="Update" type="submit" value="Update">'
print '</form>'

for fo in bottomLines:
	print fo[:-1]




