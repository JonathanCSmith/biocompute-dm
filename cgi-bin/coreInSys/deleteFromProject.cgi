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
from API.masterProjSeq2array import masterProjSeq2array
from API.getProjByID import getProjByID


for ch in cgiHeader:
        print ch[:-1]
for he in topLines:
	print he[:-1]

projectType=link["projType"].value


projectName = link["projName"].value

projectID=link["projID"].value


#id="P211"
i=getProjByID(projectID,projectType)
#print i.seqProjs[0].seqRunID
if i.seqProjs[0].seqRunID=="NULL":
	print "<p>This project isn't associated with any experiment.</p>"	
	print '<p><a href="/cgi-bin/coreInSys/linkProject?projName='+projectName+'&projType='+projectType+'&projID='+projectID+'">Click here</a> to link this project to an experiment</p>'


else:

	projArr=[]
	projArr.append(i)


	d=masterProjSeq2array(projArr)

	#print "<p>",d,"</p>"


	print "<script src='/CoreInSys/js/projectTableDelete.js'></script>"
	print "<script>"
	print d
	print "experimentTable(experArr);"
	print "</script>"

	print '<a href="/cgi-bin/coreInSys/editProjectFromDB.cgi?projName='+projectName+'&projType='+projectType+'&projID='+projectID+'"> EDIT</a>'



for fo in bottomLines:
	print fo[:-1]




