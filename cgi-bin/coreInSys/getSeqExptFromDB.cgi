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


#from API.makeExptObj import makeExptObj
from API.exptObj2array import exptObj2array
from API.exptOnly2array import exptOnly2array
from API.getExptByID import getExptByID


for ch in cgiHeader:
        print ch[:-1]
for he in topLines:
	print he[:-1]


iD = link["projID"].value
#id="P211"
#print "<p>",iD,"</p>"

i=getExptByID(iD)
"""
print "<p>",len(i.seqProjs),"</p>"
for sqPro in range(0,len(i.seqProjs)):
	print "<p>",i.seqProjs[sqPro].seqProjectName,"</p>"
	print "<p>",i.seqProjs[sqPro].seqProjectID,"</p>"
"""
expArr=[]
expArr.append(i)


#print "<p>len(i.seqProjs)",len(i.seqProjs),"</p>"
#print "<p>i.seqProjs[0].seqProjectID",i.seqProjs[0].seqProjectID,"</p>"

#print "<p>",d,"</p>"

if len(i.seqProjs)==1 and i.seqProjs[0].seqProjectID==0:

	d=exptOnly2array(expArr)
	#print "<p>",d,"</p>"
	print '<a href="/cgi-bin/coreInSys/linkSeqProj2Expt?exptID='+iD+'">Link Seq Project</a>'
	#print "<p>",expArr[0].seqExptName,"</p>" 
	print "<script src='/CoreInSys/js/experOnlyTable.js'></script>"
	print "<script>"
        print d
        print "experOnlyTable(experArr);"
        print "</script>"



else:

	d=exptObj2array(expArr)

	print '<a href="/cgi-bin/coreInSys/linkSeqProj2Expt?exptID='+iD+'">Link Seq Project</a>'
	print "<script src='/CoreInSys/js/experimentTable.js'></script>"
	print "<script>"
	print d
	print "experimentTable(experArr);"
	print "</script>"

	print '<a href="/cgi-bin/coreInSys/editSeqExptFromDB.cgi?projID='+iD+'"> EDIT</a>'
	print '<br>'
	print '<a href="/CoreInSys/staging/'+i.dataLocation+'/First_Base_Report.htm" target="_blank">SEQUENCING STATUS</a>'

print '<br>'


for fo in bottomLines:
	print fo[:-1]




