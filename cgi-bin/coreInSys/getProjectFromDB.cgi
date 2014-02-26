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

from API.masterProject import masterProject
from API.getState import getState


for ch in cgiHeader:
        print ch[:-1]
for he in topLines:
	print he[:-1]

if link.has_key('projName'):
	projectName = link["projName"].value

if link.has_key('projID'):
	projectID=link["projID"].value

i=masterProject()
i.getMasterProjectByID(projectID)
projectName=i.projectName

i.appendChildren()

#print i.seqProjs[0].seqRunID


print '<table  border="1" align="center" valign="top" cellpadding="0" cellspacing="0" width=100%>'
print '<tr><th>Project Name</th><th>Project Lead</th><th>Customer ID</th><th>Project opened</th><th>Project last updated</th><th>Project Status</th></tr>'


print '<tr><td><a href="/cgi-bin/coreInSys/getProjectFromDB.cgi?projName='+i.projectName+'&projID='+i.masterProjectID+'"> '+i.projectName+'</a></td><td>'+i.projectLead+'</td><td>Customer</td><td>'+str(i.openDate)+'</td><td>'+str(i.lastUpdate)+'</td><td>'+i.status+'</td></tr>'
print '<tr><td colspan=6><p><br>'+i.description+'</p></td></tr>'
print '<tr><td colspan=6><a href="/cgi-bin/coreInSys/editMasterProject?projID='+i.masterProjectID+'">Edit Master Project Details</a></td></tr>'
print '<tr><td colspan=6><a href="/cgi-bin/coreInSys/createMasterProjectLink?projID='+i.masterProjectID+'">Create Master Project Link</a></td></tr>'

print '<tr><td colspan=6>'

print "<table border=1>"
print '<form  id="unlinkst" name="entitiesToUnlink" method="post" action="/cgi-bin/coreInSys/removeLinks?masterProjID='+projectID+'&projName='+projectName+'">'
print '<tr><td>Status</td><td>Link Type</td><td>Link Name</td><td>Link ID</td><td></td><td><input id="unLinkButton" type="submit" value="unlink"></td>'
for ch in range(0,len(i.children)):
	print "<tr>"
	status=getState(i.childTypes[ch],i.children[ch].seqProjectID)
	print "<td align='middle'><img src='/CoreInSys/images/status_"+status+".gif' alt='status_"+status+"' ></td>"

	print "<td>",i.childTypes[ch],"</td>"
	if i.childTypes[ch]=="sequencing":
		print "<td>"+i.children[ch].seqProjectName+"</td>"
		print "<td>"+str(i.children[ch].seqProjectID)+"</td>"		
		print '<td><a href=/cgi-bin/coreInSys/seqProjDetails?projID='+str(i.children[ch].seqProjectID)+'>details</a></td>'
		print '<td> unlink<input type="checkbox" name="unlink'+str(ch)+'" value="seq_'+str(i.children[ch].seqProjectID)+'"></td>'
	print "</tr>"
print "</table>"


print '</td></tr>'






print "</table>"

for fo in bottomLines:
	print fo[:-1]




