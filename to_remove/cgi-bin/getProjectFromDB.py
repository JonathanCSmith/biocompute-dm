#!/usr/bin/python
import os
import cgi
import cgitb
import API.templateBuilder as writer
cgitb.enable()

link = cgi.FieldStorage()

from API.masterProject import masterProject
from API.getState import getState

writer.writeOutCGI()
writer.writeOutHeader()

if link.has_key('projName'):
    projectName = link["projName"].value

if link.has_key('projID'):
    projectID = link["projID"].value

i = masterProject()
i.getMasterProjectByID(projectID)
projectName = i.projectName

i.appendChildren()

# print i.seqProjs[0].seqRunID


print '<table  border="1" align="center" valign="top" cellpadding="0" cellspacing="0" width=100%>'
print '<tr><th>Project Name</th><th>Project Lead</th><th>Customer ID</th><th>Project opened</th><th>Project last updated</th><th>Project Status</th></tr>'

print '<tr><td><a href="/cgi-bin/getProjectFromDB.cgi?projName=' + i.projectName + '&projID=' + i.masterProjectID + '"> ' + i.projectName + '</a></td><td>' + i.projectLead + '</td><td>Customer</td><td>' + str(
    i.openDate) + '</td><td>' + str(i.lastUpdate) + '</td><td>' + i.status + '</td></tr>'
print '<tr><td colspan=6><p><br>' + i.description + '</p></td></tr>'
print '<tr><td colspan=6><a href="/cgi-bin/editMasterProject?projID=' + i.masterProjectID + '">Edit Master Project Details</a></td></tr>'
# print '<tr><td colspan=6><a href="/cgi-bin/createMasterProjectLink?projID='+i.masterProjectID+'">Create Master Project Link</a></td></tr>'
print '<tr><td colspan=6><hr></td></tr>'
print '<tr><td colspan=6>'
print "<table border=1>"
print '<form  id="unlinkst" name="entitiesToUnlink" method="post" action="/cgi-bin/removeLinks?masterProjID=' + projectID + '&projName=' + projectName + '">'
print '<tr><td colspan=6><h2>Master Project Links</h2>'
print '<p><a href="/cgi-bin/createMasterProjectLink?projID=' + i.masterProjectID + '">Create Master Project Link</a></p>'
print '</td></tr>'

print '<tr><td>Status</td><td>Link Type</td><td>Link Name</td><td>Link ID</td><td></td><td><input id="unLinkButton" type="submit" value="unlink"></td>'
for ch in range(0, len(i.children)):
    print "<tr>"
    status = getState(i.childTypes[ch], i.children[ch].seqProjectID)
    print "<td align='middle'><img src='/images/status_" + status + ".gif' alt='status_" + status + "' ></td>"

    print "<td>", i.childTypes[ch], "</td>"
    if i.childTypes[ch] == "sequencing":
        print "<td>" + i.children[ch].seqProjectName + "</td>"
        print "<td>" + str(i.children[ch].seqProjectID) + "</td>"
        print '<td><a href=/cgi-bin/seqProjDetails?projID=' + str(
            i.children[ch].seqProjectID) + '>details</a></td>'
        print '<td> unlink<input type="checkbox" name="unlink' + str(ch) + '" value="seq_' + str(
            i.children[ch].seqProjectID) + '"></td>'
    print "</tr>"
print "</table>"

print '</td></tr><hr><br><br>'
i.checkDocuments()

# print "<p>",i.documents,"</p>"


print '<tr><td colspan="6">'
print "<table border=1 width=100%>"
print '<tr><td colspan="6">'

print '<h2>Master Project Documents</h2>'
print '<p><a href="/cgi-bin/upLoadDoc.py?masterProjID=' + projectID + '">Upload a document</a></p>'
print '</td></tr>'
if len(i.documents) > 0:
    print '<tr><th>DocumentID</th><th>Document location</th><th>Description</th><th></th></tr>'

for dispDocs in range(0, len(i.documents)):
    projNum = "%05g" % (int(i.masterProjectID))
    projDir = "/link/Projects/" + projNum + "_" + i.projectName + "/Documents"
    fileName = os.path.basename(i.documents[dispDocs][2])
    linkLocation = os.path.join(projDir, fileName)
    print '<tr>'
    # print "<td>"+linkLocation+"</td>"
    # print "</tr><tr>"
    print '<td>', i.documents[dispDocs][0], '</td>'
    print '<td><a href="' + linkLocation + '">', i.documents[dispDocs][2], '</a></td>'
    print '<td>', i.documents[dispDocs][1], '</td>'
    print '<td><a href="/cgi-bin/removeMastProjDoc?masterProjID=' + str(
        i.masterProjectID) + '&documentID=' + str(i.documents[dispDocs][0]) + '">Delete Document</a></td>'
    print '</tr>'

print '</table>'



# print '</td></tr>'

print "</table>"

writer.writeOutFooter()