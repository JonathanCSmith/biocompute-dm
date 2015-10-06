#!/usr/bin/python

import os
import cgi

from API.masterProject import masterProject

docFileStuff = cgi.FieldStorage()

CH = open("../template/cgi_header", "r")
cgiHeader = CH.readlines()
CH.close()

HF = open("../template/template_top.html", "r")
topLines = HF.readlines()
HF.close()

FF = open("../template/template_bottom.html", "r")
bottomLines = FF.readlines()
FF.close()

for ch in cgiHeader:
    print
    ch[:-1]

for he in topLines:
    print
    he[:-1]

print
'<p>fileName', docFileStuff["docFile"].filename, '</p>'
# print '<p>masterProjID', docFileStuff["masterProjID"].value,'</p>'
# print '<p>value',docFileStuff["docFile"].value,'</p>'
# print '<p>docDescription',docFileStuff["docDescription"].value,'</p>'
projectID = docFileStuff["masterProjID"].value

i = masterProject()
i.getMasterProjectByID(projectID)
i.setProjectDir()
i.makeProjectDir()


# write the file to the appropriate project dir

if docFileStuff["docFile"].filename:
    print
    "<p>A</p>"
    # strip leading path from file name to avoid directory traversal attacks
    fn = os.path.basename(docFileStuff["docFile"].filename)
    fileToWriteName = os.path.join(i.projectDir, "Documents", fn)
    print
    "<p>" + fileToWriteName + "</p>"
    open(fileToWriteName, 'wb').write(docFileStuff["docFile"].file.read())

    i.addDocToDB(docFileStuff["docDescription"].value, fileToWriteName)

# send it back to the project page

print
'<script language="javascript" type="text/javascript">'
print
'window.location.href="/cgi-bin/getProjectFromDB?projID=' + projectID + '";'
print
'</script>'

for fo in bottomLines:
    print
    fo[:-1]
