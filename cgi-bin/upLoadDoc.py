#!/usr/bin/python


import cgi

link = cgi.FieldStorage()

CH = open("../template/cgi_header", "r")
cgiHeader = CH.readlines()
CH.close()

HF = open("../template/template_top.html", "r")
topLines = HF.readlines()
HF.close()

FF = open("../template/template_bottom.html", "r")
bottomLines = FF.readlines()
FF.close()

# masterProjID

if link.has_key('masterProjID'):
    masterProjID = link["masterProjID"].value

for ch in cgiHeader:
    print
    ch[:-1]

for he in topLines:
    print
    he[:-1]

print
"<p>Upload master project document:</p>"
print
'<form name="input" enctype="multipart/form-data" action="/cgi-bin/upLoaderDocStage1.py?masterProjID=' + masterProjID + '" method="post">'
print
'<p>Document Description<br>'
print
'<textarea name="docDescription" width="200"  height="5" > </textarea></p>'

print
'<p><input type="file" name="docFile" size="40"><br> <input type="submit" Value="Upload" ></p>'
print
'</form>'

for fo in bottomLines:
    print
    fo[:-1]