#!/usr/bin/python




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
'<p><a href="/cgi-bin/sftpSetup?seqProjectID=89">Setup sftp transfer.</a></p>'
# print '<p><a href="/cgi-bin/sftpSetup?seqProjectID=890">Setup sftp transfer.</a></p>'







for fo in bottomLines:
    print
    fo[:-1]
