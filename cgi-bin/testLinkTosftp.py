#!/usr/bin/python




CH = open("./cgi_header", "r")
cgiHeader = CH.readlines()
CH.close()

HF = open("./template_top.html", "r")
topLines = HF.readlines()
HF.close()

FF = open("./template_bottom.html", "r")
bottomLines = FF.readlines()
FF.close()

for ch in cgiHeader:
    print
    ch[:-1]

for he in topLines:
    print
    he[:-1]

print
'<p><a href="/cgi-bin/coreInSys/sftpSetup?seqProjectID=89">Setup sftp transfer.</a></p>'
# print '<p><a href="/cgi-bin/coreInSys/sftpSetup?seqProjectID=890">Setup sftp transfer.</a></p>'







for fo in bottomLines:
    print
    fo[:-1]
