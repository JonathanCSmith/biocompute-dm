#!/usr/bin/python

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
from API.runObj2array import runObj2array

for ch in cgiHeader:
	print ch[:-1]

for he in topLines:
	print he[:-1]


i=makeRunObj()
expArr=[]
expArr.append(i)

d=runObj2array(expArr)

print "<script src='/CoreInSys/js/experimentTable.js'></script>"

print "<script>"
print d
print "experimentTable(experArr);"
print "</script>"


for fo in bottomLines:
	print fo[:-1]




