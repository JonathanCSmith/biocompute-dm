#!/usr/bin/python

import os
import cgi
# import re
# import cgitb; cgitb.enable()
from API.runObj2array import runObj2array
from API.getRunByID import getRunByID
from API.seqRun import seqRun
from API.getTagSequence import getTagSequence

csvFile = cgi.FieldStorage()

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

# print '<p>fileName',csvFile["seqDataFile"].filename,'</p>'

dataArr = []
csvForm = {}

if csvFile.has_key("seqDataFile"):
    CF = csvFile["seqDataFile"]
    lines = CF.file.readlines()

elif csvFile.has_key("csv"):
    # print csvFile["csv"].value
    lin = csvFile["csv"].value
    lines = lin.split("\r\n")
# print lines

for delLin in range(len(lines) - 1, -1, -1):
    lines[delLin] = lines[delLin].strip()
    if lines[delLin] == "":
        del lines[delLin]


# Put the csv data into an array. We'll need to come back and apply some sanity checks to the data later.

# li=lines[0][:-1].split(",")
li = lines[0].split(",")

if lines[0][0] != "#":
    for h in range(0, len(li)):
        li[h] = li[h].strip()
        csvForm[li[h]] = h
    default_form = 0


# dataArr.append(li)

for t in range(1, len(lines)):
    # li=lines[t][:-1].split(",")
    li = lines[t].split(",")
    if lines[t][0] != "#":
        for h in range(0, len(li)):
            # li[h]=li[h].replace(" ", "")
            li[h] = li[h].strip()
            if li[h] == "" and len(
                    dataArr) > 0:  # Here we try and fill in any blanks in the csv file with the entry in the cell above.

                li[h] = dataArr[-1][h]
        dataArr.append(li)
# print "<p>dataArr[0]",dataArr[0],"</p>"
# print "<p>dataArr[1]",dataArr[1],"</p>"

"""
0       Well ,
1       Customer Sample ID,
2       P.I.,
3       Index kit I.D.,
4       Index ID,
5       Pool I.D.,
6       Hi-Seq/Flow Cell Project number,
7       Library Prep Operator,
8       Library Prep Start Date,
9       Library Prep End Date,
10      Samples sent for Sequencing,
11      Sequencing Start Date,
12      Sequencing End Date,
13      Sequencing Operator,
14      Lane,
15      Flow Cell I.D.,
16      # of Index Tag Cycles,
17      Sequencing Concentration
"""

# headerList=["Ignore","Well","Customer Sample ID","Index sequence","Sequencing Project Title","Index kit I.D.","Index ID","Pool I.D.","Hi-Seq/Flow Cell Project number","Library Prep Operator","Library Prep Start Date","Library Prep End Date","Samples sent for Sequencing","Sequencing Start Date","Sequencing End Date","Sequencing Operator","Lane","Flow Cell I.D.","# of Index Tag Cycles","Sequencing Concentration"]

headerList = ["Ignore", "Customer Sample ID", "Index sequence", "Index 2 sequence", "Sequencing Project Title",
              "Index kit I.D.", "Index ID", "Lane", "Sequencing Concentration", "Read1 Cluster Density", "PhiXspiked",
              "spike", "spikeRatio", "adaptorSequence"]


def optionsList(columnNumber, header, csvForm, headerList):
    if header == "P.I.":
        header = "Sequencing Project Title"

    possHeaders = headerList[:]
    chosen = "Ignore"
    print
    "<select name='colHeader_" + str(columnNumber) + "'>"
    for l in range(0, len(possHeaders)):
        if possHeaders[l] == header:
            chosen = possHeaders[l]
            del possHeaders[l]
            break
    print
    "<option value='" + chosen + "'>" + chosen + "</option>"
    for m in range(0, len(possHeaders)):
        print
        "<option  value='" + possHeaders[m] + "'>" + possHeaders[m] + "</option>"
        # print "<option name='colHeader_"+str(columnNumber)+"'>"+possHeaders[m]+"</option>"
    print
    "</select>"


print
"<br>"
print
"<h3>Please select table headings and press confirm.</h3>"
print
'<form id="headerConfirm" name="headerConfirm" method="post" action="/cgi-bin/upLoaderSeqProjStage2">'
print
'<input id="Insert" type="submit" value="Confirm">'
print
"<table border=1>"
print
"<tr>"
li = lines[0].split(",")

for n in range(0, len(li)):
    print
    "<td>" + li[n] + "</td>"
print
"</tr>"
print
"<tr>"

for n in range(0, len(li)):
    print
    "<td>"
    optionsList(n, li[n], csvForm, headerList)
    print
    "</td>"
print
"</tr>"

for o in range(0, len(dataArr)):
    print
    "<tr>"
    for p in range(0, len(dataArr[o])):
        print
        "<td><input type='hidden' name='dataArr_" + str(o) + "_" + str(p) + "' value='" + dataArr[o][p] + "'>" + \
        dataArr[o][p] + "</td>"
    print
    "</tr>"
print
"</table>"
print
'<input id="Insert" type="submit" value="Confirm">'
print
"</form>"

for fo in bottomLines:
    print
    fo[:-1]