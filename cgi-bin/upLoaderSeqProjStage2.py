#!/usr/bin/python

import os
import cgi
import sys

# import re


# import cgitb; cgitb.enable()
from API.runObj2array import runObj2array
from API.getRunByID import getRunByID
from API.seqRun import seqRun
from API.getTagSequence import getTagSequence

form = cgi.FieldStorage()

# Illegal characters ? ( ) [ ] / \ = + < > : ; " ' , * ^ | & .
illChars = set("?()[]/\=+<>:;\"',*^|&.")

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
"<p>upLoaderSeqProjStage2</p>"

formKeys = form.keys()

dataArr = []
dAA = []
dAB = []

for u in range(0, len(formKeys)):
    if formKeys[u][:7] == "dataArr":
        dA = formKeys[u].split("_")
        dAA.append(int(dA[1]))
        dAB.append(int(dA[2]))

dAA.sort()
dAB.sort()

dimA = dAA[-1] + 1
dimB = dAB[-1] + 1

for v in range(0, dimA):
    dataArr.append([])

for v in range(0, dimA):
    for w in range(0, dimB):
        dataArr[v].append(0)

for u in range(0, len(formKeys)):
    if formKeys[u][:7] == "dataArr":
        dA = formKeys[u].split("_")
        A = int(dA[1])
        B = int(dA[2])
        dataArr[A][B] = form[formKeys[u]].value

try:
    csvForm = {}
    for w in range(0, len(formKeys)):
        if formKeys[w][:9] == "colHeader":
            dA = formKeys[w].split("_")
            csvForm[form[formKeys[w]].value] = int(dA[1])

except Exception:
    print
    "<p>Unable to read form data.</p>"
    for fo in bottomLines:
        print
        fo[:-1]
    sys.exit()

try:
    # Create a seqRun instance
    i = []

    i.append(seqRun())
    i[-1].addSeqProject()
    i[-1].seqProjs[0].addLaneData()
    i[-1].seqProjs[0].lanes[0].addSamples()

except Exception:
    print
    "<p>Unable to make seqRun instance seqExpt object tree</p>"
    for fo in bottomLines:
        print
        fo[:-1]
    sys.exit()

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
# Put the data into the seqRun instance
# default_form=1


# csvForm={}
"""
if default_form==1:
        csvForm["Well"]=0
        csvForm["Customer Sample ID"]=1                         #
        csvForm["P.I."]=2                                       #
        csvForm["Index kit I.D."]=3                             #
        csvForm["Index ID"]=4                                   #
        csvForm["Pool I.D."]=5
        csvForm["Hi-Seq/Flow Cell Project number"]=6            #
        csvForm["Library Prep Operator"]=7
        csvForm["Library Prep Start Date"]=8
        csvForm["Library Prep End Date"]=9
        csvForm["Samples sent for Sequencing"]=10
        csvForm["Sequencing Start Date"]=11                     #
        csvForm["Sequencing End Date"]=12                       #
        csvForm["Sequencing Operator"]=13                       #
        csvForm["Lane"]=14                                      #
        csvForm["Flow Cell I.D."]=15                            #
        csvForm["# of Index Tag Cycles"]=16
        csvForm["Sequencing Concentration"]=17                  #
	#csvForm["Index sequence"]=18
"""

try:
    for t in range(0, len(dataArr)):
        # for t in range(0,1):	#take only the first seqProject
        if csvForm.has_key("Sequencing Project Title"):
            if i[-1].seqProjs[-1].seqProjectName == "SeqProjName":  # Has the seqProjectID been initialized?
                i[-1].seqProjs[-1].seqProjectName = dataArr[t][csvForm["Sequencing Project Title"]]
                if csvForm.has_key("Customer Sample ID"):
                    i[-1].projectName = dataArr[t][csvForm["Customer Sample ID"]]

                    # Illegal characters ? ( ) [ ] / \ = + < > : ; " ' , * ^ | & .
                    illChars = set("?()[]/\=+<>:;\"',*^|&.")
                    if any((c in illChars) for c in dataArr[t][csvForm["Customer Sample ID"]]):
                        print
                        "<p>Illegal character in sample name: " + dataArr[t][csvForm["Customer Sample ID"]] + "</p>"
                        raise Exception

                    else:
                        i[-1].projectName = dataArr[t][csvForm["Customer Sample ID"]]



            elif i[-1].seqProjs[-1].seqProjectName != dataArr[t][csvForm[
                "Sequencing Project Title"]]:  # Do we need to make a new seqProject object within this seqRun?
                # if we have a more than one seqProject here we want to take only the first so we break here.
                break
                """
                i[-1].addSeqProject()
                i[-1].seqProjs[-1].addLaneData()
                i[-1].seqProjs[-1].lanes[-1].addSamples()
                i[-1].seqProjs[-1].seqProjectName=dataArr[t][csvForm["Sequencing Project Title"]]
                """


                # else:
                # print "<p>No data column named 'Sequencing Project Title'</p>"
                # break

        if csvForm.has_key("Lane") == 0:
            print
            "<p>No data column named 'Lane'</p>"
            # break

        if i[-1].seqProjs[-1].lanes[-1].laneNumber == "NULL":  # Has the laneData been initialized?
            i[-1].seqProjs[-1].lanes[-1].laneNumber = dataArr[t][csvForm["Lane"]]

            if csvForm.has_key("Sequencing Concentration"):
                dataArr[t][csvForm["Sequencing Concentration"]].strip()
                if dataArr[t][csvForm["Sequencing Concentration"]][-2:] == "pM":
                    dataArr[t][17] = dataArr[t][csvForm["Sequencing Concentration"]][:-2]
                i[-1].seqProjs[-1].lanes[-1].sequencingConc = dataArr[t][csvForm["Sequencing Concentration"]]

                # if i[-1].seqProjs[-1].lanes[-1].PhiXspiked=="NULL":
            if csvForm.has_key("PhiXspiked"):
                try:
                    float(dataArr[t][csvForm["PhiXspiked"]])
                    i[-1].seqProjs[-1].lanes[-1].PhiXspiked = dataArr[t][csvForm["PhiXspiked"]]
                except ValueError:
                    i[-1].seqProjs[-1].lanes[-1].PhiXspiked = "NULL"


                    # if i[-1].seqProjs[-1].lanes[-1].read1ClusterDensity==0:
            if csvForm.has_key("Read1 Cluster Density"):
                i[-1].seqProjs[-1].lanes[-1].read1ClusterDensity = dataArr[t][csvForm["Read1 Cluster Density"]]

                # if i[-1].seqProjs[-1].lanes[-1].spike=="NULL":
            if csvForm.has_key("spike"):
                i[-1].seqProjs[-1].lanes[-1].spike = dataArr[t][csvForm["spike"]]

                # if i[-1].seqProjs[-1].lanes[-1].spikeRatio==0.0:
            if csvForm.has_key("spikeRatio"):
                i[-1].seqProjs[-1].lanes[-1].spikeRatio = dataArr[t][csvForm["spikeRatio"]]




                # i.seqProjs[-1].lanes[-1].PhiXspiked=dataArr[t][4]

        elif i[-1].seqProjs[-1].lanes[-1].laneNumber != dataArr[t][csvForm["Lane"]]:
            i[-1].seqProjs[-1].addLaneData()
            i[-1].seqProjs[-1].lanes[-1].addSamples()
            i[-1].seqProjs[-1].lanes[-1].laneNumber = dataArr[t][csvForm["Lane"]]
            # Nasty. Need to find a better way to extract the value of this in case the units/number format is changed
            if csvForm.has_key("Sequencing Concentration"):
                dataArr[t][csvForm["Sequencing Concentration"]].strip()
                if dataArr[t][csvForm["Sequencing Concentration"]][-2:] == "pM":
                    dataArr[t][csvForm["Sequencing Concentration"]] = dataArr[t][csvForm["Sequencing Concentration"]][
                                                                      :-2]
                i[-1].seqProjs[-1].lanes[-1].sequencingConc = dataArr[t][csvForm["Sequencing Concentration"]]

                # if i[-1].seqProjs[-1].lanes[-1].PhiXspiked=="NULL":
            if csvForm.has_key("PhiXspiked"):
                try:
                    float(dataArr[t][csvForm["PhiXspiked"]])
                    i[-1].seqProjs[-1].lanes[-1].PhiXspiked = dataArr[t][csvForm["PhiXspiked"]]
                except ValueError:
                    i[-1].seqProjs[-1].lanes[-1].PhiXspiked = "NULL"

                    # if i[-1].seqProjs[-1].lanes[-1].read1ClusterDensity==0:
            if csvForm.has_key("Read1 Cluster Density"):
                i[-1].seqProjs[-1].lanes[-1].read1ClusterDensity = dataArr[t][csvForm["Read1 Cluster Density"]]

                # if i[-1].seqProjs[-1].lanes[-1].spike=="NULL":
            if csvForm.has_key("spike"):
                i[-1].seqProjs[-1].lanes[-1].spike = dataArr[t][csvForm["spike"]]

                # if i[-1].seqProjs[-1].lanes[-1].spikeRatio==0.0:
            if csvForm.has_key("spikeRatio"):
                i[-1].seqProjs[-1].lanes[-1].spikeRatio = dataArr[t][csvForm["spikeRatio"]]






                # i.seqProjs[-1].lanes[-1].PhiXspiked=dataArr[t][4]

        if csvForm.has_key("Customer Sample ID") == 0:
            print
            "<p>No column named 'Customer Sample ID'</p>"
            # break

        if i[-1].seqProjs[-1].lanes[-1].samples[-1].sampleName == "NULL":

            # Stop if sample ID has any illegal characters ? ( ) [ ] / \ = + < > : ; " ' , * ^ | & .
            if any((c in illChars) for c in dataArr[t][csvForm["Customer Sample ID"]]):
                print
                "<p><b>Illegal character in sample name: " + dataArr[t][csvForm["Customer Sample ID"]] + "</b></p>"
                raise Exception

            else:
                i[-1].seqProjs[-1].lanes[-1].samples[-1].sampleName = dataArr[t][csvForm["Customer Sample ID"]]


                # i[-1].seqProjs[-1].lanes[-1].samples[-1].sampleName=dataArr[t][csvForm["Customer Sample ID"]]

            if csvForm.has_key("Index sequence"):
                indexSequence = dataArr[t][csvForm["Index sequence"]]
                if csvForm.has_key("Index 2 sequence"):
                    if dataArr[t][csvForm[
                        "Index 2 sequence"]] != "Ignore":  ########################################Rev comp this is certain circumstances
                        indexSequence = indexSequence + "-" + dataArr[t][csvForm["Index 2 sequence"]]
                i[-1].seqProjs[-1].lanes[-1].samples[-1].tagSequence = indexSequence
                # i[-1].seqProjs[-1].lanes[-1].samples[-1].tagSequence=dataArr[t][csvForm["Index sequence"]]




            else:
                if csvForm.has_key("Index ID"):
                    i[-1].seqProjs[-1].lanes[-1].samples[-1].tagID = dataArr[t][csvForm["Index ID"]]
                if csvForm.has_key("Index kit I.D."):
                    i[-1].seqProjs[-1].lanes[-1].samples[-1].tagKit = dataArr[t][csvForm["Index kit I.D."]]
                    if csvForm.has_key("Index ID"):
                        seq = getTagSequence(dataArr[t][csvForm["Index kit I.D."]], dataArr[t][csvForm["Index ID"]])
                        i[-1].seqProjs[-1].lanes[-1].samples[-1].tagSequence = seq
                if csvForm.has_key("adaptorSequence"):
                    i[-1].seqProjs[-1].lanes[-1].samples[-1].adaptorSequence = dataArr[t][csvForm["adaptorSequence"]]

        # i.seqProjs[-1].lanes[-1].samples[-1].analysisID=dataArr[t][8]
        else:
            i[-1].seqProjs[-1].lanes[-1].addSamples()

            # Stop if sample ID has any illegal characters ? ( ) [ ] / \ = + < > : ; " ' , * ^ | & .
            if any((c in illChars) for c in dataArr[t][csvForm["Customer Sample ID"]]):
                print
                "<p><b>Illegal character in sample name: " + dataArr[t][csvForm["Customer Sample ID"]] + "</b></p>"
                raise Exception

            else:
                i[-1].seqProjs[-1].lanes[-1].samples[-1].sampleName = dataArr[t][csvForm["Customer Sample ID"]]
            # i[-1].seqProjs[-1].lanes[-1].samples[-1].sampleName=dataArr[t][csvForm["Customer Sample ID"]]
            if csvForm.has_key("Index sequence"):
                indexSequence = dataArr[t][csvForm["Index sequence"]]
                if csvForm.has_key("Index 2 sequence"):
                    if dataArr[t][csvForm[
                        "Index 2 sequence"]] != "Ignore":  ########################################Rev comp this is certain circumstances
                        indexSequence = indexSequence + "-" + dataArr[t][csvForm["Index 2 sequence"]]
                i[-1].seqProjs[-1].lanes[-1].samples[-1].tagSequence = indexSequence
                # i[-1].seqProjs[-1].lanes[-1].samples[-1].tagSequence=dataArr[t][csvForm["Index sequence"]]

            else:
                if csvForm.has_key("Index ID"):
                    i[-1].seqProjs[-1].lanes[-1].samples[-1].tagID = dataArr[t][csvForm["Index ID"]]
                if csvForm.has_key("Index kit I.D."):
                    i[-1].seqProjs[-1].lanes[-1].samples[-1].tagKit = dataArr[t][csvForm["Index kit I.D."]]
                    if csvForm.has_key("Index ID"):
                        seq = getTagSequence(dataArr[t][csvForm["Index kit I.D."]], dataArr[t][csvForm["Index ID"]])
                        i[-1].seqProjs[-1].lanes[-1].samples[-1].tagSequence = seq
                if csvForm.has_key("adaptorSequence"):
                    i[-1].seqProjs[-1].lanes[-1].samples[-1].adaptorSequence = dataArr[t][csvForm["adaptorSequence"]]

                    # i.seqProjs[-1].lanes[-1].samples[-1].analysisID=dataArr[t][8]

except Exception:
    print
    "<p>Unable to read data into seqExpt object tree</p>"
    for fo in bottomLines:
        print
        fo[:-1]
    sys.exit()



# expArr=[]
# expArr.append(projs)
# d=runObj2array(expArr)
try:
    d = runObj2array(i)
except Exception:
    print
    "<p>Unable to make javaScript data array</p>"
    for fo in bottomLines:
        print
        fo[:-1]
    sys.exit()



# print "<p>",d,"</p>"
print
'<form id="newSeqRun" name="newSeqRun" method="post" action="/cgi-bin/insertSeqProject">'
print
"<script src='/js/projectUploadTable.js'></script>"
print
"<script>"
print
d
print
"experimentTable(runArr);"
print
"</script>"
print
'<input id="Insert" type="submit" value="Insert">'
print
'</form>'

for fo in bottomLines:
    print
    fo[:-1]
