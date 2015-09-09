#!/usr/bin/python

import seqRun


#Create a seqRun instance


i=seqRun.seqRun()
i.addSeqProject()
i.seqProjs[0].addLaneData()
i.seqProjs[0].lanes[0].addSamples()


#Read the data


IF=open("./testData.csv","r")
lines=IF.readlines()
IF.close()

dataArr=[]

for t in lines:
	li=t[:-1].split(",")
	if t[0]!="#":
		for h in range(0,len(li)):
			#li[h]=li[h].replace(" ", "")
			li[h]=li[h].strip()
			#print li[h]
		dataArr.append(li)

#Initialize dictionaries for the insertion methods
seqRunRecord= {'seqRunID':'NULL', 'flowcellID':'NULL', 'startDate':'NULL', 'completionDate':'NULL', 'genomicsLead':'NULL', 'dataLocation':'NULL', 'indexTagCycles':'NULL', 'readCycles':'NULL'}
seqProjRecord= {'seqProjectName':'NULL', 'seqRunID':'NULL', 'customerID':'NULL', 'exptType':'NULL'}
laneRecord= {'seqProjectID':'NULL', 'laneNumber':'NULL', 'sequencingConc':'NULL', 'read1ClusterDensity':'NULL', 'PhiXspiked':'NULL' ,'spike':'NULL', 'spikeRatio':'NULL'}
sampleRecord= {'sampleName':'NULL', 'tagID':'NULL', 'laneID':'NULL', 'tagSequence':'NULL', 'analysisID':'NULL', 'adaptorSequence':'NULL'}		

#Put the data into the seqRun instance

for t in range(0,len(dataArr)):
	if i.seqRunID=="NULL":			#Has the run ID been initialized?	
		i.seqRunID=dataArr[t][0]
			
	if i.seqProjs[-1].seqProjectName=="":			#Has the seqProjectID been initialized?
		i.seqProjs[-1].seqProjectName=dataArr[t][1]
	elif i.seqProjs[-1].seqProjectName!=dataArr[t][1]:	#Do we need to make a new seqProject object within this seqRun?		
		i.addSeqProject()
		i.seqProjs[-1].addLaneData()
		i.seqProjs[-1].lanes[-1].addSamples()
		i.seqProjs[-1].seqProjectName=dataArr[t][1]

	if i.seqProjs[-1].lanes[-1].laneNumber=="NULL":		#Has the laneData been initialized?
		i.seqProjs[-1].lanes[-1].laneNumber=dataArr[t][2]
		i.seqProjs[-1].lanes[-1].sequencingConc=dataArr[t][3]
		i.seqProjs[-1].lanes[-1].PhiXspiked=dataArr[t][4]

	elif i.seqProjs[-1].lanes[-1].laneNumber!=dataArr[t][2]:
		i.seqProjs[-1].addLaneData()
		i.seqProjs[-1].lanes[-1].addSamples()
		i.seqProjs[-1].lanes[-1].laneNumber=dataArr[t][2]
		i.seqProjs[-1].lanes[-1].sequencingConc=dataArr[t][3]
                i.seqProjs[-1].lanes[-1].PhiXspiked=dataArr[t][4]

	if i.seqProjs[-1].lanes[-1].samples[-1].sampleName=="NULL":
		i.seqProjs[-1].lanes[-1].samples[-1].sampleName=dataArr[t][5]
		i.seqProjs[-1].lanes[-1].samples[-1].tagID=dataArr[t][6]
		i.seqProjs[-1].lanes[-1].samples[-1].tagSequence=dataArr[t][7]
		i.seqProjs[-1].lanes[-1].samples[-1].analysisID=dataArr[t][8]
	else:
		i.seqProjs[-1].lanes[-1].addSamples()
		i.seqProjs[-1].lanes[-1].samples[-1].sampleName=dataArr[t][5]
                i.seqProjs[-1].lanes[-1].samples[-1].tagID=dataArr[t][6]
                i.seqProjs[-1].lanes[-1].samples[-1].tagSequence=dataArr[t][7]
                i.seqProjs[-1].lanes[-1].samples[-1].analysisID=dataArr[t][8]


# Convert run objects to an array and send a javascript var as a return.
def expt2Arr (expts):
	jsArray="var runArr=["
	
	for l in range(0,len(expts)):
		jsArray=jsArray+"[["
		jsArray=jsArray+"['seqRunID','"+str(expts[l].seqRunID)+"'],"
		jsArray=jsArray+"['flowcellID','"+str(expts[l].flowcellID)+"'],"
		jsArray=jsArray+"['startDate','"+str(expts[l].startDate)+"'],"
		jsArray=jsArray+"['completionDate','"+str(expts[l].completionDate)+"'],"
		jsArray=jsArray+"['genomicsLead','"+str(expts[l].genomicsLead)+"'],"
		jsArray=jsArray+"['dataLocation','"+str(expts[l].dataLocation)+"'],"
		jsArray=jsArray+"['indexTagCycles','"+str(expts[l].indexTagCycles)+"'],"
		jsArray=jsArray+"['readCycles','"+str(expts[l].readCycles)+"']"
		jsArray=jsArray+"],["

		for m in range(0,len(i.seqProjs)):
			jsArray=jsArray+"["
			jsArray=jsArray+"['seqProjectID','"+str(expts[l].seqProjs[m].seqProjectID)+"'],"
			jsArray=jsArray+"['seqProjectName','"+str(expts[l].seqProjs[m].seqProjectName)+"'],"
			jsArray=jsArray+"['customerID','"+str(expts[l].seqProjs[m].customerID)+"'],"
			jsArray=jsArray+"['exptType','"+str(expts[l].seqProjs[m].exptType)+"']"
			jsArray=jsArray+"],["


			for n in range(0,len(i.seqProjs[m].lanes)):
				jsArray=jsArray+"["
				jsArray=jsArray+"['laneID','"+str(expts[l].seqProjs[m].lanes[n].laneID)+"'],"  
				jsArray=jsArray+"['laneNumber','"+str(expts[l].seqProjs[m].lanes[n].laneNumber)+"'],"
				jsArray=jsArray+"['sequencingConc','"+str(expts[l].seqProjs[m].lanes[n].sequencingConc)+"'],"
				jsArray=jsArray+"['read1ClusterDensity','"+str(expts[l].seqProjs[m].lanes[n].read1ClusterDensity)+"'],"
				jsArray=jsArray+"['PhiXspiked','"+str(expts[l].seqProjs[m].lanes[n].PhiXspiked)+"'],"
				jsArray=jsArray+"['spike','"+str(expts[l].seqProjs[m].lanes[n].spike)+"'],"
				jsArray=jsArray+"['spikeRatio','"+str(expts[l].seqProjs[m].lanes[n].spikeRatio)+"']"			
				jsArray=jsArray+"],["


				for o in range(0,len(i.seqProjs[m].lanes[n].samples)):
					jsArray=jsArray+"["
					jsArray=jsArray+"['sampleID','"+str(expts[l].seqProjs[m].lanes[n].samples[o].sampleID)+"'],"
					jsArray=jsArray+"['sampleName','"+str(expts[l].seqProjs[m].lanes[n].samples[o].sampleName)+"'],"
					jsArray=jsArray+"['tagID','"+str(expts[l].seqProjs[m].lanes[n].samples[o].tagID)+"'],"
					jsArray=jsArray+"['tagSequence','"+str(expts[l].seqProjs[m].lanes[n].samples[o].tagSequence)+"'],"
					jsArray=jsArray+"['analysisID','"+str(expts[l].seqProjs[m].lanes[n].samples[o].analysisID)+"'],"
					jsArray=jsArray+"['adaptorSequence','"+str(expts[l].seqProjs[m].lanes[n].samples[o].adaptorSequence)+"']"
					jsArray=jsArray+"],"


				jsArray=jsArray[:-1]+"],"	#Close the samples part and remove the ","
			jsArray=jsArray[:-1]+"],"		#Close the lane part and remove the ","
		jsArray=jsArray[:-1]+"]]"			#Close the project part and remove the ","
	jsArray=jsArray+"];"			 	#Close the run part and remove the ","
	return(jsArray)

m=[]
m.append(i)

d=expt2Arr(m)

print d

"""
numLeft=0
numRight=0
for lr in range(0,len(d)):
	if d[lr]=="[":
		numLeft=numLeft+1
	if d[lr]=="]":
		numRight=numRight+1
print numLeft, numRight
"""

"""
#seqRunID,	seqProjID,	laneID,	sequencingConc,	PhiXspiked,	SampleID,	tagID,	tagSequence,	analysisID
P211, 		AlProj,		1, 	0.5, 		0.5,		A1,		1,	AAAAAA,		1
P211,           AlProj,         1,      0.5,            0.5,            A2,             2,      GAAAAA,         1
"""






