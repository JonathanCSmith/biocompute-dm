#!/usr/bin/python

import seqRun
import masterProject


#Create a container for the masterProjects
maPr=[]

#Create a seqRun instance


i=seqRun.seqRun()
i.addSeqProject()
i.seqProjs[0].addLaneData()
i.seqProjs[0].lanes[0].addSamples()


#Read the data


IF=open("./testData_1.csv","r")
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
		i.projectName=dataArr[t][1]
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



i.insertSeqRun()
for m in range(0,len(i.seqProjs)):
	#print i.seqProjs[m].seqProjectID		
	i.seqProjs[m].insertSeqProj()
	# At this stage we need to add to the master Project object and then add to the DB
	i.seqProjs[m].getLastSeqProjectID()
	maPr.append(masterProject.masterProject())
        maPr[-1].localProjectID=i.seqProjs[m].seqProjectID
	maPr[-1].projectType="sequencing"
	maPr[-1].projectName=i.seqProjs[m].seqProjectName
	maPr[-1].insertMasterProject()
	#Now we need to add the master project ID to the seqProject object

	#masterProjectID=maPr[-1].getLastMasterProjectID()
	maPr[-1].getLastMasterProjectID()	
	i.seqProjs[m].masterProjectID=maPr[-1].masterProjectID

	#Now we need to update the seqProject database entry
	i.seqProjs[m].updateDB()



	for n in range(0,len(i.seqProjs[m].lanes)):
		#print "\t",i.seqProjs[m].lanes[n].laneNumber
		#print "\t",i.seqProjs[m].lanes[n].sequencingConc
		#print "\t",i.seqProjs[m].lanes[n].PhiXspiked
		i.seqProjs[m].lanes[n].insertLane()


		for o in range(0,len(i.seqProjs[m].lanes[n].samples)):
			#print "\t\t",i.seqProjs[m].lanes[n].samples[o].sampleID
			#print "\t\t",i.seqProjs[m].lanes[n].samples[o].tagID
			#print "\t\t",i.seqProjs[m].lanes[n].samples[o].tagSequence
			#print "\t\t",i.seqProjs[m].lanes[n].samples[o].analysisID
			i.seqProjs[m].lanes[n].samples[o].insertSample()


"""
#seqRunID,	seqProjID,	laneID,	sequencingConc,	PhiXspiked,	SampleID,	tagID,	tagSequence,	analysisID
P211, 		AlProj,		1, 	0.5, 		0.5,		A1,		1,	AAAAAA,		1
P211,           AlProj,         1,      0.5,            0.5,            A2,             2,      GAAAAA,         1
"""






