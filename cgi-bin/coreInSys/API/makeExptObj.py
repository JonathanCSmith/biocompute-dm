# file makeExptObj.py


import seqExpt


#print "Running"


#Create a seqExperiment instance

def makeExptObj():

	i=seqExpt.seqExpt()
	i.addSeqProject()
	i.seqProjs[0].addLaneData()
	i.seqProjs[0].lanes[0].addSamples()


	#Read the data


	IF=open("/var/www/cgi-bin/coreInSys/API/testData.csv","r")
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
	seqExptRecord= {'seqExptID':'NULL', 'flowcellID':'NULL', 'startDate':'NULL', 'completionDate':'NULL', 'genomicsLead':'NULL', 'dataLocation':'NULL', 'indexTagCycles':'NULL', 'readCycles':'NULL'}
	seqProjRecord= {'seqProjectName':'NULL', 'seqExptID':'NULL', 'customerID':'NULL'}
	laneRecord= {'seqProjectID':'NULL', 'laneNumber':'NULL', 'sequencingConc':'NULL', 'read1ClusterDensity':'NULL', 'PhiXspiked':'NULL' ,'spike':'NULL', 'spikeRatio':'NULL'}
	sampleRecord= {'sampleName':'NULL', 'tagID':'NULL', 'laneID':'NULL', 'tagSequence':'NULL', 'analysisID':'NULL', 'adaptorSequence':'NULL'}		

	#Put the data into the seqExperiment instance

	for t in range(0,len(dataArr)):
		if i.seqExptID=="NULL":			#Has the experiment ID been initialized?	
			i.seqExptID=dataArr[t][0]
			
		if i.seqProjs[-1].seqProjectName=="":			#Has the seqProjectID been initialized?
			i.seqProjs[-1].seqProjectName=dataArr[t][1]
		elif i.seqProjs[-1].seqProjectName!=dataArr[t][1]:	#Do we need to make a new seqProject object within this seqExperiment?		
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

	return(i)






