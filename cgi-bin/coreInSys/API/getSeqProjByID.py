#file getSeqProjByID.py



def getSeqProjByID(parent,projectID):
	import seqProject
	i=seqProject.seqProject(parent)
	i.getSeqProjByID(projectID)

	#Get the IDs of all of the lanes in this project
	laneDataIDs=i.getLaneDataIDs()
	i.addLaneData()


	#Work down to the lane layer
	for la in range(0,len(laneDataIDs)):
		#add a lane instance
		if la>0 or len(i.lanes)==0:
			if i.lanes[-1].laneID!="NULL":
                              	i.addLaneData()

		i.lanes[-1].getLaneDataByID(laneDataIDs[la][0])

		sampleIDs=i.lanes[-1].getSampleIDs()
		#print sampleIDs
		i.lanes[-1].addSamples()	
		#Now we move into the sample layer
		for sa in range(0,len(sampleIDs)):
			#add a sample instance
			if sa>0 or len(i.lanes[-1].samples)==0:
				if i.lanes[-1].samples[-1].sampleID!=0:
                               		i.lanes[-1].addSamples()

			i.lanes[-1].samples[-1].getSampleDataByID(sampleIDs[sa][0])


	return(i)





