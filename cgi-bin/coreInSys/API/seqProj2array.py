#file seqProj2array.py



# Convert run objects to an array and send a javascript var as a return.
def seqProj2array(expts):
	jsArray="var runArr=["
	for m in range(0,len(expts)):

		jsArray=jsArray+"[["
		jsArray=jsArray+"['seqProjectID','"+str(expts[m].seqProjectID)+"'],"
		jsArray=jsArray+"['seqProjectName','"+str(expts[m].seqProjectName)+"'],"
		jsArray=jsArray+"['masterProjectID','"+str(expts[m].masterProjectID)+"'],"
		jsArray=jsArray+"['seqRunID','"+str(expts[m].seqRunID)+"'],"
		jsArray=jsArray+"['customerID','"+str(expts[m].customerID)+"'],"
		jsArray=jsArray+"['exptType','"+str(expts[m].exptType)+"']"
		jsArray=jsArray+"],["


		for n in range(0,len(expts[m].lanes)):
			jsArray=jsArray+"[["
			jsArray=jsArray+"['laneID','"+str(expts[m].lanes[n].laneID)+"'],"  
			jsArray=jsArray+"['laneNumber','"+str(expts[m].lanes[n].laneNumber)+"'],"
			jsArray=jsArray+"['sequencingConc','"+str(expts[m].lanes[n].sequencingConc)+"'],"
			jsArray=jsArray+"['read1ClusterDensity','"+str(expts[m].lanes[n].read1ClusterDensity)+"'],"
			jsArray=jsArray+"['PhiXspiked','"+str(expts[m].lanes[n].PhiXspiked)+"'],"
			jsArray=jsArray+"['spike','"+str(expts[m].lanes[n].spike)+"'],"
			jsArray=jsArray+"['spikeRatio','"+str(expts[m].lanes[n].spikeRatio)+"']"			
			jsArray=jsArray+"],["


			for o in range(0,len(expts[m].lanes[n].samples)):
				jsArray=jsArray+"["
				jsArray=jsArray+"['sampleID','"+str(expts[m].lanes[n].samples[o].sampleID)+"'],"
				jsArray=jsArray+"['sampleName','"+str(expts[m].lanes[n].samples[o].sampleName)+"'],"
				jsArray=jsArray+"['tagID','"+str(expts[m].lanes[n].samples[o].tagID)+"'],"
				jsArray=jsArray+"['laneID','"+str(expts[m].lanes[n].samples[o].laneID)+"'],"
				jsArray=jsArray+"['tagSequence','"+str(expts[m].lanes[n].samples[o].tagSequence)+"'],"
				jsArray=jsArray+"['tagKit','"+str(expts[m].lanes[n].samples[o].tagKit)+"'],"
				jsArray=jsArray+"['analysisID','"+str(expts[m].lanes[n].samples[o].analysisID)+"'],"
				jsArray=jsArray+"['adaptorSequence','"+str(expts[m].lanes[n].samples[o].adaptorSequence)+"']"
				jsArray=jsArray+"],"


			jsArray=jsArray[:-1]+"]],"	#Close the samples part and remove the ","
		jsArray=jsArray[:-1]+"]],"		#Close the lane part and remove the ","
	jsArray=jsArray+"];"			 	#Close the run part and remove the ","

	return(jsArray)





