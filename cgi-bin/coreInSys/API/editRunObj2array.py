#file exptObj2array.py

import seqRun


# Convert run objects to an array and send a javascript var as a return.
def runObj2array (expts):
	jsArray="var runArr=["
	
	for l in range(0,len(expts)):
		jsArray=jsArray+"[["
		jsArray=jsArray+"['seqRunID','"+str(expts[l].seqRunID)+"'],"
		jsArray=jsArray+"['seqRunName','"+str(expts[l].seqRunName)+"'],"
		jsArray=jsArray+"['flowcellID','"+str(expts[l].flowcellID)+"'],"
		jsArray=jsArray+"['startDate','"+str(expts[l].startDate)+"'],"
		jsArray=jsArray+"['completionDate','"+str(expts[l].completionDate)+"'],"
		jsArray=jsArray+"['genomicsLead','"+str(expts[l].genomicsLead)+"'],"
		jsArray=jsArray+"['dataLocation','"+str(expts[l].dataLocation)+"'],"
		jsArray=jsArray+"['indexTagCycles','"+str(expts[l].indexTagCycles)+"'],"
		jsArray=jsArray+"['readCycles','"+str(expts[l].readCycles)+"']"
		jsArray=jsArray+"],["

		for m in range(0,len(expts[l].seqProjs)):
			jsArray=jsArray+"[["
			jsArray=jsArray+"['seqProjectID','"+str(expts[l].seqProjs[m].seqProjectID)+"'],"
			jsArray=jsArray+"['seqProjectName','"+str(expts[l].seqProjs[m].seqProjectName)+"'],"
			jsArray=jsArray+"['masterProjectID','"+str(expts[l].seqProjs[m].masterProjectID)+"'],"
			jsArray=jsArray+"['seqRunID','"+str(expts[l].seqProjs[m].seqRunID)+"'],"
			jsArray=jsArray+"['customerID','"+str(expts[l].seqProjs[m].customerID)+"'],"
			jsArray=jsArray+"['exptType','"+str(expts[l].seqProjs[m].exptType)+"']"
			jsArray=jsArray+"],["


			for n in range(0,len(expts[l].seqProjs[m].lanes)):
				jsArray=jsArray+"[["
				jsArray=jsArray+"['laneID','"+str(expts[l].seqProjs[m].lanes[n].laneID)+"'],"  
				jsArray=jsArray+"['laneNumber','"+str(expts[l].seqProjs[m].lanes[n].laneNumber)+"'],"
				jsArray=jsArray+"['seqProjectID','"+str(expts[l].seqProjs[m].lanes[n].seqProjectID)+"'],"
				jsArray=jsArray+"['sequencingConc','"+str(expts[l].seqProjs[m].lanes[n].sequencingConc)+"'],"
				jsArray=jsArray+"['read1ClusterDensity','"+str(expts[l].seqProjs[m].lanes[n].read1ClusterDensity)+"'],"
				jsArray=jsArray+"['PhiXspiked','"+str(expts[l].seqProjs[m].lanes[n].PhiXspiked)+"'],"
				jsArray=jsArray+"['spike','"+str(expts[l].seqProjs[m].lanes[n].spike)+"'],"
				jsArray=jsArray+"['spikeRatio','"+str(expts[l].seqProjs[m].lanes[n].spikeRatio)+"']"			
				jsArray=jsArray+"],["


				for o in range(0,len(expts[l].seqProjs[m].lanes[n].samples)):
					jsArray=jsArray+"["
					jsArray=jsArray+"['sampleID','"+str(expts[l].seqProjs[m].lanes[n].samples[o].sampleID)+"'],"
					jsArray=jsArray+"['sampleName','"+str(expts[l].seqProjs[m].lanes[n].samples[o].sampleName)+"'],"
					jsArray=jsArray+"['tagID','"+str(expts[l].seqProjs[m].lanes[n].samples[o].tagID)+"'],"
					jsArray=jsArray+"['laneID','"+str(expts[l].seqProjs[m].lanes[n].samples[o].laneID)+"'],"
					jsArray=jsArray+"['tagSequence','"+str(expts[l].seqProjs[m].lanes[n].samples[o].tagSequence)+"'],"
					jsArray=jsArray+"['analysisID','"+str(expts[l].seqProjs[m].lanes[n].samples[o].analysisID)+"'],"
					jsArray=jsArray+"['adaptorSequence','"+str(expts[l].seqProjs[m].lanes[n].samples[o].adaptorSequence)+"']"
					jsArray=jsArray+"],"


				jsArray=jsArray[:-1]+"]],"	#Close the samples part and remove the ","
			jsArray=jsArray[:-1]+"]],"		#Close the lane part and remove the ","
		jsArray=jsArray[:-1]+"]]"			#Close the project part and remove the ","
	jsArray=jsArray+"];"			 	#Close the run part and remove the ","
	return(jsArray)





