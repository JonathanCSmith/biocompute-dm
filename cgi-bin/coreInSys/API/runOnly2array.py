#file exptObj2array.py

import seqRun


# Convert run objects to an array and send a javascript var as a return.
def runOnly2array(expts):

	#print "<p>expts[0].seqRunName",expts[0].seqRunName,"</p>"

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
		jsArray=jsArray+"]]"
	jsArray=jsArray+"]"

	return(jsArray)





