#file exptObj2array.py

import seqExpt


# Convert experiment objects to an array and send a javascript var as a return.
def exptOnly2array(expts):

	#print "<p>expts[0].seqExptName",expts[0].seqExptName,"</p>"

	jsArray="var experArr=["
	for l in range(0,len(expts)):
		jsArray=jsArray+"[["
		jsArray=jsArray+"['seqExptID','"+str(expts[l].seqExptID)+"'],"
		jsArray=jsArray+"['seqExptName','"+str(expts[l].seqExptName)+"'],"
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





