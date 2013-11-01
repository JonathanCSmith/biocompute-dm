#file getAllUnlinkedSeqProjects.py 



def getAllUnlinkedSeqProjects():
	import runQuery


	DBquery="select seqProjectID, seqProjectName, seqExptID, customerID from seqProject where masterProjectID=0 order by seqProjectID"
			#print "<p>",DBquery,"</p>"
	localRes=runQuery.runQuery(DBquery)

	return(localRes)





