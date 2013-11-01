#file getAllProjects.py 



def getAllProjects():
	import runQuery


	DBquery='select  masterProjectID, projectName, projectLead, status, description, openDate, lastUpdate from masterProject  order by masterProjectID desc'

	#print "<p>",DBquery,"</p>"

	res=runQuery.runQuery(DBquery)
	resList=[]
	for g in range(0,len(res)):
		tmpResList=[]

		for h in range(0,len(res[g])):
			tmpResList.append(res[g][h])
		DBquery="select childID,type from typeLinker where parentID="+str(res[g][0])
		chRes=runQuery.runQuery(DBquery)	
	
		for chPr in range(0,len(chRes)):

			if chRes[chPr][1]=="sequencing":
				DBquery="select seqProjectName from seqProject where seqProjectID="+str(chRes[chPr][0])			

				resC=runQuery.runQuery(DBquery)

				tmpResList.append(resC[0][0])	
			elif chRes[chPr][1]=="flow cytometry":
				DBquery="select flowCytProjectName from flowCytProject where flowCytProjectID="+str(chRes[chPr][0])                 
                                resC=runQuery.runQuery(DBquery)
                                tmpResList.append(resC[0][0])
			elif chRes[chPr][1]=="analysis":
				tmpResList.append(chRes[chPr])
			else:
				tmpResList.append(chRes[chPr])	


		resList.append(tmpResList)	
		
	#print res
	return(resList)





