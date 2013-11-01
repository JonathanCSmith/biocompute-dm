

def getState(itemType,itemID):
	from runQuery import runQuery

	DBquery="select state from state where type='"+itemType+"' and itemID='"+str(itemID)+"'"

	res=runQuery(DBquery)
	state=res[0][0]
	return state





