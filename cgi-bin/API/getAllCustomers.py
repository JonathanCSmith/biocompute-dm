#file getAllCustomers.py 



def getAllCustomers():
	import runQuery


	DBquery='select  customerID, name from customer order by customerID desc'
	custList=[]
	try:
		res=runQuery.runQuery(DBquery)

		for cu in range(0,len(res)):
			custList.append([res[cu][0],res[cu][1],0])
	except Exception:
		print "<p>DB lookup failed</p>"
		print "<p>",DBquery,"</p>"
	for cl in range(0,len(custList)):
			DBquery='select contactID,name, email,tel from customerContact where customerID='+str(custList[cl][0])
			try:
				res=runQuery.runQuery(DBquery)
				tmpCon=[]
				for con in range(0,len(res)):
					tmpCon.append([res[con][0],res[con][1],res[con][2],res[con][3]])

				custList[cl][2]=tmpCon[:]
			except Exception:
				print "<p>DB contact lookup failed</p>"
				print "<p>",DBquery,"</p>"




		
	return(custList)





