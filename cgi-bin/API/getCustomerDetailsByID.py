#file getCustomerDetailsByID.py


def getCustomerDetailsByID(customerID):

	from runQuery import runQuery

	DBquery="select name,email,tel from customer where customerID="+str(customerID)

	res=runQuery(DBquery)
	return(res[0])






