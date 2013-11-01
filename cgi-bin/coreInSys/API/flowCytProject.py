#file flowCytProject.py


class flowCytProject:
	def __init__(self,parent):
		self.parent=parent
		try:
			self.flowCytExptID=self.parent.flowCytExptID
                except AttributeError:
			self.flowCytExptID="NULL"
	
		self.masterProjectID=0
		self.flowCytProjectID=0
		self.flowCytProjectName=""
		self.customerID=0


		self.lanes=[]



"""
	
	def getLastSeqProjectID(self):
		#Return the number of the last lane in order to update the sampleData object just after data has been inserted
		import runQuery

		DBquery="select seqProjectID from seqProject order by seqProjectID desc limit 1"

		seqProjectID=runQuery.runQuery(DBquery)

		self.seqProjectID=seqProjectID[0][0]


	def getSeqProjByID(self,projID):
		import runQuery

		DBquery="select seqProjectName, masterProjectID, customerID, seqExptID from seqProject where seqProjectID="+str(projID)
		
		res=runQuery.runQuery(DBquery)

		self.seqProjectID=int(projID)
                self.seqProjectName=res[0][0]
		self.masterProjectID=res[0][1]
                self.customerID=res[0][2]
		self.seqExptID=res[0][3]

	def getLaneDataIDs(self):
		import runQuery

		DBquery="select laneID from laneData where seqProjectID="+str(self.seqProjectID)
		res=runQuery.runQuery(DBquery)

                return(res)



        #insertSeqProj method
        def insertSeqProj(self):


		try:
                        self.seqExptID=self.parent.seqExptID
                except AttributeError:
                        self.seqExptID="NULL"



		import runQuery

		insQuery="INSERT INTO seqProject (seqProjectName,masterProjectID, seqExptID,customerID) "
                vals=" VALUES('"+self.seqProjectName+"','"+str(self.masterProjectID)+"','"+self.seqExptID+"',"+str(self.customerID)+")"

                DBins=insQuery+vals

                inser=runQuery.runQuery(DBins)
        
                self.getLastSeqProjectID()


	#Update the database with the current object
	def updateDB(self):
		import runQuery

		updateQuery="UPDATE seqProject SET seqProjectName='"+self.seqProjectName+"', masterProjectID="+str(self.masterProjectID)+", customerID="+str(self.customerID)+" where seqProjectID="+str(self.seqProjectID)
		#print updateQuery 
		update=runQuery.runQuery(updateQuery)






        #setSeqProjByDict method expects a dictionary called seqProjRecord as argument
        def setSeqProjByDict(self,seqProjRecord):


		if seqProjRecord.has_key('masterProjectID'):
                        self.masterProjectID=seqProjRecord['masterProjectID']

		if seqProjRecord.has_key('seqProjectName'):
			self.seqProjectName=seqProjRecord['seqProjectName']

		if seqProjRecord.has_key('seqProjectID'):
			self.seqProjectID=seqProjRecord['seqProjectID']

		if seqProjRecord.has_key('customerID'):
			self.customerID=seqProjRecord['customerID']

"""

	







