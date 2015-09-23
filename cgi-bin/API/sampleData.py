#file sampleData.py


class sampleData:
        def __init__(self,parent):
                self.parent=parent
		self.sampleID=0
		self.sampleName="NULL"
		self.tagID=0
		self.tagKit="NULL"
		self.tagSequence="NULL"
		self.analysisID=0
		self.adaptorSequence="NULL"
	

                try:
                        self.laneID=str(self.parent.laneID)      
                except AttributeError:
                        self.laneID="NULL"



        def getLastSampleID(self):
                #Return the number of the last sample
                import runQuery

                DBquery="select sampleID from sampleData order by sampleID desc limit 1"

                sampleID=runQuery.runQuery(DBquery)

                self.sampleID=int(sampleID[0][0])

	def getSampleDataByID(self,sampleID):
		import runQuery

		DBquery="select sampleName, tagID, laneID, tagSequence, analysisID, adaptorSequence, tagKit from sampleData where sampleID="+str(sampleID)

		res=runQuery.runQuery(DBquery)

		self.sampleID=int(sampleID)
                self.sampleName=res[0][0]
                self.tagID=res[0][1]
		self.laneID=res[0][2]
                self.tagSequence=res[0][3]
                self.analysisID=res[0][4]
                self.adaptorSequence=res[0][5]
		self.tagKit=res[0][6]


	def deleteSample(self):
		import runQuery

		DBquery="delete from sampleData where sampleID="+str(self.sampleID)

		res=runQuery.runQuery(DBquery)	
	
	#insertSample method
        def insertSample(self):


		try:
                        self.laneID=str(self.parent.laneID)
                except AttributeError:
                        self.laneID="NULL"



                import runQuery

		insQuery="INSERT INTO sampleData (sampleName, tagID, laneID, tagSequence, analysisID,adaptorSequence, tagKit)"
                vals=" VALUES( '"+self.sampleName+"',"+str(self.tagID)+","+str(self.laneID)+",'"+str(self.tagSequence)+"',"+str(self.analysisID)+",'"+self.adaptorSequence+"','"+self.tagKit+"')"

                DBins=insQuery+vals
		#print "<p>",DBins,"</p>"
                inser=runQuery.runQuery(DBins)

		#self.getLastSampleID()


	#setSampleByDict method expects a dictionary called sampleRecord as argument
	def setSampleByDict(self,sampleRecord):


		if sampleRecord.has_key('sampleName'):
			self.sampleName=sampleRecord['sampleName']

		if sampleRecord.has_key('sampleID'):
                        self.sampleID=sampleRecord['sampleID']

		if sampleRecord.has_key('laneID'):
                        self.laneID=str(sampleRecord['laneID'])

		if sampleRecord.has_key('tagID'):
			self.tagID=str(sampleRecord['tagID'])

		if sampleRecord.has_key('tagSequence'):
			self.tagSequence=sampleRecord['tagSequence']

		if sampleRecord.has_key('tagKit'):
                        self.tagKit=sampleRecord['tagKit']

		if sampleRecord.has_key('analysisID'):
                        self.analysisID=sampleRecord['analysisID']

		if sampleRecord.has_key('adaptorSequence'):
			self.adaptorSequence=sampleRecord['adaptorSequence']

	
#Update the database with the current object
        def updateDB(self):
                import runQuery

                updateQuery="UPDATE sampleData SET sampleName='"+str(self.sampleName)+"', laneID="+str(self.laneID)+", tagID="+str(self.tagID)+", tagSequence='"+str(self.tagSequence)+"', analysisID="+str(self.analysisID)+", adaptorSequence='"+str(self.adaptorSequence)+"', tagKit='"+str(self.tagKit)+"' where sampleID="+str(self.sampleID)
                #print updateQuery
                update=runQuery.runQuery(updateQuery)







