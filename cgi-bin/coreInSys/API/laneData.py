#file laneData.py

import sampleData

class laneData:
	def __init__(self,parent):
		self.parent=parent
		self.laneID="NULL"  
		self.laneNumber="NULL"
		try:
			self.seqProjectID=self.parent.seqProjectID
		except AttributeError:
			self.seqProjectID="NULL"

		self.sequencingConc="NULL"
		self.read1ClusterDensity=0
		self.PhiXspiked="NULL"
		self.spike="NULL"
		self.spikeRatio=0.0			

		self.samples=[]

	def addSamples(self):
		#add a sampleData instance and pass it a refernce to the current object
		self.samples.append(sampleData.sampleData(self))

	def getLastLaneID(self):
		#Return the number of the last lane in order to update the sampleData object just after data has been inserted
		import runQuery

		DBquery="select laneID from laneData order by laneID desc limit 1"

		laneID=runQuery.runQuery(DBquery)

		self.laneID=int(laneID[0][0])


	def getSampleIDs(self):
		import runQuery
		
		DBquery="select sampleID from sampleData where laneID="+str(self.laneID)

		res=runQuery.runQuery(DBquery)

		return(res)

	#Update the database with the current object
	def updateDB(self):
		import runQuery

		updateQuery="UPDATE laneData SET laneNumber='"+str(self.laneNumber)+"', seqProjectID="+str(self.seqProjectID)+", sequencingConc='"+str(self.sequencingConc)+"', read1ClusterDensity="+str(self.read1ClusterDensity)+", PhiXspiked='"+str(self.PhiXspiked)+"', spike='"+str(self.spike)+"', spikeRatio="+str(self.spikeRatio)+" where laneID="+str(self.laneID)
		#print updateQuery 
		update=runQuery.runQuery(updateQuery)



	def getLaneDataByID(self,laneID):
		import runQuery

		DBquery="select laneNumber, seqProjectID, sequencingConc, read1ClusterDensity, PhiXspiked, spike, spikeRatio from laneData where laneID="+str(laneID)

		res=runQuery.runQuery(DBquery)

		self.laneID=int(laneID)
                self.laneNumber=res[0][0]
                self.seqProjectID=res[0][1]
                self.sequencingConc=res[0][2]
                self.read1ClusterDensity=res[0][3]
                self.PhiXspiked=res[0][4]
                self.spike=res[0][5]
                self.spikeRatio=res[0][6]



	#Deletes this lane and all samples associated with it
        def deleteLaneTree(self):

		#first delete all child samples
		sampleIDs=self.getSampleIDs()
		for sa in range(0,len(sampleIDs)):
			#add a sample instance
                        self.addSamples()
                        self.samples[-1].getSampleDataByID(sampleIDs[sa][0])	
			self.samples[-1].deleteSample()	

		#Now remove this lane
                import runQuery
                DBquery="delete from laneData where laneID="+str(self.laneID)
                res=runQuery.runQuery(DBquery)

		



        #insertLane method 
        def insertLane(self):
		try:
                        self.seqProjectID=self.parent.seqProjectID
                except AttributeError:
                        self.seqProjectID="NULL"

	
                import runQuery

                insQuery="INSERT INTO laneData (seqProjectID,laneNumber,sequencingConc,read1ClusterDensity,PhiXspiked,spike,spikeRatio) "
                vals=" VALUES('"+str(self.seqProjectID)+"',"+str(self.laneNumber)+","+str(self.sequencingConc)+",'"+str(self.read1ClusterDensity)+"',"+str(self.PhiXspiked)+",'"+str(self.spike)+"',"+str(self.spikeRatio)+")"

                DBins=insQuery+vals
                inser=runQuery.runQuery(DBins)

                self.getLastLaneID()


        #setLane method expects a dictionary called laneRecord as argument
        def setLaneByDict(self,laneRecord):

		if laneRecord.has_key('seqProjectID'):
			self.seqProjectID=laneRecord['seqProjectID']
		
		if laneRecord.has_key('laneID'):
                        self.laneID=laneRecord['laneID']

		if laneRecord.has_key('laneNumber'):
                        self.laneNumber=laneRecord['laneNumber']

		if laneRecord.has_key('sequencingConc'):
			self.sequencingConc=laneRecord['sequencingConc']

		if laneRecord.has_key('read1ClusterDensity'):
                        self.read1ClusterDensity=laneRecord['read1ClusterDensity']

		if laneRecord.has_key('PhiXspiked'):
                        self.PhiXspiked=laneRecord['PhiXspiked']

		if laneRecord.has_key('spike'):
			self.spike=laneRecord['spike']

		if laneRecord.has_key('spikeRatio'):
			self.spikeRatio=laneRecord['spikeRatio']


	




