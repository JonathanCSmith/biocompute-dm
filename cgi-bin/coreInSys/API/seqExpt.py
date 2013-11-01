#file seqExpt.py
#The sequencing experiment object

import seqProject

class seqExpt():
	def __init__(self):
		self.seqExptID="0"
		self.seqExptName="SeqExptName"
		self.flowcellID="NULL"
		self.startDate="0000-00-00"
		self.completionDate="0000-00-00"
		self.genomicsLead="NULL"
		self.dataLocation="NULL"
		self.indexTagCycles="0"
		self.readCycles="0"

		self.seqProjs=[]
		#self.laneData=[[]]
		#self.sampleData=[[[]]]

	def addSeqProject(self):
		#add a seqence project instance and pass it a refernce to the current object
		self.seqProjs.append(seqProject.seqProject(self))

	
	def getLastSeqExptID(self):
		#Return the number of the last lane in order to update the sampleData object just after data has been inserted
		import runQuery

		DBquery="select seqExptID from seqExperiment order by seqExptID desc limit 1"

		seqExptID=runQuery.runQuery(DBquery)
		#print seqExptID
		#print seqExptID[0][0]
		try: 
			self.seqExptID=seqExptID[0][0]
		except IndexError:
			self.seqExptID="0"




	def getSeqExptByID(self,exptID):
		import runQuery
		query='select flowcellID, startDate, completionDate, genomicsLead, dataLocation, indexTagCycles, readCycles,seqExptName from seqExperiment where seqExptID like "'+exptID+'"'

        	res=runQuery.runQuery(query)

        	self.seqExptID=exptID
        	self.flowcellID=res[0][0]
        	self.startDate=res[0][1]
        	self.completionDate=res[0][2]
        	self.genomicsLead=res[0][3]
        	self.dataLocation=res[0][4]
        	self.indexTagCycles=res[0][5]
        	self.readCycles=res[0][6]
		self.seqExptName=res[0][7]

		if not self.flowcellID:
			self.flowcellID="NULL"
		if not self.startDate:
			self.startDate="NULL"
		if not self.completionDate:
			self.completionDate="NULL"
		if not self.genomicsLead:
			self.genomicsLead="NULL"
		if  not self.dataLocation:
			self.dataLocation="NULL"
		if not self.indexTagCycles:
			self.indexTagCycles="NULL"
		if  not self.readCycles:
			self.readCycles="NULL"
		if not self.seqExptName:
			self.seqExptName="NULL"


	
	def getSeqProjects(self):
		import runQuery
		query='select seqProjectID from seqProject where seqExptID like "'+self.seqExptID+'"'
		res=runQuery.runQuery(query)
		
		return(res)



	#insertSeqExpt method
        def insertSeqExpt(self):

                import runQuery

		insQuery="INSERT INTO seqExperiment (seqExptName,flowcellID,startDate,completionDate,genomicsLead,dataLocation,indexTagCycles,readCycles) "
                vals=" VALUES ('"+self.seqExptName+"','"+str(self.flowcellID)+"','"+str(self.startDate)+"','"+str(self.completionDate)+"','"+self.genomicsLead+"','"+self.dataLocation+"',"+str(self.indexTagCycles)+","+str(self.readCycles)+")"

                DBins=insQuery+vals
		#print "<p>",DBins,"</p>"

                inser=runQuery.runQuery(DBins)

                #self.getLastSeqExptID()


	#update the seqExperiment
	def updateDB(self):
                import runQuery

                updateQuery="UPDATE seqExperiment SET seqExptName='"+self.seqExptName+"', flowcellID='"+str(self.flowcellID)+"', startDate='"+str(self.startDate)+"', completionDate='"+str(self.completionDate)+"', genomicsLead='"+str(self.genomicsLead)+"', dataLocation='"+str(self.dataLocation)+"', indexTagCycles="+str(self.indexTagCycles)+", readCycles="+str(self.readCycles)+" where seqExptID='"+str(self.seqExptID)+"'"

                #print updateQuery 
                update=runQuery.runQuery(updateQuery)



        #setSeqExptByDict method expects a dictionary called seqExptRecord as argument
        def setSeqExptByDict(self,seqExptRecord):


		if seqExptRecord.has_key('seqExptName'):
			self.seqExptName=seqExptRecord['seqExptName']

		if seqExptRecord.has_key('flowcellID'):
			self.flowcellID=seqExptRecord['flowcellID']

		if seqExptRecord.has_key('startDate'):
			self.startDate=seqExptRecord['startDate']

		if seqExptRecord.has_key('completionDate'):
			self.completionDate=seqExptRecord['completionDate']

		if seqExptRecord.has_key('genomicsLead'):
			self.genomicsLead=seqExptRecord['genomicsLead']

		if seqExptRecord.has_key('dataLocation'):
			self.dataLocation=seqExptRecord['dataLocation']

		if seqExptRecord.has_key('indexTagCycles'):
			self.indexTagCycles=seqExptRecord['indexTagCycles']

		if seqExptRecord.has_key('readCycles'):
			self.readCycles=seqExptRecord['readCycles']			

		if seqExptRecord.has_key('seqExptID'):
                        self.seqExptID=seqExptRecord['seqExptID']


		if self.seqExptID=="NULL" or self.seqExptID=="0":	
			self.getLastSeqExptID()



			

	

