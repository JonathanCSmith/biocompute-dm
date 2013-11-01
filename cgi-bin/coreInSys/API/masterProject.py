#file masterProject.py
#The master project object

import seqProject
import getSeqProjByID

import flowCytProject
#import getFlowCytProjByID

class masterProject:
	def __init__(self):

		self.masterProjectID="NULL"
		self.status="analysis"
		self.projectName=""
		self.projectLead=""
		self.description=""
		self.openDate=""
		self.lastUpdate=""
		self.children=[]
		self.childTypes=[]

		#self.seqProjs=[]
		#self.flowCytProjs=[]




	def addChild(self,childType,projectID):

		if childType=='sequencing':
			self.children.append(getSeqProjByID.getSeqProjByID(self,projectID))
                	self.childTypes.append(childType)
		elif childType=='flow cytometry':
			#self.children.append(getFlowCytProjByID.getflowCytProjByID(projectID))
			self.children.append([])
                	self.childTypes.append(childType)
		elif childType=='analysis':
			self.children.append([])
                        self.childTypes.append(childType)
		else:
			self.children.append([])
			self.childTypes.append(childType)	


	def setLastUpdate(self):
		import datetime
		today=datetime.date.today()
		self.lastUpdate=str(today.year)+"-"+str(today.month)+"-"+str(today.day)



	
	def getLastMasterProjectID(self):
		import runQuery
		DBquery="select masterProjectID from masterProject order by masterProjectID desc limit 1"

		masterProjectID=runQuery.runQuery(DBquery)
		self.masterProjectID=masterProjectID[0][0]
		#return self.masterProjectID



	def getMasterProjectByID(self,masterProjectID):
		import runQuery
		query='select status, projectName, projectLead, description,openDate, lastUpdate from masterProject where masterProjectID='+str(masterProjectID)



        	res=runQuery.runQuery(query)
		self.masterProjectID=masterProjectID
                self.status=res[0][0]
                self.projectName=res[0][1]
                self.projectLead=res[0][2]
		self.description=res[0][3]
		self.openDate=res[0][4]
		self.lastUpdate=res[0][5]

	def appendChildren(self):
		import runQuery

		query='select type, childID from typeLinker where parentID='+self.masterProjectID
		#print "<p>",query,"</p>"
		res=runQuery.runQuery(query)
		for li in range(0,len(res)):
			
			if res[li][0]=="sequencing":
				self.addChild(res[li][0],res[li][1])     
			elif res[li][0]=="flow cytometry":
				self.addChild(res[li][0],res[li][1])
				#self.children[-1].getFlowCytProjByID(res[li][1])	
			else:
				self.addChild(res[li][0],res[li][1])



	
	def getSeqProjects(self):
		import runQuery
		query='select seqProjectID from seqProject where masterProjectID like "'+self.masterProjectID+'"'
		res=runQuery.runQuery(query)
		
		return(res)

	def getFLowCytProjects(self):
                import runQuery
                query='select flowCytProjectID from flowCytProject where masterProjectID like "'+self.masterProjectID+'"'
                res=runQuery.runQuery(query)

                return(res)




	def insertMasterProject(self):
		import runQuery

		insQuery="INSERT INTO masterProject (projectName,projectLead,status,description, openDate,lastUpdate) "
		vals=" VALUES ('"+self.projectName+"','"+self.projectLead+"','"+str(self.status)+"','"+str(self.description)+"','"+self.openDate+"','"+self.lastUpdate+"')"

		DBins=insQuery+vals
		#print "<p>",DBins,"</p>"
                inser=runQuery.runQuery(DBins)

                self.getLastMasterProjectID()


	#Update the database with the current object
	def updateDB(self):
		import runQuery

		updateQuery="UPDATE masterProject SET projectName='"+self.projectName+"', projectLead='"+self.projectLead+"', status='"+str(self.status)+"' ,description='"+self.description+"', openDate='"+self.openDate+"', lastUpdate='"+self.lastUpdate+"' where masterProjectID="+str(self.masterProjectID)
		#print updateQuery 
		update=runQuery.runQuery(updateQuery)




        #setMasterProjByDict method expects a dictionary called seqProjRecord as argument
        def setMasterProjByDict(self,masterProjRecord):

	        if masterProjRecord.has_key('masterProjectID'):
			self.masterProjectID=masterProjRecord['masterProjectID']

                if masterProjRecord.has_key('projectName'):
                        self.projectName=masterProjRecord['projectName']

		if masterProjRecord.has_key('projectLead'):
                        self.projectLead=masterProjRecord['projectLead']

		if masterProjRecord.has_key('status'):
                        self.status=masterProjRecord['status']

		if masterProjRecord.has_key('description'):
                        self.description=masterProjRecord['description']

		if masterProjRecord.has_key('openDate'):
                        self.openDate=masterProjRecord['openDate']

		if masterProjRecord.has_key('lastUpdate'):
                        self.lastUpdate=masterProjRecord['lastUpdate']







			

	

