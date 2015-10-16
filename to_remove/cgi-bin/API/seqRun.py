# file seqRun.py
# The sequencing run object

import seqProject


class seqRun():
    def __init__(self):
        self.seqRunID = 0
        self.seqRunName = "SeqRunName"
        self.flowcellID = "NULL"
        self.startDate = "0000-00-00"
        self.completionDate = "0000-00-00"
        self.genomicsLead = "NULL"
        self.dataLocation = "NULL"
        self.indexTagCycles = "0"
        self.readCycles = "0"

        self.seqProjs = []

        # self.laneData=[[]]
        # self.sampleData=[[[]]]

    def addSeqProject(self):
        # add a seqence project instance and pass it a refernce to the current object
        self.seqProjs.append(seqProject.seqProject(self))

    def getLastSeqRunID(self):
        # Return the number of the last lane in order to update the sampleData object just after data has been inserted
        import runQuery

        DBquery = "select seqRunID from seqRun order by seqRunID desc limit 1"

        seqRunID = runQuery.runQuery(DBquery)
        # print seqRunID
        # print seqRunID[0][0]
        try:
            self.seqRunID = seqRunID[0][0]
        except IndexError:
            self.seqRunID = "0"

    def getSeqNameByID(self, exptID):
        import runQuery
        # query='select flowcellID, startDate, completionDate, genomicsLead, dataLocation, indexTagCycles, readCycles,seqRunName from seqRun where seqRunID like "'+exptID+'"'
        query = 'select flowcellID, startDate, completionDate, genomicsLead, dataLocation, indexTagCycles, readCycles,seqRunName from seqRun where seqRunID=' + str(
            exptID)

        res = runQuery.runQuery(query)

        self.seqRunID = exptID
        self.flowcellID = res[0][0]
        self.startDate = res[0][1]
        self.completionDate = res[0][2]
        self.genomicsLead = res[0][3]
        self.dataLocation = res[0][4]
        self.indexTagCycles = res[0][5]
        self.readCycles = res[0][6]
        self.seqRunName = res[0][7]

        if not self.flowcellID:
            self.flowcellID = "NULL"
        if not self.startDate:
            self.startDate = "NULL"
        if not self.completionDate:
            self.completionDate = "NULL"
        if not self.genomicsLead:
            self.genomicsLead = "NULL"
        if not self.dataLocation:
            self.dataLocation = "NULL"
        if not self.indexTagCycles:
            self.indexTagCycles = "NULL"
        if not self.readCycles:
            self.readCycles = "NULL"
        if not self.seqRunName:
            self.seqRunName = "NULL"

    def getSeqProjects(self):
        import runQuery
        # query='select seqProjectID from seqProject where seqRunID like "'+self.seqRunID+'"'
        query = 'select seqProjectID from seqProject where seqRunID=' + str(self.seqRunID)
        res = runQuery.runQuery(query)

        return (res)

        # insertSeqRun method

    def insertSeqRun(self):

        import runQuery

        insQuery = "INSERT INTO seqRun (seqRunName,flowcellID,startDate,completionDate,genomicsLead,dataLocation,indexTagCycles,readCycles) "
        vals = " VALUES ('" + self.seqRunName + "','" + str(self.flowcellID) + "','" + str(
            self.startDate) + "','" + str(
            self.completionDate) + "','" + self.genomicsLead + "','" + self.dataLocation + "'," + str(
            self.indexTagCycles) + "," + str(self.readCycles) + ")"

        DBins = insQuery + vals
        # print "<p>",DBins,"</p>"

        inser = runQuery.runQuery(DBins)

        # self.getLastSeqRunID()

        # update the seqRun

    def updateDB(self):
        import runQuery

        updateQuery = "UPDATE seqRun SET seqRunName='" + self.seqRunName + "', flowcellID='" + str(
            self.flowcellID) + "', startDate='" + str(self.startDate) + "', completionDate='" + str(
            self.completionDate) + "', genomicsLead='" + str(self.genomicsLead) + "', dataLocation='" + str(
            self.dataLocation) + "', indexTagCycles=" + str(self.indexTagCycles) + ", readCycles=" + str(
            self.readCycles) + " where seqRunID=" + str(self.seqRunID)

        # print updateQuery
        update = runQuery.runQuery(updateQuery)

        # setSeqRunByDict method expects a dictionary called seqRunRecord as argument

    def setSeqRunByDict(self, seqRunRecord):

        if seqRunRecord.has_key('seqRunName'):
            self.seqRunName = seqRunRecord['seqRunName']

        if seqRunRecord.has_key('flowcellID'):
            self.flowcellID = seqRunRecord['flowcellID']

        if seqRunRecord.has_key('startDate'):
            self.startDate = seqRunRecord['startDate']

        if seqRunRecord.has_key('completionDate'):
            self.completionDate = seqRunRecord['completionDate']

        if seqRunRecord.has_key('genomicsLead'):
            self.genomicsLead = seqRunRecord['genomicsLead']

        if seqRunRecord.has_key('dataLocation'):
            self.dataLocation = seqRunRecord['dataLocation']

        if seqRunRecord.has_key('indexTagCycles'):
            self.indexTagCycles = seqRunRecord['indexTagCycles']

        if seqRunRecord.has_key('readCycles'):
            self.readCycles = seqRunRecord['readCycles']

        if seqRunRecord.has_key('seqRunID'):
            self.seqRunID = int(seqRunRecord['seqRunID'])

        if self.seqRunID == "NULL" or self.seqRunID == 0:
            self.getLastSeqRunID()
