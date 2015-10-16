# file seqProject.py

import laneData


class seqProject:
    def __init__(self, parent):
        self.parent = parent
        try:
            self.seqRunID = self.parent.seqRunID
        except AttributeError:
            self.seqRunID = "NULL"

        self.masterProjectID = 0
        self.seqProjectID = 0
        self.seqProjectName = "SeqProjName"
        self.customerID = 0
        self.exptType = "other"
        self.state = "RED"

        self.lanes = []

    def addLaneData(self):
        # add a laneData instance and pass it a refernce to the current object
        self.lanes.append(laneData.laneData(self))

    def getStateDB(self):
        import runQuery

        DBquery = "select state from state where itemID='" + str(self.seqProjectID) + "' and type='sequencing'"

        res = runQuery.runQuery(DBquery)
        self.state = res[0][0]

    def setStateDB(self):
        import runQuery

        DBquery = "update state set state='" + self.state + "'  where itemID='" + str(
            self.seqProjectID) + "' and type='sequencing'"
        # print "<p>",DBquery,"</p>"
        res = runQuery.runQuery(DBquery)

    def getLastSeqProjectID(self):
        # Return the number of the last lane in order to update the sampleData object just after data has been inserted
        import runQuery

        DBquery = "select seqProjectID from seqProject order by seqProjectID desc limit 1"

        seqProjectID = runQuery.runQuery(DBquery)

        self.seqProjectID = seqProjectID[0][0]

    def getSeqProjByID(self, projID):
        import runQuery

        DBquery = "select seqProjectName, masterProjectID, customerID, seqRunID, exptType from seqProject where seqProjectID=" + str(
            projID)

        res = runQuery.runQuery(DBquery)

        # if len(res[0])>2:
        self.seqProjectID = int(projID)
        self.seqProjectName = res[0][0]
        self.masterProjectID = res[0][1]
        self.customerID = res[0][2]
        self.seqRunID = res[0][3]
        if not res[0][4]:
            self.exptType = "NULL"
        else:
            self.exptType = res[0][4]

        self.getStateDB()

    def getLaneDataIDs(self):
        import runQuery

        DBquery = "select laneID from laneData where seqProjectID=" + str(self.seqProjectID) + " order by laneNumber"
        res = runQuery.runQuery(DBquery)

        return (res)

        # Deletes this sequence project and all lanes and samples associated with it

    def deleteSeqProjectTree(self):

        # first delete all child lanes
        laneDataIDs = self.getLaneDataIDs()
        for la in range(0, len(laneDataIDs)):
            # add a sample instance
            self.addLaneData()
            self.lanes[-1].getLaneDataByID(laneDataIDs[la][0])
            self.lanes[-1].deleteLaneTree()

            # Now remove this seqRun
        import runQuery

        DBquery = "delete from seqProject where seqProjectID=" + str(self.seqProjectID)
        res = runQuery.runQuery(DBquery)

        DBquery = "delete from state where itemID=" + str(self.seqProjectID) + " and type='sequencing'"
        res = runQuery.runQuery(DBquery)

        DBquery = "delete from typeLinker where childID=" + str(self.seqProjectID) + " and type='sequencing'"
        res = runQuery.runQuery(DBquery)

        DBquery = "delete from demultiplex where seqProjectID=" + str(self.seqProjectID)
        res = runQuery.runQuery(DBquery)

        DBquery = "delete from fastQC where seqProjectID=" + str(self.seqProjectID)
        res = runQuery.runQuery(DBquery)

        # insertSeqProj method

    def insertSeqProj(self):

        try:
            self.seqRunID = self.parent.seqRunID
        except AttributeError:
            self.seqRunID = "NULL"

        import runQuery

        insQuery = "INSERT INTO seqProject (seqProjectName,masterProjectID, seqRunID,customerID,exptType) "
        vals = " VALUES('" + self.seqProjectName + "','" + str(self.masterProjectID) + "','" + str(
            self.seqRunID) + "','" + str(self.customerID) + "','" + self.exptType + "')"

        DBins = insQuery + vals

        inser = runQuery.runQuery(DBins)

        self.getLastSeqProjectID()

        DBquery = "insert into state (state,itemID,type) values('" + self.state + "','" + str(
            self.seqProjectID) + "','sequencing')"
        # print "<p>",DBquery,"</p>"
        res = runQuery.runQuery(DBquery)

        # Update the database with the current object

    def updateDB(self):
        import runQuery

        updateQuery = "UPDATE seqProject SET seqRunID='" + str(
            self.seqRunID) + "', seqProjectName='" + self.seqProjectName + "', masterProjectID=" + str(
            self.masterProjectID) + ", customerID=" + str(
            self.customerID) + ", exptType='" + self.exptType + "' where seqProjectID=" + str(self.seqProjectID)
        # print "<p>",updateQuery,"</p>"
        update = runQuery.runQuery(updateQuery)

        # setSeqProjByDict method expects a dictionary called seqProjRecord as argument

    def setSeqProjByDict(self, seqProjRecord):

        if seqProjRecord.has_key('masterProjectID'):
            self.masterProjectID = seqProjRecord['masterProjectID']

        if seqProjRecord.has_key('seqProjectName'):
            self.seqProjectName = seqProjRecord['seqProjectName']

        if seqProjRecord.has_key('seqProjectID'):
            self.seqProjectID = seqProjRecord['seqProjectID']

        if seqProjRecord.has_key('customerID'):
            self.customerID = seqProjRecord['customerID']

        if seqProjRecord.has_key('exptType'):
            self.exptType = seqProjRecord['exptType']

        if seqProjRecord.has_key('seqRunID'):
            self.seqRunID = seqProjRecord['seqRunID']
