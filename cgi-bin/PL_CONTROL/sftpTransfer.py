# file sftpTransfer.py

class sftpTransfer:
    def __init__(self):
        import os

        self.exec_root = "/home/biocis/www/cgi-bin/PL_CONTROL"
        self.APIpath = "/var/www/biocis/cgi-bin/coreInSys/API"
        self.sftp_root = "/gpfs/home/transfer"

        self.seqRunName = ""
        self.username = ""
        self.password = ""
        self.userContact = ""
        self.creationDate = ""
        self.accountLocation = ""
        self.accountStatus = ""

        self.projectID = 0
        self.seqProjectID = 0
        self.transLocation = ""

        self.sftpAccountID = 0
        self.transferID = 0
        self.dataStatus = ""

    def setTransfer(self):
        import sys

        sys.path.append(self.APIpath)

        from runQuery import runQuery

        if self.transLocation == "":
            self.setTransLocation()

        DBquery = "select transferID,transLocation,sftpAccountID,dataStatus from transfer where seqProjectID=" + str(
            self.seqProjectID)

        # print "<p>",DBquery,"</p>"
        res = runQuery(DBquery)

        if len(res) > 0:
            self.transferID = res[0][0]
            self.transLocaton = res[0][1]
            self.sftpAccountID = res[0][2]
            if res[0][3]:
                self.dataStatus = res[0][3]
            else:
                self.dataStatus = ""
                # print res
                # print "<p>Transfer already exists in database.</p>"
            return

    # This should happen when the transfer happens
    # DBquery="insert into transfer(seqProjectID, transLocation, sftpAccountID) values("+str(self.seqProjectID)+",'"+self.transLocation+"',"+str(self.sftpAccountID)+")"
    # res=runQuery(DBquery)

    def setTransLocation(self):
        import sys

        import os
        # sys.exit()
        # os.environ['PATH']=os.environ['PATH']+':/var/www/cgi-bin/API'
        sys.path.append(self.APIpath)
        from seqProject import seqProject

        from getRunByID import getRunByID

        try:
            i = seqProject(0)
            i.getSeqProjByID(int(self.seqProjectID))

        except Exception:
            print
            "Unable to load seq project with ID " + str(self.seqProjectID)
            sys.exit()

        try:
            j = getRunByID(int(i.seqRunID))
        except Exception:
            print
            "Unable to get sequencing run" + str(i.seqRunID)
            sys.exit()

        self.seqRunName = j.seqRunName

        sftpPath = os.path.join(self.sftp_root, self.username, "data/")
        self.transLocation = os.path.join(sftpPath, j.seqRunName)

    def restore2demuxDir(self):
        import subprocess as sub
        import os
        import sys

        sys.path.append(self.APIpath)

        from seqProject import seqProject
        from runQuery import runQuery

        try:
            i = seqProject(0)
            i.getSeqProjByID(int(self.seqProjectID))

        except Exception:
            print
            "Unable to load seq project with ID " + str(self.seqProjectID)
            sys.exit()

        DBquery = "select transLocation, sftpAccountID from transfer where seqProjectID=" + str(self.seqProjectID)

        # print "<p>"+DBquery+"</p>"

        res = runQuery(DBquery)
        # print "<p>",res,"</p>"
        try:
            self.transLocation = res[0][0]
            self.sftpAccountID = res[0][1]
        except Exception:
            print
            "<p>Unable to find transfer details in database.</p>"
            sys.exit()

        DBquery = "select accountLocation from sftpAccount where sftpAccountID=" + str(self.sftpAccountID)
        # print "<p>"+DBquery+"</p>"
        res = runQuery(DBquery)
        try:
            self.accountLocation = res[0][0]
        except Exception:
            # print "<p>",res,"</p>"
            print
            "<p>Unable to find location of origin.</p>"
            sys.exit()
        sourcePATH = self.transLocation
        qcSourcePATH = sourcePATH + "/QCreports"
        dataSourcePATH = sourcePATH + "/processed"

        DBquery = "select location from demultiplex where seqProjectID=" + str(self.seqProjectID)
        # print "<p>"+DBquery+"</p>"

        res = runQuery(DBquery)
        print
        "<p>", res, "</p>"

        destPATH = "/gpfs" + res[0][0]
        print
        "<p>", destPATH, "</p>"

        dataDestPATH = destPATH + "/processed"
        qcDestPATH = destPATH + "/QCreports"

        sshComm = "ssh root@159.92.115.5 '"

        if os.path.islink(dataDestPATH):

            # linkCommand=sshComm+"rm "+dataDestPATH+"'"
            # print "<p>",linkCommand,"</p>"

            command = os.path.join(self.exec_root, "remLink") + " " + dataDestPATH
            p = sub.Popen(command, shell=True, stdout=sub.PIPE, stdin=sub.PIPE).stdout
            linkRemoved = p.readlines()
            print
            "<p>", command, "</p>"

            # command=sshComm+"mv "+dataSourcePATH+" "+dataDestPATH+"'"
            # print "<p>",command,"</p>"

            command = os.path.join(self.exec_root, "moveData") + " " + dataSourcePATH + " " + dataDestPATH
            p = sub.Popen(command, shell=True, stdout=sub.PIPE, stdin=sub.PIPE).stdout
            dataMoved = p.readlines()
            print
            "<p>", command, "</p>"





        else:
            print
            "<p>Data folder in demux is not a symbolic link. Aborting.</p>"





            # sys.exit()

        if os.path.islink(qcDestPATH):

            # linkCommand=sshComm+"rm "+qcDestPATH+"'"
            # print "<p>",linkCommand,"</p>"

            command = os.path.join(self.exec_root, "remLink") + " " + qcDestPATH
            p = sub.Popen(command, shell=True, stdout=sub.PIPE, stdin=sub.PIPE).stdout
            linkRemoved = p.readlines()
            print
            "<p>", command, "</p>"


            # command=sshComm+"mv "+qcSourcePATH+" "+qcDestPATH+"'"
            # print "<p>",command,"</p>"

            command = os.path.join(self.exec_root, "moveData") + " " + qcSourcePATH + " " + qcDestPATH
            p = sub.Popen(command, shell=True, stdout=sub.PIPE, stdin=sub.PIPE).stdout
            dataMoved = p.readlines()
            print
            "<p>", command, "</p>"



        else:
            print
            "<p>QC folder in demux is not a symbolic link. Aborting.</p>"

        DBquery = "update  transfer set dataStatus='origin' where seqProjectID=" + str(self.seqProjectID)
        res = runQuery(DBquery)

    # sys.exit()





    def tranfer2sftp(self):
        import subprocess as sub
        import os
        import sys

        sys.path.append(self.APIpath)

        from seqProject import seqProject
        from runQuery import runQuery
        from getRunByID import getRunByID

        try:
            i = seqProject(0)
            i.getSeqProjByID(int(self.seqProjectID))

        except Exception:
            print
            "Unable to load seq project with ID " + str(self.seqProjectID)
            sys.exit()

        DBquery = "select location from demultiplex where seqProjectID=" + str(self.seqProjectID)
        res = runQuery(DBquery)
        # print res[0][0]
        sourcePATH = "/gpfs" + res[0][0]
        qcSourcePATH = sourcePATH + "/QCreports"
        dataSourcePATH = sourcePATH + "/processed"

        # print sourcePATH
        # print i.seqProjectName
        # print i.seqRunID

        j = getRunByID(int(i.seqRunID))

        # print j.seqRunName

        # sftpPath=sftpPath+j.seqRunName

        DBquery = "select dataStatus from transfer where seqProjectID=" + str(self.seqProjectID)

        res = runQuery(DBquery)

        print
        "<p>res=", res, "</p>"

        if len(res) > 0:
            print
            "<p>res=", res[0][0], "</p>"
            if res[0][0] == 'sftp':
                print
                "<p>Data has already been moved to sftp area. Aborting transfer</p>"
                return
            else:
                DBquery = "update transfer set dataStatus='sftp' where seqProjectID=" + str(self.seqProjectID)
                print
                "<p>", DBquery, "</p>"
                res = runQuery(DBquery)


        else:
            DBquery = "insert into transfer(seqProjectID, transLocation, sftpAccountID,dataStatus) values(" + str(
                self.seqProjectID) + ",'" + self.transLocation + "'," + str(self.sftpAccountID) + ",'sftp')"
            res = runQuery(DBquery)

        sshComm = "ssh root@159.92.115.5 '"

        if os.path.islink(dataSourcePATH):
            print
            "<p>Data source directory is a symbolic link. Aborting!<br>"
            print
            "(This probably means the data have already been moved to the sftp area)</p>"
        else:

            # if 1:
            # command=sshComm+"mkdir "+self.transLocation+"'"
            # print "<p>",command,"</p>"

            command = os.path.join(self.exec_root, "makeDir") + " " + self.transLocation
            p = sub.Popen(command, shell=True, stdout=sub.PIPE, stdin=sub.PIPE).stdout
            dirMade = p.readlines()
            print
            "<p>", command, "</p>"


            # command=sshComm+"mv "+dataSourcePATH+" "+self.transLocation+"'"
            # print "<p>",command,"</p>"

            command = os.path.join(self.exec_root, "moveData") + " " + dataSourcePATH + " " + self.transLocation
            p = sub.Popen(command, shell=True, stdout=sub.PIPE, stdin=sub.PIPE).stdout
            linkRemoved = p.readlines()
            print
            "<p>", command, "</p>"

            # linkCommand=sshComm+"ln -s "+self.transLocation+"/processed "+dataSourcePATH+"'"
            # print "<p>",linkCommand,"</p>"


            command = os.path.join(self.exec_root,
                                   "makeLink") + " " + self.transLocation + "/processed " + dataSourcePATH
            p = sub.Popen(command, shell=True, stdout=sub.PIPE, stdin=sub.PIPE).stdout
            linkRemoved = p.readlines()
            print
            "<p>", command, "</p>"

        if os.path.islink(qcSourcePATH):
            print
            "<p>QC source directory is a symbolic link. Aborting!<br>"
            print
            "(This probably means the QC stats have already been moved to the sftp area)</p>"
        else:
            # if 1:
            # command=sshComm+"mv "+qcSourcePATH+" "+self.transLocation+"'"
            # print "<p>",command,"</p>"

            command = os.path.join(self.exec_root, "moveData") + " " + qcSourcePATH + " " + self.transLocation
            p = sub.Popen(command, shell=True, stdout=sub.PIPE, stdin=sub.PIPE).stdout
            dataMoved = p.readlines()
            print
            "<p>", command, "</p>"

            # linkCommand=sshComm+"ln -s "+self.transLocation+"/QCreports "+qcSourcePATH+"'"
            # print "<p>",linkCommand,"</p>"

            command = os.path.join(self.exec_root, "makeLink") + " " + self.transLocation + "/QCreports " + qcSourcePATH
            p = sub.Popen(command, shell=True, stdout=sub.PIPE, stdin=sub.PIPE).stdout
            linkRemoved = p.readlines()
            print
            "<p>", command, "</p>"

    # sys.exit()


    def getTransStat(self):
        import runQuery

        DBquery = "select transferID, transLocation, sftpAccountID from transfer where seqProjectID=" + str(
            self.seqProjectID)
        res = runQuery.runQuery(DBquery)
        if len(res) > 0:
            self.transferID = res[0][0]
            self.transLocation = res[0][1]
            self.sftpAccountID = res[0][2]

            DBquery = "select accountLocation, username, userContact, creationDate, accountStatus from sftpAccount where sftpAccountID=" + str(
                self.sftpAccountID)
            res = runQuery.runQuery(DBquery)
            if len(res) > 0:
                self.accountLocation = res[0][0]
                self.username = res[0][1]
                self.userContact = res[0][2]
                self.creationDate = res[0][3]
                self.accountStatus = res[0][4]

    def loadSftpByUsername(self):
        import runQuery

        DBquery = "select accountLocation, username, userContact, creationDate, accountStatus, sftpAccountID from sftpAccount where username='" + self.username + "'"
        res = runQuery.runQuery(DBquery)
        if len(res) > 0:
            self.accountLocation = res[0][0]
            self.username = res[0][1]
            self.userContact = res[0][2]
            self.creationDate = res[0][3]
            self.accountStatus = res[0][4]
            self.sftpAccountID = res[0][5]

    def genPassword(self):
        import string, random

        s = string.lowercase + string.digits + string.uppercase
        self.password = ''.join(random.sample(s, 8))

    def genUsername(self):
        if self.seqProjectID == 0:
            self.username = "NULL"
        else:
            self.username = "brc" + "%05i" % (int(self.seqProjectID))

    def resetPassword(self):
        self.genPassword()
        import subprocess as sub
        import os

        command = os.path.join(self.exec_root, "updatePassword") + " " + self.username + " " + self.password
        p = sub.Popen(command, shell=True, stdout=sub.PIPE, stdin=sub.PIPE).stdout
        passwordUpdated = p.readlines()

    def addSftpAccount(self):
        import runQuery
        import subprocess as sub
        import os

        import time

        self.creationDate = time.strftime('%Y-%m-%d')
        # if self.accountStatus!="restore":
        #	self.accountStatus="created"


        DBquery = "select creationDate, accountLocation, accountStatus,userContact from sftpAccount where username='" + self.username + "'"
        # print DBquery
        res = runQuery.runQuery(DBquery)
        if len(res) > 0:
            if res[0][2] != "restore":
                # print "<p>Account with username "+self.username+" already exists</p>"
                self.password = ""
                self.userContact = res[0][3]
                self.accountStatus = "existing"
                self.accountLocation = res[0][1]
                return
                # else:
            #	self.accountStatus="restore"


            # for now we want to generate these from elsewhere so commented
            # self.genUsername()
            # self.genPassword()

        if self.username == "NULL":
            print
            "<p>Need to generate a username before adding an account</p>"
            return 0

        if len(self.accountLocation) == 0:
            self.accountLocation = "/home/transfer/" + self.username
            # print command
        command = os.path.join(self.exec_root,
                               "runSftpSetup") + " " + self.username + " " + self.password + " " + self.accountLocation
        # print "<p>",command,"</p>"
        p = sub.Popen(command, shell=True, stdout=sub.PIPE, stdin=sub.PIPE).stdout
        accountMade = p.readlines()

        if self.accountStatus == "restore":

            DBquery = "update sftpAccount set accountStatus='created', creationDate='" + self.creationDate + "' where username='" + self.username + "'"
        else:
            self.accountStatus = 'created'
            DBquery = "insert into sftpAccount(username, creationDate, accountLocation, userContact, accountStatus) values('" + self.username + "','" + self.creationDate + "','" + self.accountLocation + "','" + self.userContact + "','" + self.accountStatus + "')"

        print
        DBquery
        res = runQuery.runQuery(DBquery)

        print
        "<p>Account created, password=" + self.password + "</p>"

    def delSftpAccount(self):

        import subprocess as sub
        import runQuery
        import os

        DBquery = "select sftpAccountID, creationDate, accountLocation, userContact, accountStatus from sftpAccount where username='" + self.username + "'"

        res = runQuery.runQuery(DBquery)
        print
        "<p>", res, "</p>"
        self.accountStatus = res[0][4]
        self.sftpAccountID = int(res[0][0])

        if self.accountStatus == 'created':
            command = os.path.join(self.exec_root, "runUserDel") + " " + self.username
            p = sub.Popen(command, shell=True, stdout=sub.PIPE, stdin=sub.PIPE).stdout
            accountDeleted = p.readlines()

        DBquery = "update sftpAccount set accountStatus='deleted' where sftpAccountID=" + str(self.sftpAccountID)
        res = runQuery.runQuery(DBquery)

    def restoreSftpAccount(self):

        import runQuery

        DBquery = "update sftpAccount set accountStatus='restore' where username='" + self.username + "'"
        res = runQuery.runQuery(DBquery)
        self.accountStatus = 'restore'

        # self.addSftpAccount()
