# file demux.py


class demux:
    def __init__(self, i, j):

        import os, stat
        import re

        # i is a seqProject object
        # j is a seqExpr object

        self.i = i
        self.j = j
        # print "<p>In demux object instance</p>"
        # print "<p>"+j.dataLocation+"</p>"
        self.sourceLocation = "/gpfs/stagingTMP/" + j.dataLocation
        self.seqProjectID = i.seqProjectID
        self.JID = "0"
        self.status = "setup"
        self.projectsRoot = "/var/www/biocis/link/demux"
        # print "<p>"+i.seqProjectName+"___</p>"
        self.projectDir = str(i.seqProjectID) + "_" + re.sub('\s+', '_', i.seqProjectName)
        # print "<p>"+self.projectDir+"</p>"
        dirToMake = os.path.join(self.projectsRoot, self.projectDir)
        self.location = dirToMake
        # print"<p>location", self.location,"</p>"
        self.demuxType = "casava"  # This is currently either casava or bcl2fastq2
        try:
            os.mkdir(dirToMake)
            # os.chmod(dirToMake,stat.S_IWRITE)

            # os.chmod(dirToMake,stat.S_IRWXO)
        except os.error:
            print
            "<p>", os.path.join(self.projectsRoot, self.projectDir), "</p>"

        os.chmod(dirToMake, 511)
        self.params = {}

        self.params["--use-bases-mask"] = ""
        self.params["--input-dir"] = "/gpfs/stagingTMP/" + j.dataLocation + "/Data/Intensities/BaseCalls/"
        self.params["--intensities-dir"] = "/gpfs/stagingTMP/" + j.dataLocation + "/Data/Intensities/"
        # self.params["--positions-dir"]="/gpfs/stagingTMP/"+j.dataLocation+"/Data/Intensities/"
        self.params["--positions-format"] = ".clocs"
        self.params["--output-dir"] = os.path.join(self.projectsRoot, self.projectDir, "processed")
        self.params["--fastq-cluster-count"] = "500000000"
        self.params["--flowcell-id"] = j.flowcellID
        self.params["--tiles"] = ""
        self.params["--sample-sheet"] = os.path.join(self.projectsRoot, self.projectDir, "sampleSheet.csv")
        self.params["--mismatches"] = "0"
        self.options = {}
        self.options["--no-eamss"] = "Y"
        self.options["--force"] = "Y"
        self.options["--ignore-missing-stats"] = "Y"
        self.options["--ignore-missing-bcl"] = "Y"
        self.options["--ignore-missing-control"] = "Y"

        self.Bcl2FQparams = {}
        self.Bcl2FQparams["--runfolder-dir"] = "/gpfs/stagingTMP/" + j.dataLocation
        self.Bcl2FQparams["--use-bases-mask"] = ""
        self.Bcl2FQparams["--output-dir"] = os.path.join(self.projectsRoot, self.projectDir, "processed")
        self.Bcl2FQparams["--loading-threads"] = "2"
        self.Bcl2FQparams["--writing-threads"] = "2"
        self.Bcl2FQparams["--demultiplexing-threads"] = "1"
        self.Bcl2FQparams["--processing-threads"] = "7"
        self.Bcl2FQparams["--tiles"] = ""
        self.Bcl2FQparams["--sample-sheet"] = os.path.join(self.projectsRoot, self.projectDir, "sampleSheet.csv")
        self.Bcl2FQparams["--stats-dir"] = os.path.join(self.projectsRoot, self.projectDir, "processed", "Stats")
        self.Bcl2FQparams["--reports-dir"] = os.path.join(self.projectsRoot, self.projectDir, "processed", "Reports")
        self.Bcl2FQparams["--interop-dir"] = os.path.join(self.projectsRoot, self.projectDir, "processed", "Interop")
        self.Bcl2FQoptions = {}
        self.Bcl2FQoptions["--create-fastq-for-index-reads"] = "Y"
        self.Bcl2FQoptions["--ignore-missing-bcl"] = "Y"
        self.Bcl2FQoptions["--ignore-missing-positions"] = "Y"
        self.Bcl2FQoptions["--ignore-missing-filter"] = "Y"

    def setDemuxType(self):
        bits = self.j.dataLocation.split("_")
        machineName = bits[1]
        if machineName == "J00129":
            self.demuxType = "bcl2fastq2"
        else:
            self.demuxType = "casava"

    def insertDB(self):

        # need to make sure that the same seqproject isn't already there. If so just perform an update command.


        import runQuery

        insQuery = "INSERT INTO demultiplex (seqProjectID,status,sourceLocation,location,JID) "
        vals = " VALUES('" + str(self.seqProjectID) + "','" + str(self.status) + "','" + str(
            self.sourceLocation) + "','" + str(self.location) + "','" + str(self.JID) + "')"

        DBins = insQuery + vals
        inser = runQuery.runQuery(DBins)

    def updateDBStatus(self, demuxID, status):
        import runQuery

        DBquery = "update demultiplex set status='" + status + "' where demuxID=" + str(demuxID)
        # print "<p>",DBquery,"</p>"
        res = runQuery.runQuery(DBquery)

    def setTiles(self):
        # get the lanes from the seqRun object
        lanes = []
        for la in range(0, len(self.i.lanes)):
            lanes.append(self.i.lanes[la].laneNumber)
        if self.demuxType == "casava":
            self.params["--tiles"] = "s_["
            for g in range(0, len(lanes)):
                self.params["--tiles"] = self.params["--tiles"] + str(lanes[g]) + ","
            self.params["--tiles"] = self.params["--tiles"][:-1] + "]"
        if self.demuxType == "bcl2fastq2":
            self.Bcl2FQparams["--tiles"] = "s_" + str(lanes[0])
            for g in range(1, len(lanes)):
                self.Bcl2FQparams["--tiles"] = self.Bcl2FQparams["--tiles"] + ",s_" + str(lanes[g])

    def setUseBasesMask(self):
        # get the length of the barcodes
        BClen = len(self.i.lanes[0].samples[0].tagSequence)
        tagMask = "I" + str(BClen)
        N = int(self.j.indexTagCycles) - BClen
        for l in range(0, N):
            tagMask = tagMask + "N"
        if self.demuxType == "casava":
            self.params["--use-bases-mask"] = "Y" + str(self.j.readCycles) + "," + tagMask + "," + "Y" + str(
                self.j.readCycles)
        if self.demuxType == "bcl2fastq2":
            self.Bcl2FQparams["--use-bases-mask"] = "Y" + str(self.j.readCycles) + "," + tagMask + "," + "Y" + str(
                self.j.readCycles)

    def makeComFile(self, projectName, projectID):
        # write the command file
        commandFile = ""
        parList = self.params.keys()
        print
        '<form method="post" action="writeDemuxCommand?projName=' + projectName + '&projID=' + projectID + '">'
        print
        '<table border=0 >'
        for co in range(0, len(parList)):
            print
            "<tr><td>" + parList[co] + '</td><td><input type="text" name="' + parList[co] + '" value="' + self.params[
                parList[co]] + '" size="80"></td></tr>'
        optList = self.options.keys()
        for opts in range(0, len(optList)):
            print
            "<tr><td>" + optList[opts] + '</td><td><select name="' + optList[opts] + '" value="' + self.options[
                optList[opts]] + '" >'
            print
            '<option value="Y">Y</option>'
            print
            '<option value="N">N</option>'
            print
            "</select></td></tr>"

        print
        "</table>"
        print
        '<input type="submit" value="Submit">'
        print
        '</form>'

        for co in range(0, len(parList) - 1):
            commandFile = commandFile + parList[co] + "\t" + self.params[parList[co]] + " \\\n<br>"
        commandFile = commandFile + parList[-1] + "\t" + self.params[parList[-1]] + " \\\n<br>"

        # print "<p>",commandFile,"</p>"

    def makeComFileBcl2fastq2(self, projectName, projectID):
        # write the command file
        commandFile = ""
        parList = self.Bcl2FQparams.keys()

        print
        '<form method="post" action="writeDemuxCommand?projName=' + projectName + '&projID=' + projectID + '">'
        print
        '<table border=0 >'
        for co in range(0, len(parList)):
            print
            "<tr><td>" + parList[co] + '</td><td><input type="text" name="' + parList[co] + '" value="' + \
            self.Bcl2FQparams[parList[co]] + '" size="80"></td></tr>'
        optList = self.Bcl2FQoptions.keys()
        for opts in range(0, len(optList)):
            print
            "<tr><td>" + optList[opts] + '</td><td><select name="' + optList[opts] + '" value="' + self.Bcl2FQoptions[
                optList[opts]] + '" >'
            print
            '<option value="Y">Y</option>'
            print
            '<option value="N">N</option>'
            print
            "</select></td></tr>"

        print
        "</table>"
        print
        '<input type="submit" value="Submit">'
        print
        '</form>'

        for co in range(0, len(parList) - 1):
            commandFile = commandFile + parList[co] + "\t" + self.Bcl2FQparams[parList[co]] + " \\\n<br>"
        commandFile = commandFile + parList[-1] + "\t" + self.Bcl2FQparams[parList[-1]] + " \\\n<br>"

    # print "<p>",commandFile,"</p>"









    def queueDemuxJob(self):
        print
        "<p>In the queue method!</p>"
        import os, sys
        import subprocess as sub
        import re
        # print "<p>self.projectDir"+self.projectDir+"</p>"
        # casavaComFile=os.path.join(self.projectDir,"casava_Bcl2FastQ.sh")
        if self.demuxType == "casava":
            command = 'ssh biocis@athena "source /etc/profile; cd /var/www/biocis/link/demux/' + self.projectDir + '; qsub -q shortterm.q,longterm.q casava_Bcl2FastQ.sh"'
        if self.demuxType == "bcl2fastq2":
            command = 'ssh biocis@athena "source /etc/profile; cd /var/www/biocis/link/demux/' + self.projectDir + '; qsub -q shortterm.q,longterm.q bcl2fastq2.sh"'

        print
        "<p>" + command + "</p>"

        p = sub.Popen(command, shell=True, stdout=sub.PIPE, stdin=sub.PIPE).stdout
        line = p.readlines()

        print
        "<p>", line[0], "</p>"

        # sys.exit()
        JIDfile = '/var/www/biocis/link/demux/' + self.projectDir + '/JID'
        jid = "0"
        if len(line[0]) > 0:
            bits = line[0].split()
            jid = bits[2]
        print
        "Job ID=" + jid

        JIDF = open(JIDfile, "w")
        # This needs re-doing
        if re.match("[0-9]", jid):
            if re.match("[A-Z]|[a-z]", jid):
                JIDF.write("0")
            else:
                JIDF.write(jid)
                self.JID = jid
        else:
            JIDF.write("0")

        JIDF.close()

    def getDemuxStatus(self):
        import os
        import subprocess as sub
        # print "/home/biocis/demux/"+self.projectDir
        stuff = os.listdir("/var/www/biocis/link/demux/" + self.projectDir)

        # print "<p>",stuff,"</p>"
        DBquery = "select demuxID, status, location, sourceLocation, JID from demultiplex where seqProjectID='" + str(
            self.seqProjectID) + "'"
        # print "<p>"+DBquery+"</p>"

        import runQuery

        res = runQuery.runQuery(DBquery)

        # print res

        command = 'ssh biocis@apollo "source /etc/profile; qstat"'
        p = sub.Popen(command, shell=True, stdout=sub.PIPE, stdin=sub.PIPE).stdout
        qstatLine = p.readlines()

        print
        "<table border='1'>"
        print
        "<tr><th>demuxID</th><th>status</th><th>location</th><th>source location</th><th>JID</th></tr>"
        for arr in range(0, len(res)):
            runStat = ""
            running = "N"
            JID = str(res[arr][4])
            for lineOut in range(0, len(qstatLine)):
                m = qstatLine[lineOut].split()
                # print "<p>",m,"</p>"
                try:
                    jobID = int(m[0])
                except Exception:
                    jobID = 0
                if jobID == int(JID):
                    runStat = qstatLine[lineOut]
                    running = "Y"

            demuxID = res[arr][0]
            status = res[arr][1]
            location = res[arr][2]
            sourceLocation = res[arr][3]
            # self.updateDBStatus(demuxID,status)
            if running == "Y":
                status = "running"
            else:
                status = "finished"

            print
            "<tr><td>" + str(demuxID) + "</td><td>"
            if status == "finished":

                pa = location.split("/")
                # print pa
                projDir = pa[-1]
                if self.demuxType == "casava":
                    print
                    "<a href='/link/demux/" + projDir + "/processed/Basecall_Stats_" + self.j.flowcellID + "/Demultiplex_Stats.htm'>" + status + "</a>"

                if self.demuxType == "bcl2fastq2":
                    print
                    "<a href='/link/demux/" + projDir + "/processed/Reports/html/index.html'>" + status + "</a>"




            else:
                print
                status

            print
            "</td><td>" + location + "</td><td>" + sourceLocation + "</td><td>" + str(
                JID) + "</td><td>", runStat, "</td>"

            print
            "</tr>"
            self.updateDBStatus(demuxID, status)

        print
        "</table>"

    def writeCommandFile(self):
        import os
        import re

        if self.demuxType == "casava":
            fileName = os.path.join(self.projectsRoot, self.projectDir, "casava_Bcl2FastQ.sh")
        if self.demuxType == "bcl2fastq2":
            fileName = os.path.join(self.projectsRoot, self.projectDir, "bcl2fastq2.sh")

        try:
            OF = open(fileName, "w")
        except Exception, e:
            print
            "<p>", e, "</p>"

        header = "#!/bin/bash\n#$ -cwd\n#$ -j y\n#$ -S /bin/bash\n#$ -pe threaded 8\n\n"

        OF.write(header)
        # commands="/apps/casava/1.8.2/bin/configureBclToFastq.pl \\\n"
        if self.demuxType == "casava":
            commands = "/apps/casava/1.8.3/bin/configureBclToFastq.pl \\\n"
            parList = self.params.keys()
            for co in range(0, len(parList)):
                paramToWrite = re.sub('\s+', '_', self.params[parList[co]])
                commands = commands + "\t\t" + parList[co] + "\t" + paramToWrite + " \\\n"
                print
                "<p>", commands, "</p>"

                # commands=commands+"\t\t"+parList[-1]+"\t"+self.params[parList[-1]]+" \\\n"
                # OF.write(commands)
            optsText = ""
            optList = self.options.keys()
            # numbOptions=0
            for opts in range(0, len(optList)):
                if self.options[optList[opts]] == "Y":
                    optsText = optsText + "\t\t" + optList[opts] + " \\\n"
                    # numbOptions=numbOptions+1
            if len(optsText) > 0:
                optsText = optsText[:-2]
            else:
                commands = commands[:-2]
        if self.demuxType == "bcl2fastq2":
            commands = "/apps/bcl2fastq2/2.16.0.10/bin/bcl2fastq \\\n"
            parList = self.Bcl2FQparams.keys()
            for co in range(0, len(parList)):
                paramToWrite = re.sub('\s+', '_', self.Bcl2FQparams[parList[co]])
                commands = commands + "\t\t" + parList[co] + "\t" + paramToWrite + " \\\n"
                print
                "<p>", commands, "</p>"

                # commands=commands+"\t\t"+parList[-1]+"\t"+self.params[parList[-1]]+" \\\n"
                # OF.write(commands)
            optsText = ""
            optList = self.Bcl2FQoptions.keys()
            # numbOptions=0
            for opts in range(0, len(optList)):
                if self.Bcl2FQoptions[optList[opts]] == "Y":
                    optsText = optsText + "\t\t" + optList[opts] + " \\\n"
                    # numbOptions=numbOptions+1
            if len(optsText) > 0:
                optsText = optsText[:-2]
            else:
                commands = commands[:-2]

        OF.write(commands)
        OF.write(optsText)
        if self.demuxType == "casava":
            makeComm = "\n\nmake -j 8 -C " + self.params["--output-dir"] + "\n"
            OF.write(makeComm)

        OF.close()

    def makeSampleSheet(self):
        import re

        if self.demuxType == "casava":
            sampleSheet = "FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject\n"
        if self.demuxType == "bcl2fastq2":
            print
            "<p>", self.i.lanes[0].samples[0].tagSequence, "</p>"
            d = self.i.lanes[0].samples[0].tagSequence.split("-")
            # if len(d)>1:
            sampleSheet = "[data]\nFCID,Lane,SampleID,SampleRef,index,index2,Description,Control,Recipe,Operator,SampleProject\n"
            # else:
        #	 sampleSheet="[data]\nFCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject\n"

        # go through each lane
        FCID = self.j.flowcellID
        Control = "N"
        SampleRef = ""
        Description = ""
        Recipe = ""
        Operator = ""
        SampleProject = self.i.seqProjectName
        # go through each lane
        for la in range(0, len(self.i.lanes)):
            # Now each sample
            Lane = str(self.i.lanes[la].laneNumber)
            for sa in range(0, len(self.i.lanes[la].samples)):
                SampleID = self.i.lanes[la].samples[sa].sampleName
                Index = self.i.lanes[la].samples[sa].tagSequence

                if Index == "NULL":
                    Index = ""
                if self.demuxType == "casava":
                    sampleSheet = sampleSheet + re.sub('\s+', '_',
                                                       FCID + "," + Lane + "," + SampleID + "," + SampleRef + "," + Index + "," + Description + "," + Control + "," + Recipe + "," + Operator + "," + SampleProject) + "\n"

                if self.demuxType == "bcl2fastq2":
                    d = self.i.lanes[la].samples[sa].tagSequence.split("-")
                    if len(d) > 1:
                        Index = d[0]
                        Index2 = d[1]

                        sampleSheet = sampleSheet + re.sub('\s+', '_',
                                                           FCID + "," + Lane + "," + SampleID + "," + SampleRef + "," + Index + "," + Index2 + "," + Description + "," + Control + "," + Recipe + "," + Operator + "," + SampleProject) + "\n"
                    else:
                        sampleSheet = sampleSheet + re.sub('\s+', '_',
                                                           FCID + "," + Lane + "," + SampleID + "," + SampleRef + "," + Index + ",," + Description + "," + Control + "," + Recipe + "," + Operator + "," + SampleProject) + "\n"


                        # sampleSheet=re.sub('\s+', '_',sampleSheet)
                        # print "<p>"+sampleSheet+"</p>"
        SF = open(self.params["--sample-sheet"], "w")
        SF.write(sampleSheet)
        SF.close()


"""
FCID		,Lane	,SampleID	,SampleRef	,Index	,Description	,Control,Recipe	,Operator	,SampleProject
AD1TMBACXX	,6	,215A002	,		,CGATGT	,		,N	,	,		,
AD1TMBACXX,6,329A004,,TGACCA,,N,,,
AD1TMBACXX,6,156A007,,CAGATC,,N,,,



#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe threaded 8


/apps/casava/1.8.0/bin/configureBclToFastq.pl \
		--use-bases-mask Y100,I6N,Y100 \
		--no-eamss \
                --input-dir /gpfs/stagingTMP/130417_SN1014_0162_AD1TMBACXX/Data/Intensities/BaseCalls/ \
                --intensities-dir /gpfs/stagingTMP/130417_SN1014_0162_AD1TMBACXX/Data/Intensities/ \
                --positions-dir /gpfs/stagingTMP/130417_SN1014_0162_AD1TMBACXX/Data/Intensities/ \
                --positions-format _pos.txt \
		--output-dir /home/alan/runs/plagnol/processed_B \
                --force \
		--fastq-cluster-count 500000000 \
                --ignore-missing-stats \
		--ignore-missing-bcl \
                --ignore-missing-control \
		--flowcell-id AD1TMBACXX \
                --tiles s_[6] \
		--sample-sheet /home/alan/runs/plagnol/SampleSheet.csv \

make -j 8 -C /home/alan/runs/plagnol/processed_B
"""
