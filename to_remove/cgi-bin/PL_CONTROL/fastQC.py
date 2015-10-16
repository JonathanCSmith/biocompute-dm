# file runFastQC.py

class fastQC:
    def __init__(self, i, j):

        import os
        import re

        self.i = i  # This is the seqProject
        self.j = j  # This is the run

        # establish the paths that we'll be using
        self.projectsRoot = "/var/www/biocis/link/demux"
        self.projectDir = str(i.seqProjectID) + "_" + re.sub('\s+', '_', i.seqProjectName)
        self.sampleDirs = os.path.join
        self.QCreportsDir = os.path.join(self.projectsRoot, self.projectDir, "QCreports")
        self.fastQCdir = os.path.join(self.QCreportsDir, "FastQC")

        # If they don't exist, make the relevant dirs

        self.makeDir(self.QCreportsDir)
        self.makeDir(self.fastQCdir)
        self.samples = []
        self.currentSample = "NULL"
        self.currentSampleNumber = 0
        self.lanes = []
        self.currentLane = "NULL"
        self.outDirNames = []
        self.currentOutDirName = "NULL"
        self.samplePaths = []
        self.currentSamplePath = "NULL"
        self.sampleNumbers = []
        self.JIDs = []
        self.currentJID = "NULL"

        self.demuxType = ""

    def establishPaths(self):
        import os
        import re
        # print "<p>Establishing Paths</p>"
        for la in self.i.lanes:
            self.currentLane = la

            sampleNumber = 1
            for sa in la.samples:

                self.currentSample = sa
                if self.demuxType == "casava":
                    self.currentSamplePath = os.path.join(self.projectsRoot, self.projectDir, 'processed',
                                                          'Project_' + re.sub('\s+', '_', self.i.seqProjectName),
                                                          'Sample_' + self.currentSample.sampleName)
                if self.demuxType == "bcl2fastq2":
                    self.currentSamplePath = os.path.join(self.projectsRoot, self.projectDir, 'processed')

                laneNumber = "%03g" % (self.currentLane.laneNumber)
                self.currentOutDirName = os.path.join(self.fastQCdir, self.currentSample.sampleName) + "_L" + laneNumber

                self.lanes.append(self.currentLane)
                self.samples.append(self.currentSample)
                self.samplePaths.append(self.currentSamplePath)
                self.outDirNames.append(self.currentOutDirName)
                self.sampleNumbers.append(sampleNumber)
                sampleNumber = sampleNumber + 1

    def setDemuxType(self):
        bits = self.j.dataLocation.split("_")
        machineName = bits[1]
        if machineName == "J00129":
            self.demuxType = "bcl2fastq2"
        else:
            self.demuxType = "casava"

    def getFastQCStatus(self):
        import runQuery
        import subprocess as sub
        import os
        # print "<p>Getting FastQC status.</p>"
        print
        "<table border='1'>"
        print
        "<tr><th>sample name</th><th>status</th><th>location</th><th>source location</th><th>JID</th><th>Read1</th><th>Read2</th></tr>"

        command = 'ssh biocis@apollo "source /etc/profile; qstat"'
        p = sub.Popen(command, shell=True, stdout=sub.PIPE, stdin=sub.PIPE).stdout
        qstatLine = p.readlines()

        for sa in range(0, len(self.samples)):
            DBquery = "select fastQCID,status,location,sourceLocation, JID from fastQC where sampleID=" + str(
                self.samples[sa].sampleID)
            res = runQuery.runQuery(DBquery)
            # print "<p>",DBquery,"</p>"
            # print "<p>",res,"</p>"
            laneNumber = "%03g" % (self.lanes[sa].laneNumber)
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

            fastQCID = res[arr][0]
            status = res[arr][1]
            location = res[arr][2]
            sourceLocation = res[arr][3]
            loc = location.split("/")
            loca = os.path.join("/link/demux/", loc[4], loc[5], loc[6], loc[7])

            read1StatsLocation = os.path.join(loca, self.samples[sa].sampleName + "_" + self.samples[
                sa].tagSequence + "_L" + laneNumber + "_R1_001_fastqc", "fastqc_report.html")
            read2StatsLocation = os.path.join(loca, self.samples[sa].sampleName + "_" + self.samples[
                sa].tagSequence + "_L" + laneNumber + "_R2_001_fastqc", "fastqc_report.html")

            if running == "Y":
                status = "running"
                print
                "<tr><td>" + self.samples[
                    sa].sampleName + "</td><td>" + status + "</td><td>" + location + "</td><td>" + sourceLocation + "</td><td>" + JID + "</td><td>--</td><td>--</td></tr>"
            else:
                status = "finished"

                print
                "<tr><td>" + self.samples[
                    sa].sampleName + "</td><td>" + status + "</td><td>" + location + "</td><td>" + sourceLocation + "</td><td>" + JID + "</td>"

                try:
                    if self.demuxType == "casava":
                        fastQChtmlFile1 = os.path.join(location, self.samples[sa].sampleName + "_" + self.samples[
                            sa].tagSequence + "_L" + laneNumber + "_R1_001_fastqc", "fastqc_report.html")
                    if self.demuxType == "bcl2fastq2":
                        fastQChtmlFile1 = os.path.join(location, self.samples[sa].sampleName + "_S" + str(
                            self.sampleNumbers[sa]) + "_L" + laneNumber + "_R1_001_fastqc", "fastqc_report.html")

                    os.stat(fastQChtmlFile1)
                    print
                    "<td><a href='" + read1StatsLocation + "'>read1</a></td>"
                except Exception:
                    try:
                        if self.demuxType == "casava":
                            read1StatsLocation = os.path.join(loca, self.samples[sa].sampleName + "_" + self.samples[
                                sa].tagSequence + "_L" + laneNumber + "_R1_001_fastqc.html")
                            fastQChtmlFile1 = os.path.join(location, self.samples[sa].sampleName + "_" + self.samples[
                                sa].tagSequence + "_L" + laneNumber + "_R1_001_fastqc.html")
                        if self.demuxType == "bcl2fastq2":
                            read1StatsLocation = os.path.join(loca, self.samples[sa].sampleName + "_S" + str(
                                self.sampleNumbers[sa]) + "_L" + laneNumber + "_R1_001_fastqc.html")
                            fastQChtmlFile1 = os.path.join(location, self.samples[sa].sampleName + "_S" + str(
                                self.sampleNumbers[sa]) + "_L" + laneNumber + "_R1_001_fastqc.html")
                        os.stat(fastQChtmlFile1)
                        print
                        "<td><a href='" + read1StatsLocation + "'>read1</a></td>"
                    except Exception:
                        print
                        "<td>_____</td>"

                try:
                    if self.demuxType == "casava":
                        fastQChtmlFile2 = os.path.join(location, self.samples[sa].sampleName + "_" + self.samples[
                            sa].tagSequence + "_L" + laneNumber + "_R2_001_fastqc", "fastqc_report.html")
                    if self.demuxType == "bcl2fastq2":
                        fastQChtmlFile2 = os.path.join(location, self.samples[sa].sampleName + "_S" + str(
                            self.sampleNumbers[sa]) + "_L" + laneNumber + "_R2_001_fastqc", "fastqc_report.html")
                        # fastQChtmlFile2=os.path.join(location,self.samples[sa].sampleName+"_S"+str(sa+1)+"_L"+laneNumber+"_R2_001_fastqc","fastqc_report.html")
                    os.stat(fastQChtmlFile2)
                    print
                    "<td><a href='" + read2StatsLocation + "'>read2</a></td>"
                except Exception:
                    try:
                        if self.demuxType == "casava":
                            read2StatsLocation = os.path.join(loca, self.samples[sa].sampleName + "_" + self.samples[
                                sa].tagSequence + "_L" + laneNumber + "_R2_001_fastqc.html")
                            fastQChtmlFile2 = os.path.join(location, self.samples[sa].sampleName + "_" + self.samples[
                                sa].tagSequence + "_L" + laneNumber + "_R2_001_fastqc.html")
                        if self.demuxType == "bcl2fastq2":
                            read2StatsLocation = os.path.join(loca, self.samples[sa].sampleName + "_S" + str(
                                self.sampleNumbers[sa]) + "_L" + laneNumber + "_R2_001_fastqc.html")
                            fastQChtmlFile2 = os.path.join(location, self.samples[sa].sampleName + "_S" + str(
                                self.sampleNumbers[sa]) + "_L" + laneNumber + "_R2_001_fastqc.html")

                        os.stat(fastQChtmlFile2)
                        print
                        "<td><a href='" + read2StatsLocation + "'>read2</a></td>"
                    except Exception:
                        print
                        "<td>_____</td>"

                print
                "</tr>"
            """
            if len(res)>0:
                print "<tr><td>"+self.samples[sa].sampleName+"</td><td>"+status+"</td><td>"+location+"</td><td>"+sourceLocation+"</td><td>"+JID+"</td></tr>"
            else:
                print "<tr><td>"+self.samples[sa].sampleName+"</td><td>QC not run</td><td>--</td><td>--</td><td>--</td><td>--</td><td>--</td></tr>"
            """

            # print "</tr>"
        print
        "</table>"

    def makeDir(self, dirToMake):
        import os

        try:
            if not os.path.exists(dirToMake):
                os.mkdir(dirToMake)
                os.chmod(dirToMake, 0777)
        except OSError:
            print
            "<p>Failed to make directory ", dirToMake, ". </p>"

    def goAllSamples(self):
        for la in self.i.lanes:
            self.currentLane = la
            sampleNumber = 1
            for sa in la.samples:
                self.currentSampleNumber = sampleNumber
                sampleNumber = sampleNumber + 1
                self.currentSample = sa
                # print "<p>",self.currentSample.sampleName,"</p>"
                self.makeComFile()
                self.queueFastQC()

    def makeComFile(self):
        import os
        import re

        if self.demuxType == "casava":
            self.currentSamplePath = os.path.join(self.projectsRoot, self.projectDir, 'processed',
                                                  'Project_' + re.sub('\s+', '_', self.i.seqProjectName),
                                                  'Sample_' + self.currentSample.sampleName)

            laneNumber = "%03g" % (self.currentLane.laneNumber)
            read1sampleName = self.currentSample.sampleName + "_" + self.currentSample.tagSequence + "_L" + laneNumber + "_R1_001.fastq.gz"
            read2sampleName = self.currentSample.sampleName + "_" + self.currentSample.tagSequence + "_L" + laneNumber + "_R2_001.fastq.gz"

            fileName = os.path.join(self.fastQCdir,
                                    "FastQC_" + self.currentSample.sampleName + "_L" + laneNumber + ".sh")

        if self.demuxType == "bcl2fastq2":
            self.currentSamplePath = os.path.join(self.projectsRoot, self.projectDir, 'processed')

            laneNumber = "%03g" % (self.currentLane.laneNumber)
            read1sampleName = self.currentSample.sampleName + "_S" + str(
                self.currentSampleNumber) + "_L" + laneNumber + "_R1_001.fastq.gz"
            read2sampleName = self.currentSample.sampleName + "_S" + str(
                self.currentSampleNumber) + "_L" + laneNumber + "_R2_001.fastq.gz"

            fileName = os.path.join(self.fastQCdir,
                                    "FastQC_" + self.currentSample.sampleName + "_L" + laneNumber + ".sh")

        OF = open(fileName, "w")

        header = "#!/bin/bash\n#$ -cwd\n#$ -j y\n#$ -S /bin/bash\n#$ -pe threaded 2\n#$ -N FastQC_" + self.currentSample.sampleName + "_L" + laneNumber + "\n\n"
        OF.write(header)

        self.currentOutDirName = os.path.join(self.fastQCdir, self.currentSample.sampleName) + "_L" + laneNumber
        self.makeDir(self.currentOutDirName)

        read1 = os.path.join(self.currentSamplePath, read1sampleName)
        read2 = os.path.join(self.currentSamplePath, read2sampleName)

        command = "/apps/FastQC/0.11.3/fastqc --threads 2 -o " + self.currentOutDirName + " " + read1 + " " + read2 + "\n\n"

        OF.write(command)
        OF.close()

    def queueFastQC(self):

        import subprocess as sub
        import re
        import os

        laneNumber = "%03g" % (self.currentLane.laneNumber)
        command = 'ssh biocis@athena "source /etc/profile; cd ' + self.fastQCdir + '; qsub -q shortterm.q,longterm.q,Admintest.q FastQC_' + self.currentSample.sampleName + '_L' + laneNumber + '.sh"'

        print
        "<p>" + command + "</p>"

        p = sub.Popen(command, shell=True, stdout=sub.PIPE, stdin=sub.PIPE).stdout
        line = p.readlines()
        print
        "<p>", line[0], "</p>"

        JIDfile = os.path.join(self.currentOutDirName, 'JID')

        jid = "0"
        if len(line[0]) > 0:
            bits = line[0].split()
            jid = bits[2]
            self.currentJID = jid
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

        # Once we have a JID, update the database

        # Check whether DBentry already exists for this sample and then either update the existing entry or create a new one.
        alreadyDone = self.checkDB()

        if alreadyDone == 1:
            self.updateDBentry()
        elif alreadyDone == 0:
            self.createDBentry()

    def checkDB(self):

        import runQuery

        DBquery = "select fastQCID from fastQC where sampleID=" + str(self.currentSample.sampleID)
        res = runQuery.runQuery(DBquery)

        return len(res)

    def updateDBentry(self):

        import runQuery

        DBquery = "update fastQC set location='" + self.currentOutDirName + "' , sourceLocation='" + self.currentSamplePath + "', JID=" + self.currentJID + ", status='running' where sampleID=" + str(
            self.currentSample.sampleID)

        # print "<p>DBquery "+DBquery+"</p>"
        res = runQuery.runQuery(DBquery)

    def createDBentry(self):

        import runQuery

        DBinsert = "insert into fastQC (seqProjectID,sampleID,location,sourceLocation,JID,status)"
        # print "<p>",str(self.i.seqProjectID),"</p>"

        DBvals = " values(" + str(self.i.seqProjectID) + "," + str(
            self.currentSample.sampleID) + ",'" + self.currentOutDirName + "','" + self.currentSamplePath + "'," + self.currentJID + ",'running')"

        DBquery = DBinsert + DBvals
        # print "<p>DBquery "+DBquery+"</p>"
        res = runQuery.runQuery(DBquery)
