# file postAlignQC.py

class postAlignQC:
    def __init__(self, i, j):

        import os
        import re

        self.i = i  # This is the seqProject
        self.j = j  # This is the run

        # establish the paths that we'll be using
        self.projectsRoot = "/home/biocis/demux"
        self.projectDir = str(i.seqProjectID) + "_" + re.sub('\s+', '_', i.seqProjectName)
        self.sampleDirs = os.path.join
        self.QCreportsDir = os.path.join(self.projectsRoot, self.projectDir, "QCreports")
        self.postAlignQCdir = os.path.join(self.QCreportsDir, "postAlignQC")



        # If they don't exist, make the relevant dirs

        self.makeDir(self.QCreportsDir)
        self.makeDir(self.postAlignQCdir)

        self.samples = []
        self.currentSample = "NULL"
        self.lanes = []
        self.currentLane = "NULL"
        self.outDirNames = []
        self.currentOutDirName = "NULL"
        self.samplePaths = []
        self.currentSamplePath = "NULL"
        self.JIDs = []
        self.currentJID = "NULL"

        # files that can be uploaded via biocis
        self.PAQfiles = []

    def establishPaths(self):
        import os
        import re
        # print "<p>Establishing Paths</p>"
        for la in self.i.lanes:
            self.currentLane = la

            for sa in la.samples:
                self.currentSample = sa
                self.currentSamplePath = os.path.join(self.projectsRoot, self.projectDir, 'processed',
                                                      'Project_' + re.sub('\s+', '_', self.i.seqProjectName),
                                                      'Sample_' + self.currentSample.sampleName)
                laneNumber = "%03g" % (self.currentLane.laneNumber)
                self.currentOutDirName = os.path.join(self.postAlignQCdir,
                                                      self.currentSample.sampleName) + "_L" + laneNumber

                self.lanes.append(self.currentLane)
                self.samples.append(self.currentSample)
                self.samplePaths.append(self.currentSamplePath)
                self.outDirNames.append(self.currentOutDirName)

    def setupAllSamples(
            self):  # This can be altered to be goAllSamples like in the fastQC things once we're automating things.
        for la in self.i.lanes:
            self.currentLane = la
            for sa in la.samples:
                self.currentSample = sa
                print
                "<p>", self.currentSample.sampleName, "</p>"
                self.currentJID = "0"
                self.createDBentry()

                # ToDo, set up automated version at this stage.
                # self.makeComFile()
                # self.queueFastQC()

    def makeDir(self, dirToMake):
        import os

        try:
            if not os.path.exists(dirToMake):
                os.mkdir(dirToMake)
                os.chmod(dirToMake, 0777)
        except OSError:
            print
            "<p>Failed to make directory ", dirToMake, ". </p>"

    def createDBentry(self):
        from runQuery import runQuery

        DBinsert = "insert into postAlignQC (seqProjectID,sampleID,location,sourceLocation,JID,status)"
        print
        "<p>", str(self.i.seqProjectID), "</p>"

        DBvals = " values(" + str(self.i.seqProjectID) + "," + str(
            self.currentSample.sampleID) + ",'" + self.currentOutDirName + "','" + self.currentSamplePath + "'," + self.currentJID + ",'running')"

        DBquery = DBinsert + DBvals
        print
        "<p>DBquery " + DBquery + "</p>"
        res = runQuery(DBquery)

    def makeLinksToFiles(self):
        import os

        # link all of the files in self.postAlignQCdir
        self.PAQfiles = os.listdir(self.postAlignQCdir)

        listOfLinks = []

        print
        '<div id="PAQCdetailsHidden" style="display:none">'
        print
        """<a href="javascript:swapDisplayed('PAQCdetailsShown','PAQCdetailsHidden')">hide post alignment QC files</a> """

        for t in range(0, len(self.PAQfiles)):
            loc = self.postAlignQCdir.split("/")
            # /home/biocis/demux/236_P519_MS/QCreports/postAlignQC
            loca = os.path.join("/CoreInSys/demux/", loc[4], loc[5], loc[6], self.PAQfiles[t])
            if self.PAQfiles[t][-3:] == "bam":
                print
                '<p><a href="sftp://biocis@sheba.genetics.kcl.ac.uk' + self.postAlignQCdir + '/' + self.PAQfiles[
                    t] + '">' + self.PAQfiles[t] + '</a></p>'
            else:
                print
                "<p><a href='" + loca + "'>" + self.PAQfiles[t] + "</a></p>"
        print
        '</div>'

        print
        '<div id="PAQCdetailsShown" style="display:inline">'
        print
        """<a href="javascript:swapDisplayed('PAQCdetailsHidden','PAQCdetailsShown')">show post alignment QC files</a> """
        print
        '</div>'
