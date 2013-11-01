#file demux.py       


class demux:
	def __init__(self,i,j):

		import os,stat
		import re
	
		#i is a seqProject object
		#j is a seqExpr object

		self.i=i
		self.j=j
		#print "<p>In demux object instance</p>"
		#print "<p>"+j.dataLocation+"</p>"
		self.sourceLocation="/gpfs/stagingTMP/"+j.dataLocation
		self.seqProjectID=i.seqProjectID
		self.JID="0"
		self.status="setup"
		self.projectsRoot="/home/biocis/demux"
		#print "<p>"+i.seqProjectName+"___</p>"
		self.projectDir=str(i.seqProjectID)+"_"+re.sub('\s+', '_',i.seqProjectName)
		#print "<p>"+self.projectDir+"</p>"
		dirToMake=os.path.join(self.projectsRoot,self.projectDir)			
		self.location=dirToMake
		#print"<p>location", self.location,"</p>"
		try:
			os.mkdir(dirToMake)
			#os.chmod(dirToMake,stat.S_IWRITE)
		
			#os.chmod(dirToMake,stat.S_IRWXO)
		except os.error:
			print "<p>",os.path.join(self.projectsRoot,self.projectDir),"</p>"
			
		os.chmod(dirToMake,511)
		self.params={}
	
		self.params["--use-bases-mask"]=""
		self.params["--input-dir"]="/gpfs/stagingTMP/"+j.dataLocation+"/Data/Intensities/BaseCalls/"
		self.params["--intensities-dir"]="/gpfs/stagingTMP/"+j.dataLocation+"/Data/Intensities/"
		#self.params["--positions-dir"]="/gpfs/stagingTMP/"+j.dataLocation+"/Data/Intensities/"
		self.params["--positions-format"]=".clocs"	
		self.params["--output-dir"]=os.path.join(self.projectsRoot,self.projectDir,"processed")
		self.params["--fastq-cluster-count"]="500000000"
		self.params["--flowcell-id"]=j.flowcellID
		self.params["--tiles"]=""
		self.params["--sample-sheet"]=os.path.join(self.projectsRoot,self.projectDir,"sampleSheet.csv")

		self.options={}
                self.options["--no-eamss"]="Y"
		self.options["--force"]="Y"
		self.options["--ignore-missing-stats"]="Y"
                self.options["--ignore-missing-bcl"]="Y"
                self.options["--ignore-missing-control"]="Y"


	def insertDB(self):

		#need to make sure that the same seqproject isn't already there. If so just perform an update command.


		import runQuery

                insQuery="INSERT INTO demultiplex (seqProjectID,status,sourceLocation,location,JID) "
                vals=" VALUES('"+str(self.seqProjectID)+"','"+str(self.status)+"','"+str(self.sourceLocation)+"','"+str(self.location)+"','"+str(self.JID)  +"')"

                DBins=insQuery+vals
                inser=runQuery.runQuery(DBins)


	def updateDBStatus(self,demuxID,status):
		import runQuery
		DBquery="update demultiplex set status='"+status+"' where demuxID="+str(demuxID)
		#print "<p>",DBquery,"</p>"
		res=runQuery.runQuery(DBquery)		


	def setTiles(self):
		#get the lanes from the seqExpt object
		lanes=[]
		for la in range(0,len(self.i.lanes)):
			lanes.append(self.i.lanes[la].laneNumber)
		self.params["--tiles"]="s_["
		for g in range(0,len(lanes)):
			self.params["--tiles"]=self.params["--tiles"]+str(lanes[g])+","
		self.params["--tiles"]=self.params["--tiles"][:-1]+"]"

	def setUseBasesMask(self):
		#get the length of the barcodes
		BClen=len(self.i.lanes[0].samples[0].tagSequence)
		tagMask="I"+str(BClen)
		N=int(self.j.indexTagCycles)-BClen
		for l in range(0,N):
			tagMask=tagMask+"N"

		self.params["--use-bases-mask"]="Y"+str(self.j.readCycles)+","+tagMask+","+"Y"+str(self.j.readCycles)

	def makeComFile(self,projectName,projectID):
		#write the command file
		commandFile=""
		parList=self.params.keys()
		print '<form method="post" action="writeDemuxCommand?projName='+projectName+'&projID='+projectID+'">'
		print '<table border=0 >'
		for co in range(0,len(parList)):
			print "<tr><td>"+parList[co]+'</td><td><input type="text" name="'+parList[co]+'" value="'+self.params[parList[co]]+'" size="80"></td></tr>'
		optList=self.options.keys()
		for opts in range(0,len(optList)):
			print "<tr><td>"+optList[opts]+'</td><td><select name="'+optList[opts]+'" value="'+self.options[optList[opts]]+'" >'
			print '<option value="Y">Y</option>'
			print '<option value="N">N</option>'
			print "</select></td></tr>"

		print "</table>"
		print '<input type="submit" value="Submit">'
		print '</form>'

		for co in range(0,len(parList)-1):
			commandFile=commandFile+parList[co]+"\t"+self.params[parList[co]]+" \\\n<br>"
		commandFile=commandFile+parList[-1]+"\t"+self.params[parList[-1]]+" \\\n<br>"
		#print "<p>",commandFile,"</p>"

	def queueDemuxJob(self):
		#print "<p>In the queue method!</p>"
		import os
		import subprocess as sub
		import re
		#print "<p>self.projectDir"+self.projectDir+"</p>"		
		#casavaComFile=os.path.join(self.projectDir,"casava_Bcl2FastQ.sh")
		command='ssh biocis@apollo "cd /home/biocis/demux/'+self.projectDir+'; qsub -q shortterm.q,longterm.q casava_Bcl2FastQ.sh"'

        	#print "<p>"+command+"</p>"
						
        	p = sub.Popen(command, shell=True, stdout=sub.PIPE,stdin=sub.PIPE).stdout
        	line=p.readlines()
		print "<p>",line[0],"</p>"
		JIDfile='/home/biocis/demux/'+self.projectDir+'/JID'
		jid="0"
		if len(line[0])>0:
			bits=line[0].split()
			jid=bits[2]
		print "Job ID="+jid
		JIDF=open(JIDfile,"w")
		#This needs re-doing
		if re.match("[0-9]",jid):
			if re.match("[A-Z]|[a-z]",jid):
				JIDF.write("0")
			else:
				JIDF.write(jid)
				self.JID=jid
		else:
			JIDF.write("0")

		JIDF.close()


	def getDemuxStatus(self):
		import os
		import subprocess as sub
		#print "/home/biocis/demux/"+self.projectDir	
		stuff=os.listdir("/home/biocis/demux/"+self.projectDir)
		
		#print "<p>",stuff,"</p>"
		DBquery="select demuxID, status, location, sourceLocation, JID from demultiplex where seqProjectID='"+str(self.seqProjectID)+"'"
		#print "<p>"+DBquery+"</p>"

		import runQuery

		res=runQuery.runQuery(DBquery)

		#print res
		
		command='ssh biocis@apollo "qstat"'
                p = sub.Popen(command, shell=True, stdout=sub.PIPE,stdin=sub.PIPE).stdout
                qstatLine=p.readlines()	
	
		
		print "<table border='1'>"
		print "<tr><th>demuxID</th><th>status</th><th>location</th><th>source location</th><th>JID</th></tr>"
		for arr in range(0,len(res)):
			runStat=""
			running="N"
			JID=str(res[arr][4])
			for lineOut in range(0,len(qstatLine)):
                                m=qstatLine[lineOut].split()
                                #print "<p>",m,"</p>"
                                try:
                                        jobID=int(m[0])
                                except Exception:
                                        jobID=0 
                                if jobID==int(JID):
                                        runStat=qstatLine[lineOut]
					running="Y"
			
			
			demuxID=res[arr][0]
			status=res[arr][1]
			location=res[arr][2]
			sourceLocation=res[arr][3]
			#self.updateDBStatus(demuxID,status)
			if running=="Y":
				status="running"
			else:
				status="finished"

			print "<tr><td>"+str(demuxID)+"</td><td>"
			if status=="finished":

				pa=location.split("/")
				#print pa
				projDir=pa[-1]
				print "<a href='/CoreInSys/demux/"+projDir+"/processed/Basecall_Stats_"+self.j.flowcellID+"/Demultiplex_Stats.htm'>"+status+"</a>"
			else:
				print status

			print "</td><td>"+location+"</td><td>"+sourceLocation+"</td><td>"+str(JID)+"</td><td>",runStat,"</td>"
			
			print "</tr>"
			self.updateDBStatus(demuxID,status)
	
		print "</table>"	
		
		


	def writeCommandFile(self):
		import os
		import re
		fileName=os.path.join(self.projectsRoot,self.projectDir,"casava_Bcl2FastQ.sh")
		OF=open(fileName,"w")
		header="#!/bin/bash\n#$ -cwd\n#$ -j y\n#$ -S /bin/bash\n#$ -pe threaded 8\n\n"
		OF.write(header)
		commands="/apps/casava/1.8.2/bin/configureBclToFastq.pl \\\n"
		parList=self.params.keys()
		for co in range(0,len(parList)):
			paramToWrite=re.sub('\s+', '_', self.params[parList[co]])
                        commands=commands+"\t\t"+parList[co]+"\t"+paramToWrite+" \\\n"
                #commands=commands+"\t\t"+parList[-1]+"\t"+self.params[parList[-1]]+" \\\n"
		#OF.write(commands)
		optsText=""
		optList=self.options.keys()
		#numbOptions=0
		for opts in range(0,len(optList)):
			if self.options[optList[opts]]=="Y":
				optsText=optsText+"\t\t"+optList[opts]+" \\\n"
				#numbOptions=numbOptions+1
		if len(optsText)>0:
			optsText=optsText[:-2]
		else:
			commands=commands[:-2]

		OF.write(commands)
		OF.write(optsText)
		makeComm="\n\nmake -j 8 -C "+self.params["--output-dir"]+"\n"  
		OF.write(makeComm)	

		OF.close()



	def makeSampleSheet(self):
		import re

		sampleSheet="FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject\n"
		#go through each lane
		FCID=self.j.flowcellID
		Control="N"
		SampleRef=""
		Description=""
		Recipe=""
		Operator=""
		SampleProject=self.i.seqProjectName
		#go through each lane
		for la in range(0,len(self.i.lanes)):
			#Now each sample
			Lane=str(self.i.lanes[la].laneNumber)
			for sa in range(0,len(self.i.lanes[la].samples)):	
				SampleID=self.i.lanes[la].samples[sa].sampleName	
				Index=self.i.lanes[la].samples[sa].tagSequence
				sampleSheet=sampleSheet+re.sub('\s+', '_',FCID+","+Lane+","+SampleID+","+SampleRef+","+Index+","+Description+","+Control+","+Recipe+","+Operator+","+SampleProject)+"\n"
		#sampleSheet=re.sub('\s+', '_',sampleSheet)
		#print "<p>"+sampleSheet+"</p>"
		SF=open(self.params["--sample-sheet"],"w")
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





