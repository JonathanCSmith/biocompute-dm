# file getAllSeqRuns.py

# Function to retrieve all of the sequencing runs
def getAllSeqRuns():
    import runQuery

    # Perform a db lookup of the available sequencing runs
    DBquery = 'select seqRunID, startDate, completionDate, seqRunName from seqRun order by seqRunID desc'
    res = runQuery.runQuery(DBquery)
    resList = []

    # Loop through sequencing runs
    for g in range(0, len(res)):
        tmpResList = []

        # need to get all seqProjects and associted master projects.

        # DBquery="select seqProjectName, masterProjectID from seqProject where seqRunID="+str(res[g][0])
        # If it isn't linked to a master project we still need to find any seqProjects here.

        # DBquery="select seqProject.seqProjectName, masterProject.projectName, typeLinker.parentID from  seqProject, masterProject, typeLinker where seqProject.seqRunID="+str(res[g][0])+"  and  seqProject.seqProjectID=typeLinker.childID and typeLinker.parentID=masterProject.masterProjectID"
        # print "<p>",DBquery,"</p>"

        DBquery = "select seqProjectName,seqProjectID from  seqProject where seqRunID=" + str(res[g][0])
        # print "<p>",DBquery,"</p>"
        prRes = runQuery.runQuery(DBquery)
        # print "<p>prRes",prRes,"</p>"

        if len(prRes) > 0:
            for h in range(0, len(prRes)):
                DBquery = "select typeLinker.parentID, masterProject.projectName from  masterProject, typeLinker where  typeLinker.childID=" + str(
                    prRes[h][1]) + " and typeLinker.parentID=masterProject.masterProjectID"
                # print "<p>",DBquery,"</p>"
                mprRes = runQuery.runQuery(DBquery)
                # print "<p>mprRes",mprRes,"</p>"
                if len(mprRes) > 0:
                    for m in range(0, len(mprRes)):

                        prTmpResList = []
                        for i in range(0, len(res[g])):
                            prTmpResList.append(res[g][i])
                            # for j in range(0,len(prRes[h])):
                        for j in range(0, 1):
                            prTmpResList.append(prRes[h][j])
                        for jj in range(0, len(mprRes[m])):
                            prTmpResList.append(mprRes[m][jj])
                        tmpResList.append(prTmpResList)
                        # print "<p>",g,h,m,"</p>"
                else:  # got a seqProject but not a masterProject

                    prTmpResList = []
                    for k in range(0, len(res[g])):
                        prTmpResList.append(res[g][k])
                        # for l in range(0,len(prRes[h])):
                    for l in range(0, 1):
                        prTmpResList.append(prRes[h][l])
                    for m in range(0, 2):
                        prTmpResList.append("NULL")
                    tmpResList.append(prTmpResList)

        # resList.append(tmpResList)
        else:  # If we don't get any results for the run, add the results from the seqRun query and send back "NULL" for the rest
            prTmpResList = []
            for k in range(0, len(res[g])):
                prTmpResList.append(res[g][k])
            for l in range(0, 3):
                prTmpResList.append("NULL")
            tmpResList.append(prTmpResList)
            # print "<p>tmpResList",tmpResList,"</p>"
        resList.append(tmpResList)

        # print "<p>",resList,"</p>"
    return (resList)
