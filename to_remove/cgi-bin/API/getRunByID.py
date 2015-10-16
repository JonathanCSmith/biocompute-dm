# file getRunByID.py



def getRunByID(exptID):
    import seqRun

    # Create a seqRun instance
    i = seqRun.seqRun()
    # return(i)

    i.getSeqNameByID(exptID)
    # print "<p>i.seqRunName",i.seqRunName,"</p>"
    # return(i)



    # i.addSeqProject()
    # i.seqProjs[-1].addLaneData()
    # i.seqProjs[-1].lanes[-1].addSamples()

    # Get the IDs of all the projects that belong to this run
    seqProjectIDs = i.getSeqProjects()

    i.addSeqProject()
    # print "<p>seqProjectIDs",seqProjectIDs,"</p>"
    # Work down to the Projects layer
    for pr in range(0, len(seqProjectIDs)):
        # add a project instance
        if pr > 0 or len(i.seqProjs) == 0:
            if i.seqProjs[-1].seqProjectID != 0:
                i.addSeqProject()
        i.seqProjs[-1].getSeqProjByID(seqProjectIDs[pr][0])
        # print "<p>i.seqProjs[-1].getSeqProjByID(seqProjectIDs[pr][0])",i.seqProjs[-1].getSeqProjByID(seqProjectIDs[pr][0]),"</p>"

        # Get the IDs of all of the lanes in this project
        laneDataIDs = i.seqProjs[-1].getLaneDataIDs()
        i.seqProjs[-1].addLaneData()


        # Work down to the lane layer
        for la in range(0, len(laneDataIDs)):
            # add a lane instance
            if la > 0 or len(i.seqProjs[-1].lanes) == 0:
                if i.seqProjs[-1].lanes[-1].laneID != "NULL":
                    i.seqProjs[-1].addLaneData()
            i.seqProjs[-1].lanes[-1].getLaneDataByID(laneDataIDs[la][0])

            sampleIDs = i.seqProjs[-1].lanes[-1].getSampleIDs()
            # print sampleIDs
            i.seqProjs[-1].lanes[-1].addSamples()

            # Now we move into the sample layer
            for sa in range(0, len(sampleIDs)):
                # add a sample instance
                if sa > 0 or len(i.seqProjs[-1].lanes[-1].samples) == 0:
                    if i.seqProjs[-1].lanes[-1].samples[-1].sampleID != 0:
                        i.seqProjs[-1].lanes[-1].addSamples()
                i.seqProjs[-1].lanes[-1].samples[-1].getSampleDataByID(sampleIDs[sa][0])

    return (i)
