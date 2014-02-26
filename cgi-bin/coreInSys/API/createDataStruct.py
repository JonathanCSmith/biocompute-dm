#!/usr/bin/python


import seqRun


i=seqRun.seqRun()

i.addSeqProject()

i.seqProjs[0].addLaneData()

i.seqProjs[0].lanes[0].addSamples()

i.seqProjs[0].lanes[0].samples[0].sampleName="testing!"

print i.seqProjs[0].lanes[0].samples[0].sampleName


#Down and up and down again

print i.seqProjs[0].lanes[0].samples[0].laneData.seqProject.seqRun.seqProjs[0].lanes[0].samples[0].sampleName

i.seqProjs[0].lanes[0].samples[0].laneData.seqProject.seqRun.seqProjs[0].lanes[0].samples[0].sampleName="Meh!"



print i.seqProjs[0].lanes[0].samples[0].sampleName


