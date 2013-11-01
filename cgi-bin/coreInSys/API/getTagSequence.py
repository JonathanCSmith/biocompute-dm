#!/usr/bin/python

import runQuery

def getTagSequence(kitID,tagNumber):
	DBquery="select tagSequence from IDtagLibs where libraryName='"+kitID+"' and libraryTagID='"+tagNumber+"'"
	#print "<p>",DBquery,"</p>"
	res=runQuery.runQuery(DBquery)
	#print "<p>",res[0][0],"</p>"
	return(res[0][0])
