#!/usr/bin/python

import sys
# tagFile="illuminaBC"

tagFile = sys.argv[1]

BCF = open(tagFile, "r")
lines = BCF.readlines()

import MySQLdb as msd

DBname = "coreInSys"
con = msd.connect('localhost', 'msl_user', 'msl_pass', DBname);

for g in range(0, len(lines)):
    li = lines[g][:-1].split(",")

    DBquery = "insert into IDtagLibs (libraryName,libraryTagID, tagSequence) values( '" + li[0] + "','" + li[
        1] + "','" + li[2] + "')"
    print
    DBquery
    cur = con.cursor()

    cur.execute(DBquery)

    con.commit()
    res = cur.fetchall()
    cur.close()
