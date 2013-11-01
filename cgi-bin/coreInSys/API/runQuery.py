#file runQuery.py


def runQuery(DBquery):

	#print DBquery

        import MySQLdb as msd
        #print "<p>",DBquery,"</p>"
	DBname="coreInSys"
        con = msd.connect('localhost', 'msl_user', 'msl_pass', DBname);
        cur = con.cursor()

        cur.execute(DBquery)

        con.commit()
        res=cur.fetchall()
        cur.close()

        return(res)




