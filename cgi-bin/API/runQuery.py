# file runQuery.py

# Generic function for running a query
def runQuery(DBquery):
    import MySQLdb as msd

    # Obtain the connection properties

    # Obtain a connection to the database
    DBname = "coreInSys"
    con = msd.connect('localhost', 'msl_user', 'msl_pass', DBname);
    cur = con.cursor()

    # Execute the query
    cur.execute(DBquery)

    # Commit the query and return the result
    con.commit()
    res = cur.fetchall()

    # Close the query
    cur.close()

    return (res)
