# file runQuery.py

# Generic function for running a query
def runQuery(DBquery):
    import MySQLdb as msd
    import configurationHandler as config

    # Obtain the connection properties
    user = config.getMySQLUsername()
    password = config.getMySQLPassword()

    # Obtain a connection to the database
    DBname = "coreInSys"
    con = msd.connect('localhost', user, password, DBname);
    cur = con.cursor()

    # Execute the query
    cur.execute(DBquery)

    # Commit the query and return the result
    con.commit()
    res = cur.fetchall()

    # Close the query
    cur.close()

    return (res)
