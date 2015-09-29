#!/usr/bin/python
__author__ = 'jon, alan'

import cgitb
cgitb.enable()

def run_query(query):
    import MySQLdb as msd
    import configurationhandler as config
    from sqlFormatter import stripExcess

    # Obtain the connection properties
    user = config.getMySQLUsername()
    password = config.getMySQLPassword()

    # Obtain a connection to the database
    DBname = "coreInSys"
    con = msd.connect('localhost', user, password, DBname);
    cur = con.cursor()

    # Execute the query
    clean_query = stripExcess(query)
    cur.execute(clean_query)

    # Commit the query and return the result
    con.commit()
    res = cur.fetchall()

    # Close the query
    cur.close()

    # Build a JSON object
    columns = [desc[0] for desc in cur.description]
    results = []
    for row in res:
        row = dict(zip(columns, row))
        results.append(row)

    return results
