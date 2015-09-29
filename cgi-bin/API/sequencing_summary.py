#!/usr/bin/python
__author__ = 'jon'

import cgitb
import json

cgitb.enable()

# Function to return all sequencing runs de-normalized by sequencing project and master project
def get_all():
    from databasemanager import run_query
    from jsonutils import DateTimeEncoder

    # Query to obtain all sequence runs
    query = """
        SELECT seqRun.seqRunID,
            startDate,
            completionDate,
            seqRunName,
            seqProjectName,
            masterProject.masterProjectID,
            projectName
        FROM seqRun
        LEFT JOIN seqProject
          ON seqRun.seqRunID = seqProject.seqRunID
        LEFT JOIN masterProject
          ON seqProject.masterProjectID = masterProject.masterProjectID
        ORDER BY seqRun.seqRunID DESC
    """

    # Obtain the results as a json and hand off
    results = run_query(query)
    results_json = json.dumps(results, cls=DateTimeEncoder)
    return results_json
