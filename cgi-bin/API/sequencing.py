#!/usr/bin/python
__author__ = 'jon'

import cgitb
import json

cgitb.enable()

def get_all():
    from databasemanager import run_query
    from jsonutils import DateTimeEncoder

    query = """
        SELECT seqRun.seqRunID,
            startDate,
            completionDate,
            seqRunName,
            seqProjectName,
            projectName
        FROM seqRun
        LEFT JOIN seqProject
          ON seqRun.seqRunID = seqProject.seqRunID
        LEFT JOIN masterProject
          ON seqProject.masterProjectID = masterProject.masterProjectID
        ORDER BY seqRun.seqRunID DESC
    """

    results = run_query(query)
    results_json = json.dumps(results, cls=DateTimeEncoder)

    return results_json