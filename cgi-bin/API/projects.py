#!/usr/bin/python
__author__ = 'jon'

import cgitb
import json

cgitb.enable()

def get_all():
    from databasemanager import run_query
    from jsonutils import DateTimeEncoder

    query = """
        SELECT masterProjectID,
          projectName,
          projectLead,
          status,
          description,
          openDate,
          lastUpdate
        FROM masterProject
          ORDER BY masterProjectID DESC
    """

    results = run_query(query)
    results_json = json.dumps(dict(results), cls=DateTimeEncoder)

    return results_json
