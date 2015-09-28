#!/usr/bin/python
__author__ = 'jon'

import cgitb
import json

cgitb.enable()


def get_all():
    from databasemanager import run_query
    from jsonutils import DateTimeEncoder

    query = "select  masterProjectID, projectName, projectLead, status, description, openDate, lastUpdate from masterProject  order by masterProjectID desc"
    results = run_query(query)
    results_json = json.dumps(results, cls=DateTimeEncoder)

    # for result in range(0, len(results)):

    return results_json
