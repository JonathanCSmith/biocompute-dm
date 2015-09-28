#!/usr/bin/python
__author__ = 'jon'

import json

class DateTimeEncoder(json.JSONEncoder):
    def default(self, obj):
        return str(obj)
