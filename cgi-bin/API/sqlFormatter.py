#!/usr/bin/python
__author__ = 'jon'

def stripIllegalSQLCharacters(string):
    output = string.replace("'", "")
    finalOutput = output

    return finalOutput
