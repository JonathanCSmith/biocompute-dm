#!/usr/bin/python
__author__ = 'jon'

# Variable for containing our configuration properties
MySQL_DB_USERNAME = ''
MySQL_DB_PASSWORD = ''

# Function to return the username for the mysql database
def getMySQLUsername():
    if MySQL_DB_USERNAME == '':
        getConfigurationProperties()

    return MySQL_DB_USERNAME

# Function to return the password for the mysql database
def getMySQLPassword():
    if MySQL_DB_PASSWORD == '':
        getConfigurationProperties()

    return MySQL_DB_PASSWORD

# Function to load and parse the configuration file
def getConfigurationProperties():
    import os

    # Obtain the configuration file
    file = open(os.path.expanduser("~") + os.pathsep + "biocis" + os.pathsep + "ConfigurationFile.txt")
    content = file.read()
    properties = content.split("\n")

    # Parse the configuration file
    for prop in properties:
        props = prop.split(" = ")
        if props[0] == "MySQL_Username":
            global MySQL_DB_USERNAME
            MySQL_DB_USERNAME = props[1]

        elif props[0] == "MySQL_Password":
            global MySQL_DB_PASSWORD
            MySQL_DB_PASSWORD = props[1]