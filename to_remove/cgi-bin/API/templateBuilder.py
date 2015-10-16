#!/usr/bin/python
__author__ = 'jon'

cgi = ''
header = ''
footer = ''

def writeOutCGI():
    global cgi
    if cgi == '':
        tempFile = open("../template/cgi_header", "r")
        cgi = tempFile.readlines()
        tempFile.close()

    writeOut(cgi)

    return

def writeOutHeader():
    global header
    if header == '':
        tempFile = open("../template/template_top.html", "r")
        header = tempFile.readlines()
        tempFile.close()

    writeOut(header)

    return

def writeOutFooter():
    global footer
    if footer == '':
        tempFile = open("../template/template_bottom.html", "r")
        footer = tempFile.readlines()
        tempFile.close()

    writeOut(footer)

    return

def writeOut(fileToWrite):
    for line in fileToWrite:
        print \
            line[:-1]

    return
