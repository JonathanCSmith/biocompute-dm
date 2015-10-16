#!/usr/bin/python
import os

files = os.listdir("../template")

pages = []
TF = open("../template/template_top.html", "r")
template_top = TF.readlines()
TF.close()

BF = open("../template/template_bottom.html", "r")
template_bottom = BF.readlines()
BF.close()

# Handle null dir
if not os.path.exists("../html"):
    os.makedirs("../html")

# Calculate the pages that we are building
for m in range(0, len(files)):
    if files[m][-5:] == ".html" and files[m] != "template_top.html" and files[m] != "template_bottom.html":
        pages.append(files[m])

for wr in pages:
    # Load the page specific content
    CT = open("../template/" + wr, "r")
    content = CT.readlines()
    CT.close()

    # Name the page
    pageName = "../html/" + wr
    print \
        pageName
    page = []

    # Build the page
    for top in template_top:
        page.append(top)
    for co in content:
        page.append(co)
    for bottom in template_bottom:
        page.append(bottom)

    # Write the page
    PA = open(pageName, "w")
    for p in page:
        PA.write(p)
    PA.close()
