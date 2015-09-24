# biocis
Repository for the BioCIS pipeline

To Setup:
Setup a machine with a Lamp(y) stack
Clone branch onto machine for deployment
Move all files (except db schema) into apache2 website location
Build using edit.py in the build folder
Link apache2 website root directory to /biocis/
Link apache2 directory root (for the website) to /biocis/html/index.html
Deploy mysql schema into db
Change configuration files to reflect db propertes
Profit?
