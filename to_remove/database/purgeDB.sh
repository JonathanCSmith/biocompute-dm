#!/bin/bash

mysql -u msl_user -p  coreInSys < coreInSys.sql


./addLibTags.py AgilentSSBC
./addLibTags.py illuminaBC
./addLibTags.py Nextera
./addLibTags.py agilent_8bp_96_barcodes
