#!/usr/bin/python

import fastQCres

fqcr = fastQCres.fastQCres()
fqcr.fileLocation = "/var/www/biocis/link/demux/75_Sezary_P321/QCreports/FastQC/6820b_L001/6820b_ATCACG_L001_R1_001_fastqc"
print
fqcr.fileLocation

fqcr.readInQCresult()

tabs = fqcr.QCtables.keys()
for ta in range(0, len(tabs)):
    print
    tabs[ta]
    for g in range(0, len(fqcr.QCtables[tabs[ta]])):
        print
        fqcr.QCtables[tabs[ta]][g]
