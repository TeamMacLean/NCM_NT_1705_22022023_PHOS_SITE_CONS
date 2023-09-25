#!/usr/bin/bash

#source /nbi/software/testing/bin/fasta3-36
source blast+-2.9.0
source python-3.8.3


python scripts/find_site.py ${1} ${2} ${3} ${4} ${5} > ${6}