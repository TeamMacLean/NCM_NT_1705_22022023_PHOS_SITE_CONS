#!/usr/bin/env python
'''
Gets the raw proteome fasta locally from the relevant databases
as described in the source column of the proteomes file

Uses source python-3.8.3 on HPC for a python with pandas
'''


import os
import sys
import pandas as pd
import subprocess

target = sys.argv[1]

if not os.path.isdir(target):
    print("directory {} does not exist. Please create it and try again".format(target))
    sys.exit()

df = pd.read_csv("lib/proteomes.csv")
for url in df['source']:
    cmd = "wget -P {} {}".format(target, url)
    print("getting {}".format(url))
    subprocess.call(cmd, shell=True)
