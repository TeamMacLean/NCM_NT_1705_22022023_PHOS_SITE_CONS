#!/usr/env python
import csv
import os
import sys
# Open the CSV file
with open(sys.argv[1], 'r') as file:
    # Create a CSV reader object
    reader = csv.reader(file)
    header = next(reader)
    # Loop over each row in the CSV file
    for row in reader:
        source = row[-2]
        dest = row[-1]
        if os.path.isfile(source) and not os.path.isfile(dest):
            os.rename(source, dest)




