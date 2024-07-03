"""
This script filters an AMR analytic matrix to retain only rows corresponding to specific isolates. 
It uses a list of isolate IDs from a text file to filter the input matrix, keeping only the relevant AMR data. 
The output is a filtered version of the input table with isolates and AMRs.

Usage: python filter_isolates_amr_data.py
"""

import sys
import os

import pandas as pd


isolates_file='isolates_unique.txt'
input_table='AMR_analytic_matrix.csv'

ids={}
# Read the first file with single column IDs
with open(isolates_file, 'r') as f:
    #ids = set(line.lower().strip() for line in f)
    for line in f:
        if line.lower().strip() in ids: pass
        else: 
            ids[line.lower().strip()]=line.strip()

#print (ids)
#print (ids)
# Read the second file with three columns
print("amrpp;isolates")
filtered_rows = []
with open(input_table, 'r') as f:
    for line in f:
        if '|' in  line:
            #print (line)
            columns = line.strip().split(',')
            #print (columns)
            amr=columns[0].split('|')[4]
            #print (amr.lower())
            if amr.lower() in ids:
                print (columns[0]+";"+ids[amr.lower()])
                #print ('-------')
                #print (line.strip())
                filtered_rows.append(line)
        else:
            #print (line.strip())
            filtered_rows.append(line)
#print (filtered_rows)
with open(input_table.replace(".csv","_isolates-AMRs.csv"), "w") as f:
    f.writelines(filtered_rows)


