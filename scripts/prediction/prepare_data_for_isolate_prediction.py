"""

This script filters a transformed matrix to retain only rows corresponding to specific isolates listed in a separate file. 
It reads isolate IDs from a text file and uses these to filter the input matrix, preparing it for isolate prediction tasks.
The output is a filtered matrix containing only the relevant isolates and AMRs.

Usage: python prepare_data_for_isolate_prediction.py
"""

import sys
import os

isolates_file='isolates_unique.txt'
input_table='transformed_matrix_for_norm_apr15_50p_all_mge_amr.tsv'

# Read the first file with single column IDs
with open(isolates_file, 'r') as f:
    ids = set(line.lower().strip() for line in f)


# Read the second file with three columns
filtered_rows = []
with open(input_table, 'r') as f:
    for line in f:
        #print (line)
        columns = line.strip().split('\t')
        amr=columns[1]
        if amr.lower() in ids:
            print (line.strip())
            filtered_rows.append(line)
#print (filtered_rows)

with open(input_table.replace(".tsv","-isolates-AMRs.tsv"), "w") as f:
    f.writelines(filtered_rows)
