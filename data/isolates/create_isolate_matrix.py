

import sys
import os

isolates_file='CAMDA2023_isolates.csv'

# generated melted resistome
# cd ~/Documents/estudos/UJ/phd/projects/camda_2023/publication/0_writing/data/isolates
# while read l; do echo $l; grep "$l" CAMDA2023_isolates.csv | cut -f1 -d';' | sed "s/^/"$l"\t/g"; done < isolates_unique.txt > isolates_melted_resistome.tsv

f=open(isolates_file,'r')

isolates_resistome={}
amr_list=[]
c=1
header=f.readline()
while True:
    l=f.readline()
    if not l: break
    else:
        #print (c)
        #print (l.strip().split(';'))
        content=l.strip().split(';')
        isolate=content[0]
        bacteria=content[1]
        category=content[2]
        list_amr=content[3]
        isolates_resistome[isolate]=list_amr.strip().replace(' ','').split(',')
        for amr in list_amr.strip().replace(' ','').split(','):
            if amr in amr_list: pass
            else: amr_list.append(amr)
        c=c+1
        #print ('------')

amr_list = list(set(amr_list))

header='isolate,'
for gene in amr_list:
    header=header+gene+','
header=header[:-1]
print (header)

for id,gene_list in isolates_resistome.items():
    #print (id,gene_list)
    profile=id
    for gene in amr_list:
        if gene in gene_list:
            match=',1'
        else:
            match=',0'
        #print (gene,match)
        profile=profile+match
    print (profile)

# stdout to file to get the matrix
#python3  create_isolate_matrix.py  > isolates_matrix.csv