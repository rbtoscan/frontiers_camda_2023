"""

This script links AMR gene information from sequencing reads to their corresponding contigs. 
It reads two input files: one containing AMR gene data and another mapping reads to contigs. 
The output maps each read to its respective gene(s) and contig(s).

Usage: python link_amr_genes_to_contigs.py <amrpp_file> <reads_contigs_file>
"""



import sys

amrpp_file=sys.argv[1]
reads_contigs=sys.argv[2]

dic_read_gene={}
dic_read_contig={}
with open(amrpp_file, 'r') as amr, open(reads_contigs, 'r') as map:
    print(amrpp_file, reads_contigs)


    # reads amr output, creating a dictionary
    # key = read, value = gene
    while True:
        l = amr.readline()
        if not l: 
            break
        else:
            read=l.split("\t")[0]
            gene= l.split("\t")[1].strip().split('|')[4]
            if read not in dic_read_gene:
                dic_read_gene[read]=[gene]
            else:
                dic_read_gene[read].append(gene)
    # read mapping output, creating dictionary
    # key = read, value = contig
    while True:
        l = map.readline()
        if not l: 
            break
        else:
            read=l.split("\t")[0]
            contig= l.split("\t")[1].strip().replace(" ","")
            if read not in dic_read_contig:
                dic_read_contig[read]=[contig]
            else:
                dic_read_contig[read].append(contig)


# iterate through amr dictionary
# if read in amr dictionary exists in mapping dictionary
    # fetch which contig, that is, mapps a read to 
dic_read_contig_amr={}
for k,v in dic_read_gene.items():
    if k in dic_read_contig:
        print (amrpp_file.split('/')[-1].split('_reads')[0]+','+';'.join(v)+',',';'.join(dic_read_contig[k]).replace(" ",""))
