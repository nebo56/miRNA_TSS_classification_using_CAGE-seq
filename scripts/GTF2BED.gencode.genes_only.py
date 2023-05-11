'''
Created on Oct 15, 2017

@author: Nejc Haberman

Script will convert GTF to BED and insert transcript ID and exon number which will be useful for removing first and last exons.

'''


import sys

def convert(fname_in, fname_out):
    fin = open(fname_in, "rt")
    fout = open(fname_out, "w")
    line = fin.readline()
    while line:
        #print line
        col = line.rstrip('\n').rsplit('\t')
        chr = col[0]
        start = col[3]
        end = col[4]
        strand = col[6]
        info = col[8]
        info_col = info.rsplit(';')
        if col[2] == "gene":
            gene_id = info_col[0].replace("gene_id ","").replace('"','')
            gene_name = info_col[2].replace(" gene_name ", "").replace('"','')
            fout.write(chr + '\t' + start + '\t' + end + '\t' + gene_id + '\t' + gene_name + '\t' + strand + '\n')
        line = fin.readline()
    fout.close()
    fin.close()


if sys.argv.__len__() == 3:
    fname_in = sys.argv[1]
    fname_out = sys.argv[2]
    convert(fname_in, fname_out)
else:
    print("python GTF2BED.py <input_file> <output_file>")