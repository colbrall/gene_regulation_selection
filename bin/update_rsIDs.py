#!/usr/bin/env python
# update_rsIDs.py
# Laura Colbran, 5-25-17
# updates rsIDs in a dosage file to the ones in PrediXcan model by genomic location
# also filters by whether or not a location has an rsID
#
# Python 3
#
#USAGE:
# python update_rsIDs.py FLAG RS_COL PATH/TO/PREDIXCAN/SNPS PATH/TO/FILE/*
# FLAG can be bed or dosage, depending on format of File you want to convert
# RS_COL is an integer, and is the column in the Predixcan file you want to get the rsIDs from. 0-indexed.

import sys
import string
import gzip

def dosageConvert(file_list,rs_col,ref_file):
    for fil in file_list:
        dict = {} #{loc:rsID}
        with gzip.open(fil, 'rt') as f:
            #print(f.readlines(1))
            chr = f.readlines(1)[0].split('\t')[0].split('r')[-1]
            if chr == "": chr = f.readlines(1)[0].split('\t')[0].split('r')[-1]
            f.seek(0)
            with gzip.open(ref_file,'rt') as g: #dict with PrediXcan rsIDs
                for l in g:
                    if l.startswith(chr + "\t"):
                        line = l.rstrip('\n')
                        dict[line.split("\t")[1]] = (line.split("\t")[int(rs_col)],line.split("\t")[3],line.split("\t")[4])
#           output
            name = "chr" + chr + ".rs_updated.dos"
            of = open(name,"w")
            for line in f:
                if line.startswith("#"):
                    of.write(line)
                    continue
                l = line.split("\t")
                if l[2] in dict:
                    if l[3] == dict[l[2]][1] or l[3] == dict[l[2]][2]: #check whether predixcan ref allele is present
                        if l[4] == dict[l[2]][1] or l[4] == dict[l[2]][2] or l[4] == ".": #check whether predixcan alt allele is present
                            l[1]= dict[l[2]][0] #plug in predixcan ids where applicable and write line out
                            of.write("\t".join(l))
            of.close()

def bedConvert(file_list,rs_col,ref_file):
    dict = {} #{loc:rsID}
    with gzip.open(ref_file,'r') as g: #dict with PrediXcan rsIDs
        for line in g:
            if line.startswith('#'): continue
            dict['\t'.join([line.split('\t')[0],line.split("\t")[1]])] = (line.split("\t")[rs_col])
    for fil in file_list:
        with open(fil, 'r') as inf:
            with open("%s_rsupdated.bed" % (".".join(fil.split(".")[0:-1])),'w') as out:
                out.write("#chr\tpos-1\tpos\trsID\n")
                for line in inf:
                    if line.startswith('#'): continue
                    pos = line.split("\t")[0:3]
                    pos[0] = pos[0].split("r")[1]
                    try:
                        rs = dict["\t".join([pos[0],pos[2]])]
                    except:
                        rs = "."
                    out.write("chr%s\t%s\t%s\t%s\n" % (pos[0],pos[1],pos[2],rs))

def main():
    if sys.argv[1] == "dosage":
        dosageConvert(sys.argv[4:],sys.argv[2],sys.argv[3])
    elif sys.argv[1] == "bed":
        bedConvert(sys.argv[4:],sys.argv[2],sys.argv[3])
    else:
        print("ERROR: Specify format as bed or dosage!")

if __name__ == '__main__': main()
